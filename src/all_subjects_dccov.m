%% ------------------------------------------------------------------------
%  SETTINGS & DATA‐FOLDER LAYOUT
%  [1] data/
%       ├─ distances.mat        % contains mtx_euc_dis
%       ├─ structural.mat       % contains full_connectome_no_symm
%       └─ model_outputs/       % one .mat per subject, each with A_hat, A, TR, output.eff_conn.NoiseVar, etc.
%           ├─ sub001.mat
%           ├─ sub002.mat
%           └─ …
%  ------------------------------------------------------------------------

clearvars, close all; clc

dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% get list of your subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

%% PRELOAD EVERYTHING THAT NEVER CHANGES
% Euclidean‐distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';      % make symmetric
n        = size(Dfull,1);
mask     = ~eye(n);
dvals    = Dfull(mask);

% build your Freedman–Diaconis BINS just once
logd     = log(dvals);
IQR      = prctile(logd,75) - prctile(logd,25);
bw       = 2 * IQR / numel(logd)^(1/3);
nbins    = ceil( (max(logd)-min(logd)) / bw );
edges    = linspace(min(dvals),max(dvals),nbins+1);
centers  = (edges(1:end-1)+edges(2:end))/2;

%% LOOP OVER SUBJECTS
max_lag  = 100;          % as decided
nLags    = 200;
tau_vals = linspace(0,max_lag,nLags);

alpha_mat = nan(numel(files),nLags);

for iSub = 1:numel(files)
    %— load this subject’s fitted A, TR, noise variance
    subj = load(fullfile(outDir,files{iSub}));
    A = subj.A;
    Q = eye(n)*subj.output.eff_conn.NoiseVar;
    P = lyap(A,Q);
    S = 0.5 * (A * P - P * (A.'));
    A_skew  = S / P;
    P_tau_skew = @(t) expm(A_skew * t) * P;
    
    %— stack all lagged covariances
    Ctens = nan(n,n,nLags);
    for k=1:nLags
        Ctens(:,:,k) = P_tau_skew(tau_vals(k));
    end

    %— now for each tau, bin the off‐diagonal P_ij(tau) against dvals
    a_vec = nan(1,nLags);
    for k=1:nLags
        cvals = Ctens(:,:,k);
        c_off  = cvals(mask);

        % bin‐average:
        B = nan(1,nbins);
        for b=1:nbins
            sel    = dvals>=edges(b) & dvals<edges(b+1);
            B(b)   = mean(c_off(sel));
        end

        % TODO: now only absolute value!
        % fit only in your desired r‐range:
        fitR   = [4.48, 12.18];
        idx    = centers>=fitR(1) & centers<=fitR(2);
        p      = polyfit(log(centers(idx)), log(abs(B(idx))),1);
        a_vec(k) = p(1);
    end
    
    alpha_mat(iSub,:) = a_vec;
end

%% PLOT: One tile per subject
nSubs = size(alpha_mat,1);
% choose a roughly square grid:
nCols = ceil(sqrt(nSubs));
nRows = ceil(nSubs/nCols);

figure('Renderer','painters','Position',[100 100 1200 800]);
t = tiledlayout(nRows,nCols, ...
      'TileSpacing','compact','Padding','compact');

for iSub = 1:nSubs
    ax = nexttile;
    plot(ax, tau_vals, alpha_mat(iSub,:));
    xlim(ax, [min(tau_vals) max(tau_vals)]);
    ylim(ax, [min(alpha_mat(:)) max(alpha_mat(:))]);
    xlabel(ax, '\tau');
    ylabel(ax, '\alpha');
    title(ax, sprintf('Subj %d', iSub));
end

% optionally, give the entire figure a super‐title:
sgtitle(t,'Spatial‐Decay Exponent \alpha(\tau) per Subject');
