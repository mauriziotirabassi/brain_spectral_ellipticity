% clear workspace and set paths
clearvars; close all; clc

dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% list of subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

% unique networks & their indices
rname = fullfile(dataDir, 'regions.xlsx');
T     = readtable(rname);
netLabels = T.NETWORK;
[uniqueNets, ~, ic] = unique(netLabels);

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';         % make symmetric
n        = size(Dfull,1);
fitR     = [4.48, 12.18]; % select distance range for fitting

% define signed lags
max_lag = 50;
nLags   = 200;
tau_vals= linspace(0, max_lag, nLags);

% load model for selected subject
iSub    = 1;
subj    = load(fullfile(outDir,files{iSub}));
A       = subj.A;
Q       = eye(n) * subj.output.eff_conn.NoiseVar;
P       = lyap(A, Q);
P_tau   = @(t) expm(A*t)*P;
stds    = sqrt(diag(P));     % n×1 vector of standard deviations
normMat = stds * stds.';     % n×n matrix where (i,j)=σ_i*σ_j

% calculate lagged-covariance
Ctens   = nan(n,n,nLags); % preallocate
for k = 1:nLags
    t = tau_vals(k);
    C    = P_tau(t) ./ normMat;
    C    = max(C, C.');
    Ctens(:,:,k) = abs(C);
end

%% SINGLE NETWORK TRIAL (row‐only)
netId   = 6;
netName = uniqueNets{netId};
idx     = find(ic == netId);   % rows for this network
m       = numel(idx);          % number of nodes in this net

% extract the m×n distance block (rows only)
Drows   = Dfull(idx, :);       % size m×n

% build a mask to exclude self‐distances (row i to column idx(i))
mask    = true(m, n);
for i = 1:m
    mask(i, idx(i)) = false;
end

% vectorize distances for fitting
dvals   = Drows(mask);
sel     = dvals >= fitR(1) & dvals <= fitR(2);

% preallocate slope array
alpha   = nan(1, nLags);

% loop over lags: pull row‐only covariances & regress
for t = 1:nLags
    Cfull = Ctens(:,:,t);        % 74×74 at lag tau_vals(t)
    Crows = Cfull(idx, :);       % m×n block of covariances
    cvals = abs( Crows(mask) );  % vectorized

    p     = polyfit(log(dvals(sel)), log(cvals(sel)), 1);
    alpha(t) = p(1);
end

% plot α(τ)
figure;
plot(tau_vals, alpha, 'LineWidth',1);
title(netName, 'Interpreter','none');
xlabel('\tau');
ylabel('\alpha(\tau)');
grid on;

% zero-lag covariance vs distance for this network
C0      = Ctens(:,:,1);          % tau = 0
C0rows  = C0(idx, :);            % m×n block
y_all   = abs( C0rows(mask) );

% sort & plot
[d_sorted, sidx] = sort(dvals);
y_sorted        = y_all(sidx);

figure;
plot(log(d_sorted), log(y_sorted));
hold on;
xline(log(fitR(1)),'r--');
xline(log(fitR(2)),'r--');
xlabel('log(distance)');
ylabel('log(|covariance|)');
title(['Zero-lag Covariance vs Distance — ' netName], 'Interpreter','none');
grid on;

%% FITTING FOR EACH NETWORK (rows only)
nNets = numel(uniqueNets);
nCols = ceil(sqrt(nNets));
nRows = ceil(nNets/nCols);

figure('Units','normalized','Position',[.1 .1 .8 .8]);
tiledlayout(nRows, nCols, 'Padding','compact','TileSpacing','compact');

for k = 1:nNets
    netName = uniqueNets{k};
    idx     = find(ic == k);        % row‐indices for this network
    m       = numel(idx);           % number of nodes in this net

    % extract the m×n distance block
    Drows   = Dfull(idx, :);        % size m×n
    % we’ll ignore self‐distances (diagonal) when idx==col
    mask    = true(m, n);
    for i = 1:m
        mask(i, idx(i)) = false;
    end
    dvals   = Drows(mask);          % vector of all row‐wise distances

    % preallocate slope array
    alpha   = nan(1, nLags);

    for t = 1:nLags
        % extract the corresponding covariance block
        Crows  = Ctens(idx, :, t);   % m×n
        cvals  = abs( Crows(mask) ); % vector, same length as dvals

        % threshold to fit only within [fitR(1) , fitR(2)]
        sel    = (dvals >= fitR(1) & dvals <= fitR(2));

        % linear fit: log(c) ~ α * log(d)
        p      = polyfit(log(dvals(sel)), log(cvals(sel)), 1);
        alpha(t) = p(1);
    end

    % plot
    nexttile;
    plot(tau_vals, alpha);
    title(netName, 'Interpreter','none');
    xlabel('\tau');
    ylabel('\alpha(\tau)');
    grid on;
end