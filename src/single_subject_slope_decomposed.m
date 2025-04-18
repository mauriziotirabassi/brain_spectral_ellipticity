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
logd     = log10(dvals);
IQR      = prctile(logd,75) - prctile(logd,25);
bw       = 2 * IQR / numel(logd)^(1/3);
nbins    = ceil( (max(logd)-min(logd)) / bw );
edges    = linspace(min(dvals),max(dvals),nbins+1);
centers  = (edges(1:end-1)+edges(2:end))/2;

%% LOOP OVER SUBJECTS
max_lag  = 45;          % as decided
nLags    = 100;
tau_vals = linspace(0,max_lag,nLags);

%- select subject number
iSub = 1;

%— load this subject’s fitted A, TR, noise variance
S = load(fullfile(outDir,files{iSub}));
A = S.A;
Q = eye(n)*S.output.eff_conn.NoiseVar;
P = lyap(A,Q);
P_tau = @(t) expm(A*t)*P;

% symmetric and skew-symmetric generator decomposition
A_sym   = -0.5 * Q / P;
A_skew  = (A * P) - A_sym * P;
A_skew  = A_skew / P;

P_tau_sym  = @(t) expm(A_sym * t)  * P;
P_tau_skew = @(t) expm(A_skew * t) * P;

%— stack all lagged covariances
Ctens       = nan(n,n,nLags);
Ctens_sym   = nan(n,n,nLags);
Ctens_skew  = nan(n,n,nLags);
for k=1:nLags
    Ctens(:,:,k)       = P_tau(tau_vals(k));
    Ctens_sym(:,:,k)   = P_tau_sym(tau_vals(k));
    Ctens_skew(:,:,k)  = P_tau_skew(tau_vals(k));
end

%— now for each tau, bin the off‐diagonal P_ij(tau) against dvals
a_vec       = nan(1,nLags);
a_sym_vec   = nan(1,nLags);
a_skew_vec  = nan(1,nLags);

for k=1:nLags
    % FULL dynamics
    cvals = Ctens(:,:,k);
    c_off  = cvals(mask);
    B = nan(1,nbins);
    for b=1:nbins
        sel  = dvals>=edges(b) & dvals<edges(b+1);
        B(b) = mean(c_off(sel));
    end
    fitR = [8.13 33.82];
    idx = centers>=fitR(1) & centers<=fitR(2);
    p = polyfit(log10(centers(idx)), log10(abs(B(idx))),1);
    a_vec(k) = p(1);

    % SYMMETRIC dynamics
    cvals_sym = Ctens_sym(:,:,k);
    c_off_sym = cvals_sym(mask);
    B_sym = nan(1,nbins);
    for b=1:nbins
        sel = dvals>=edges(b) & dvals<edges(b+1);
        B_sym(b) = mean(c_off_sym(sel));
    end
    p_sym = polyfit(log10(centers(idx)), log10(abs(B_sym(idx))),1);
    a_sym_vec(k) = p_sym(1);

    % SKEW-SYMMETRIC dynamics
    cvals_skew = Ctens_skew(:,:,k);
    c_off_skew = cvals_skew(mask);
    B_skew = nan(1,nbins);
    for b=1:nbins
        sel = dvals>=edges(b) & dvals<edges(b+1);
        B_skew(b) = mean(c_off_skew(sel));
    end
    p_skew = polyfit(log10(centers(idx)), log10(abs(B_skew(idx))),1);
    a_skew_vec(k) = p_skew(1);
end

% optional: plot the three together
figure;
plot(tau_vals, a_vec, 'k'); hold on
plot(tau_vals, a_sym_vec, 'b--');
plot(tau_vals, a_skew_vec, 'r--');
legend('full A', 'symmetric part', 'skew-symmetric part');
xlabel('\tau'); ylabel('Slope a(\tau)');
title('Distance-decay exponent vs lag');
grid on;