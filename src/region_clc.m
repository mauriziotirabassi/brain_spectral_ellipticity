clear; clc; %close all
rng(42);

% Data
dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat'); % distance matrix
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat'); % structural conn matrix
outDir = fullfile(dataDir,'regressed_001_01_sim62131');
d = dir(fullfile(outDir,'*.mat')); % list of subject‐model .mat files
files = {d.name};

% Subject data
iSub = 5; % subject number
subj = load(fullfile(outDir,files{iSub}));
A = subj.A; % effective connectivity
n = size(A, 1); I = eye(n);
Sigma_w = I * subj.output.eff_conn.NoiseVar; % noise covariance
Sigma = lyap(A, Sigma_w); % zero-lag covariance

% Information propagator
[U, D] = eig(Sigma);
A_pc = U' * A * U;
D_half     = sqrt(D);
D_half_inv = diag(1 ./ sqrt(diag(D)));
A_tilde_pc = D_half_inv * A_pc * D_half;

% % DYNAMIC ANISOTROPY INDEX (DAI) SPECTRUM
[~, T] = schur(A_tilde_pc, 'real');
kappa_spec = [];
i = 1;
while i <= n
    if i < n && abs(T(i + 1, i)) ~= 0 % 2x2 oscillatory Schur block
        b = T(i, i + 1); c = T(i + 1, i);
        kappa = (b^2 + c^2) / (2 * abs(b * c));
        kappa_spec(end + 1) = kappa; %#ok<SAGROW>
        i = i + 2;
    else, i = i + 1;
    end
end
figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, histogram(kappa_spec, 10); axis square;
title('$\kappa$-Spectrum', 'Interpreter', 'latex', 'FontSize', 15)

% SUBNETWORK SELECTION
rname = fullfile(dataDir, 'regions.xlsx');
Tr = readtable(rname);
netLabels = Tr.NETWORK;
[uniqueNets, ~, ic] = unique(netLabels);

netId = 1;
netName = uniqueNets{netId};
idx = find(ic == netId); % rows for this network
m = numel(idx); % number of nodes in this net

% TIME-LAGGED AUTO/CROSS-COVARIANCE & CORRELATION FUNCTIONS
lastLag = 350; numLags = 1000;
lags = linspace(0, lastLag, numLags + 1);

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
Cov = nan(n, n, numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Cov(:,:,k) = Sigma_tau(tau);
end

% Theoretical correlation
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Corr = Cov ./ normMat_th;

Corr = Corr(idx, idx, :); % Select subnetwork TLC matrices
X = reshape(Corr, [], size(Corr, 3)); % Data matrix

% CROSS-LAG COVARIANCE (COSINE SIMILARITY)
Xth_cos = X ./ vecnorm(X, 2, 1); % Theoretical
G_cos = Xth_cos' * Xth_cos;
G_cos(tril(true(size(G_cos)), -1)) = NaN;

nexttile, imagesc(lags, lags, G_cos, 'AlphaData', ~isnan(G_cos));
axis square; colormap(magma); colorbar; clim([-1, 1])
title(sprintf('%s Cross-Lag Covariance $X^\\top X$', netName), 'Interpreter', 'latex', 'FontSize', 15);
xlabel('$\tau_k$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\tau_\ell$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
