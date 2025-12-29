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
iSub = 2; % subject number
subj = load(fullfile(outDir,files{iSub}));
A = subj.A; % effective connectivity
n = size(A, 1); I = eye(n);
Sigma_w = I * subj.output.eff_conn.NoiseVar; % noise covariance
Sigma = lyap(A, Sigma_w); % zero-lag covariance
S = (A * Sigma - Sigma * (A.')); % dC-Cov
hr = subj.h; % haemodynamic response

% % DYNAMIC ANISOTROPY INDEX (DAI) SPECTRUM
[~, T] = schur(A, 'real');
kappa_spec = [];
i = 1;
while i <= n
    if i < n && abs(T(i + 1, i)) ~= 0 % 2x2 oscillatory Schur block
        b = T(i, i + 1); c = T(i + 1, i);
        kappa = (b^2 + c^2) / abs(b * c);
        kappa_spec(end + 1) = kappa; %#ok<SAGROW>
        i = i + 2;
    else, i = i + 1;
    end
end
% disp('Kappa spectrum:');
% disp(kappa_spec(:));
figure, histogram(kappa_spec, 10);
title('$\kappa$ Spectrum', 'Interpreter', 'latex')

% Unique networks & their indices
rname = fullfile(dataDir, 'regions.xlsx');
Tr = readtable(rname);
netLabels = Tr.NETWORK;
[uniqueNets, ~, ic] = unique(netLabels);

netId = 5;
netName = uniqueNets{netId};
idx = find(ic == netId); % rows for this network
m = numel(idx); % number of nodes in this net

% SIMULATION
% Simulation parameters
n_time = 1e3; %1e4; % Simulation time steps
transient_length = 1e3;
tr = 0.5; % Sampling period <---------------------------------------------
t = 0 : tr : (n_time + transient_length - 1) * tr;
t_sim = t(transient_length + 1:end);

% Simulation
sys = ss(A, eye(n), eye(n), zeros(n));
w = mvnrnd(zeros(1, n), Sigma_w, length(t)); % Generate noise input
x = lsim(sys, w, t); % Simulate LTI response to noise
x = x(transient_length + 1:end, :); % Remove transient

% BOLD response
bold = zeros(n_time, n);
for k = 1:n
    h = hr(:,k);
    h = h / sum(h); % normalize
    bold(:,k) = conv(x(:,k), h, "same");
end

% COVARIANCE MATRIX
maxLag = size(x, 1) / 2;
lags = (0:maxLag) * tr;

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
Cov_th = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Cov_th(:,:,k) = Sigma_tau(tau);
end
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Corr_th = Cov_th ./ normMat_th; % Pearson correlation

x = bold;

% Simulated covariance
Cov_sim0 = (x' * x) / size(x,1);
stds_emp = sqrt(diag(Cov_sim0));
normMat_emp = stds_emp * stds_emp';
Tm = size(x, 1);
Cov_sim = nan(n,n,numel(lags));
for k = 0:maxLag
    X_th = x(1:end-k, :);
    X_lag = x(1+k:end, :);
    Cov_sim(:,:,k+1) = (X_lag' * X_th) / (Tm - k);
end
Corr_sim = Cov_sim ./ normMat_emp;

% Select subnetwork covariance matrices
Corr_th = Corr_th(idx, idx, :);
Corr_sim = Corr_sim(idx, idx, :);

% Data matrix
X_th = reshape(Corr_th, [], size(Cov_th, 3));
X_sim = reshape(Corr_sim, [], size(Cov_sim, 3));

% CROSS-LAG COVARIANCE (COSINE SIMILARITY)
Xth_cos = X_th ./ vecnorm(X_th, 2, 1); % Theoretical
Gth_cos = Xth_cos' * Xth_cos;
Gth_cos(tril(true(size(Gth_cos)), -1)) = NaN;
Xsim_cos = X_sim ./ vecnorm(X_sim, 2, 1); % Simulated
Gsim_cos = Xsim_cos' * Xsim_cos;
Gsim_cos(tril(true(size(Gsim_cos)), -1)) = NaN;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, h = imagesc(lags, lags, Gth_cos); % Theoretical
axis square; colormap(magma); colorbar, %clim([-1 1]); % clim(clims)
xlabel('\tau_1'); set(h, 'AlphaData', ~isnan(Gth_cos)); ylabel('\tau_2');
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Theoretical CLC: %s', netName))
nexttile, h = imagesc(lags, lags, Gsim_cos); % Simulated
axis square; colormap(magma); colorbar, %clim([-1 1]); % clim(clims)
xlabel('\tau_1'); set(h, 'AlphaData', ~isnan(Gsim_cos)); ylabel('\tau_2');
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Simulated CLC: %s', netName))