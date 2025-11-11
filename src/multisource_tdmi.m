clear; clc; %close all
rng(42);

% Data
dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat'); % distance matrix
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat'); % structural conn matrix
outDir = fullfile(dataDir,'regressed_001_01_sim62131');
d = dir(fullfile(outDir,'*.mat')); % list of subject‐model .mat files
files = {d.name};

% Connectome
C = load(fullfile(structFile)).full_connectome_no_symm;
D = load(distFile,'mtx_euc_dis').mtx_euc_dis;

% Subject data
iSub = 5; % subject number
subj = load(fullfile(outDir,files{iSub}));
A = subj.A * 2; % effective connectivity
n = size(A, 1); I = eye(n);
Sigma_w = I * subj.output.eff_conn.NoiseVar; % noise covariance
Sigma = lyap(A, Sigma_w); % zero-lag covariance
S = (A * Sigma - Sigma * (A.')); % dC-Cov
hr = subj.h; % haemodynamic response

% Isolating inactive pairs
pct = 99;
upper_percentile = prctile(S(:), pct);
lower_percentile = prctile(S(:), 100 - pct);
mask = (S >= upper_percentile) | (S <= lower_percentile);

% Network topology
% S_plot = S; S_plot(~mask) = 0; showtop(S_plot);

% SIMULATION
% Simulation parameters
n_time = 2e3; %1e4; % Simulation time steps
transient_length = 1e3;
tr = 0.1; % Sampling period <---------------------------------------------
t = 0 : tr : (n_time + transient_length - 1) * tr;
t_sim = t(transient_length + 1:end);

% Simulate
sys = ss(A, eye(n), eye(n), zeros(n));
w = mvnrnd(zeros(1, n), Sigma_w, length(t)); % Generate noise input
x = lsim(sys, w, t); % Simulate LTI response to noise
x = x(transient_length + 1:end,:); % Remove transient

% COVARIANCE MATRIX
maxLag = size(x, 1) / 2;
lags = (0:maxLag) * tr;

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Cov_th = nan(n,n,numel(lags));
Corr_th = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Cov_th(:,:,k) = Sigma_tau(tau);
    Corr_th(:,:,k) = Sigma_tau(tau) ./ normMat_th; % Pearson correlation
end

% Simulated covariance
Sigma_emp0 = (x' * x) / size(x,1);
stds_emp = sqrt(diag(Sigma_emp0));
normMat_emp = stds_emp * stds_emp';
Sigma_emp = nan(n,n,numel(lags));
for k = 0:maxLag
    X = x(1:end-k, :);
    X_lag = x(1+k:end, :);
    Sigma_emp(:,:,k+1) = (X_lag' * X) / (size(x, 1) - k) ./ normMat_emp;
end

% Extend to negative lags based on property C_ij(-tau) = C_ji(tau)
lags_full = [-fliplr(lags(2:end)), lags]; % Symmetric lag vector
Cov_th_full = nan(n,n,numel(lags_full));
Corr_th_full = nan(n,n,numel(lags_full));
Sigma_emp_full = nan(n,n,numel(lags_full));
for i = 1:n
    for j = 1:n
        Covth_pos = squeeze(Cov_th(i,j,:));
        Covth_neg = squeeze(Cov_th(j,i,2:end));
        Cov_th_full(i,j,:) = [flipud(Covth_neg); Covth_pos];

        Corrth_pos = squeeze(Corr_th(i,j,:));
        Corrth_neg = squeeze(Corr_th(j,i,2:end));
        Corr_th_full(i,j,:) = [flipud(Corrth_neg); Corrth_pos];

        Cemp_pos = squeeze(Sigma_emp(i,j,:));
        Cemp_neg = squeeze(Sigma_emp(j,i,2:end));
        Sigma_emp_full(i,j,:) = [flipud(Cemp_neg); Cemp_pos];
    end
end

% TIME-DELAYED MUTUAL INFORMATION
I_tau = 0.5 * (Corr_th_full.^2 + Corr_th_full.^4 / 2); % Taylor approximation to second order
% I_tau = I_tau / log(2); % Conversion to bits
% I_tau = permute(I_tau, [2 1 3]); % I(j\to i) -> I(i \to j)

%% MULTISOURCE TIME-DELAYED MUTUAL INFORMATION
num_lags = numel(lags_full);
TDMI = zeros(n,num_lags);
red = nan(n,num_lags);
syn = nan(n,num_lags);
for j = 1:n
    inputs = find(S(:,j) > 0); % Nodes that influence target j
    if isempty(inputs)
        continue
    end

    Sigma_xx = Cov_th_full(inputs, inputs, 1); % Sources zero-lag cov
    Sigma_yy = Cov_th_full(j, j, 1); % Target zero-lag cov
    for k = 1:num_lags
        Sigma_xy = Cov_th_full(j, inputs, k)'; % Sources-target lagged cov
        TDMI(j,k) = 0.5 * log(det(Sigma_yy) / det(Sigma_yy - Sigma_xy' * pinv(Sigma_xx) * Sigma_xy)); % Gaussian multisource TDMI
        I_ind = I_tau(inputs, k); % TDMI of individual sources at lag k
        red(j,k) = min(I_ind); % minimum across sources
        syn(j,k) = TDMI(j,k) - sum(I_ind) + (numel(inputs) - 1) * red(j,k);
    end
end

%%

figure, tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags_full, lags_full, crosslagcov1(abs(TDMI)))
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical Multisource TDMI Subject %d', iSub));
nexttile, imagesc(lags_full, lags_full, crosslagcov1(abs(red)))
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical Multisource Redundancy Subject %d', iSub));
nexttile, imagesc(lags_full, lags_full, crosslagcov1(abs(syn)))
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical Multisource Synergy Subject %d', iSub));

%%

% CROSS-LAG COVARIANCE (CLC)
% triuIdx = find(triu(ones(n), 1));
triuIdx = find(triu(mask, 1)); % Isolate active pairs

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags_full, lags_full, crosslagcov2(Corr_th_full, triuIdx));
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical TLC Subject %d', iSub));
nexttile, imagesc(lags_full, lags_full, crosslagcov2(I_tau, triuIdx));
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical TDMI Subject %d', iSub));
% nexttile, imagesc(lags_full, lags_full, crosslagcov1(abs(TDMI)))
% axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
% title(sprintf('Theoretical Multisource TDMI Subject %d', iSub));
% nexttile, imagesc(lags_full, lags_full, crosslagcov(Sigma_emp_full, triuIdx));
% axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
% title(sprintf('Empirical Subject %d', iSub));

%% FUNCTIONS
function showtop(S)
    G = digraph(S);
    figure; h = plot(G, 'Layout','circle'); %, 'EdgeLabel',G.Edges.Weight);
    h.LineWidth = abs(G.Edges.Weight)/max(abs(G.Edges.Weight));
end

function clc = crosslagcov2(matrixvec, mask)
%CLC Construct cross-lag covariance matrix (CLC)
%   matrixvec : series of matrices whose CLC to calculate
%   mask      : mask for the matrices

    vecs = [];
    for k = 1:size(matrixvec, 3)
        Ck = matrixvec(:,:,k);
        vecs(:,k) = Ck(mask);
    end
    clc = corr(vecs);
end

function clc = crosslagcov1(vec, mask)
%CROSSLAGCOV Construct cross-lag covariance/correlation matrix from arrays
%   vec      : n_elements x num_lags
%   mask     : logical array selecting elements to include

    if nargin > 1
        vec = vec(mask,:);
    end

    clc = corr(vec);
end