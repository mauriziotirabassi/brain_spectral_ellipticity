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
A = subj.A; % effective connectivity
n = size(A, 1); I = eye(n);
Sigma_w = I * subj.output.eff_conn.NoiseVar; % noise covariance
Sigma = lyap(A, Sigma_w); % zero-lag covariance
S = (A * Sigma - Sigma * (A.')); % dC-Cov
hr = subj.h; % haemodynamic response

% Isolating inactive pairs (wrong)
pct = 99;
upper_percentile = prctile(S(:), pct);
lower_percentile = prctile(S(:), 100 - pct);
mask = (S >= upper_percentile) | (S <= lower_percentile);

% Network topology
S_plot = S; S_plot(~mask) = 0; % showtop(S_plot);

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

% State evolution in 2D/3D after PCA
% [coeff, score, latent, tsquared, explained, mu] = pca(x);
% x_reduced3 = score(:, 1:3); animate3(x_reduced3)
% x_reduced2 = score(:, 1:2); animate2(x_reduced2)

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

% Simulated covariance
Cov_sim0 = (x' * x) / size(x,1);
stds_emp = sqrt(diag(Cov_sim0));
normMat_emp = stds_emp * stds_emp';
T = size(x, 1);
Cov_sim = nan(n,n,numel(lags));
for k = 0:maxLag
    X = x(1:end-k, :);
    X_lag = x(1+k:end, :);
    Cov_sim(:,:,k+1) = (X_lag' * X) / (T - k);
end
Corr_sim = Cov_sim ./ normMat_emp;

% % Extend to negative lags based on property C_ij(-tau) = C_ji(tau)
% lags_full = [-fliplr(lags(2:end)), lags];
% Cov_th_neg = flip(Cov_th(:,:,2:end), 3); % Theoretical covariance
% Cov_th_neg = permute(Cov_th_neg, [2 1 3]);
% Cov_th_full = cat(3, Cov_th_neg, Cov_th);
% Corr_th_neg = flip(Corr_th(:,:,2:end), 3); % Theoretical correlation
% Corr_th_neg = permute(Corr_th_neg, [2 1 3]);
% Corr_th_full = cat(3, Corr_th_neg, Corr_th);
% Sigma_emp_neg = flip(Sigma_emp(:,:,2:end), 3); % Simulated correlation
% Sigma_emp_neg = permute(Sigma_emp_neg, [2 1 3]);
% Sigma_emp_full = cat(3, Sigma_emp_neg, Sigma_emp);

% CROSS-LAG CORRELATION (CLC)
% triuIdx = find(triu(ones(n), 1));
% triuIdx = find(triu(mask, 1)); % Isolate active pairs
% triuIdx = find(triu(mask | I, 0)); % Isolate active paris and keep diagonal

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags, lags, crosslag_corr(Corr_th));
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical CLCorr Subject %d', iSub));
nexttile, imagesc(lags, lags, crosslag_corr(Corr_sim));
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Simulated CLCorr Subject %d', iSub));

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags, lags, crosslag_cos(Corr_th));
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical CLCos Subject %d', iSub));
nexttile, imagesc(lags, lags, crosslag_cos(Corr_sim));
axis square; colorbar; colormap(magma); xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Simulated CLCos Subject %d', iSub));

%% SINGLE NODE CLC
% figure, tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
% clc_th = nan(length(lags),length(lags),numel(n));
% clc_emp = nan(length(lags),length(lags),numel(n));
% for i = 1:n
%     mask = false(size(A)); mask(:, i) = true;
%     clc_th(:,:,i) = crosslag_corr(Corr_th, mask);
%     clc_emp(:,:,i) = crosslag_corr(Corr_sim, mask);
%     nexttile(1), imagesc(mask);
%     nexttile(2), imagesc(lags, lags, clc_th(:,:,i))
%     nexttile(3), imagesc(lags, lags, clc_emp(:,:,i))
%     colorbar, colormap(magma)
%     drawnow, pause(0.2)
% end

%% META CLC
% figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
% mask = find(triu(ones(length(lags)), 0));
% nexttile, imagesc(crosslag_corr(clc_th, mask)), colorbar, colormap(magma)
% title('Theoretical Meta-CLC')
% nexttile, imagesc(crosslag_corr(clc_emp, mask)), colorbar, colormap(magma)
% title('Simulated Meta-CLC')

%% FUNCTIONS
function animate3(data)
    figure; 
    h = animatedline('LineWidth', 1);
    axis([min(data(:,1)) max(data(:,1)) ...
          min(data(:,2)) max(data(:,2)) ...
          min(data(:,3)) max(data(:,3))]);
    grid on
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    view(3)
    
    for i = 1:size(data, 1)
        addpoints(h, data(i, 1), data(i, 2), data(i, 3));
        drawnow
    end
end

function animate2(data)
    figure; 
    h = animatedline('LineWidth', 1);
    axis([min(data(:,1)) max(data(:,1)) ...
          min(data(:,2)) max(data(:,2))]);
    grid on
    xlabel('PC_1'); ylabel('PC_2');
    title('2D Trajectory Animation');
    
    for i = 1:size(data, 1)
        addpoints(h, data(i, 1), data(i, 2));
        drawnow
    end
end

function showtop(S)
    G = digraph(S);
    figure; h = plot(G, 'Layout','circle'); %, 'EdgeLabel',G.Edges.Weight);
    h.LineWidth = abs(G.Edges.Weight)/max(abs(G.Edges.Weight));
end

function C = crosslag_corr(matrixvec, mask)
%CROSSLAG_CORR Construct cross-lag covariance matrix (CLC) Pearson correaltion
%   matrixvec : series of matrices whose CLC to calculate
%   mask      : mask for the matrices
    vecs = reshape(matrixvec, [], size(matrixvec,3)); % vectorize each matrix
    if nargin >= 2 && ~isempty(mask), vecs = vecs(mask(:), :); end
    C = corr(vecs);
end

function C = crosslag_cos(matrixvec, mask)
%CROSSLAG_COS cosine similarity
    vecs = reshape(matrixvec, [], size(matrixvec,3));
    if nargin >= 2 && ~isempty(mask), vecs = vecs(mask(:), :); end
    C = (vecs' * vecs) ./ (sqrt(sum(vecs.^2,1))' * sqrt(sum(vecs.^2,1)));
end