%% SETUP
clearvars, close all; clc

dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% list of subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';         % make symmetric
n        = size(Dfull,1);
% mask     = ~eye(n);
mask      = triu(true(n),1); % working with upper triangular
dvals    = Dfull(mask);

% select subject number
iSub = 13;

% load this subject’s fitted A, TR, noise variance
subj = load(fullfile(outDir,files{iSub}));
A = subj.A;
Sigma_w = eye(n)*subj.output.eff_conn.NoiseVar;
Sigma = lyap(A,Sigma_w);
S = 0.5 * (A * Sigma - Sigma * (A.'));
Sigma_tau = @(t) expm(A*t)*Sigma;

% stds and normalization matrix
stds       = sqrt(diag(Sigma));     % n×1 vector of standard deviations
normMat    = stds * stds.';     % n×n matrix where (i,j)=σ_i*σ_j

%%
k       = 1;                         % which state to visualize (1..74)
T       = 10;                        % total time
dt      = 0.05;                      % time step
N       = 2000;                      % ensemble size (increase for smoother)
s0      = 5.0;                       % reference lag for slice plot

% Precompute
times = 0:dt:T;
M     = numel(times);
[~, refIdx] = min(abs(times - s0));  % index closest to s0

% Preallocate
X      = zeros(74, N);
X_store = zeros(M, 74, N);

% Euler-Maruyama simulation
sqrt_dt = sqrt(dt);
X_store(1, :, :) = X;
for ti = 2:M
    dW = sqrt_dt * (chol(Sigma_w) * randn(74, N));
    X  = X + A*X*dt + dW;
    X_store(ti, :, :) = X;
end

% Build covariance surface for variable k
C = zeros(M, M);
for i = 1:M
    x_i = squeeze(X_store(i, k, :));  % vector length N
    for j = 1:M
        x_j = squeeze(X_store(j, k, :));
        C(i,j) = (x_i - mean(x_i))'*(x_j - mean(x_j)) / N;
    end
end

% 3D Surface plot
[TT, SS] = meshgrid(times, times);
figure('Position',[100 100 1200 500]);
subplot(1,2,1);
surf(TT, SS, C, 'EdgeColor','none');
xlabel('t'); ylabel('s'); zlabel(sprintf('Cov[x_%d(t),x_%d(s)]',k,k));
title('Covariance Surface'); view(45,30)
colorbar

% Diagonal & Lag slice
diagCov = diag(C);
lagCov  = C(:, refIdx);

subplot(1,2,2);
plot(times, diagCov, 'LineWidth',1.8); hold on;
plot(times, lagCov,  'LineWidth',1.8);
xlabel('t'); ylabel('Covariance');
title(sprintf('Diagonal vs. Lag slice (s = %.2f)', times(refIdx)));
legend('Cov[x_k(t),x_k(t)]','Cov[x_k(t),x_k(s_0)]','Location','best');
grid on

% … your code up to computing Sigma_ss …

% Solve Lyapunov for Σ and compute AΣ
Sigma_ss = lyap(A, Sigma_w);            % solves A*Σ + Σ*A' = -Sigma_w

% Compute A*Sigma_ss and then index its (k,k) element
A_Sigma = A * Sigma_ss;                 % first assign to a variable
slope_kk = A_Sigma(k, k);               % now index into that matrix

% Annotate it on the plot
annotation('textbox',[0.65,0.2,0.2,0.1],...
           'String',sprintf('Zero‑lag slope ≈ %.3g', slope_kk),...
           'FitBoxToText','on','BackgroundColor','w');
