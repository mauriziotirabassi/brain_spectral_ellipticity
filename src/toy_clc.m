clear; clc; %close all

rng(42); % This affects simulated noise injections
n = 2; I = eye(n);
Sigma_w = 100 * I; % Uncorrelated noise

% Topology
% topology = 'Ring';
% S = buildS(n, topology);
% showtop(S)

S = zeros(n);

% Initial energy distribution
Sigma = I; % Balanced
% eps = 1e-1; Sigma = eps * I; Sigma(5,5) = 10; % Unbalanced

% Dynamics
A = (-0.5 * Sigma_w + S) / Sigma;

% % Define A
% topology = 'Random';
% S = randn(n); S = 0.5 * (S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q' * Q; % Random negative definite part
% A = -Q + S;
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * (A.'));
% showtop(S);

ev = eig(A); fprintf('max real part = %.4g\n', max(real(ev)));

% Simulation parameters
n_time = 2e3; %1e4; % Simulation time steps
transient_length = 1e3;
tr = 0.01; % Sampling period <--------------------------------------------
t = 0 : tr : (n_time + transient_length - 1) * tr;

% SIMULATE W/STOCHASTIC INPUT
sys = ss(A, eye(n), eye(n), zeros(n));
w = mvnrnd(zeros(1, n), Sigma_w, length(t)); % Generate noise input
y_stoch = lsim(sys, w, t); % Simulate LTI response to noise
y_stoch = y_stoch(transient_length + 1:end,:); % Remove transient

% COVARIANCE MATRIX
maxLag = size(y_stoch, 1) / 2; % number of lags to evaluate
lags = (0:maxLag) * tr;

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Cov_th = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Cov_th(:,:,k) = Sigma_tau(tau);
end
Corr_th = Cov_th ./ normMat_th; % Pearson correlation

% Simulated covariance
Sigma_emp0 = (y_stoch' * y_stoch) / size(y_stoch,1);
stds_emp = sqrt(diag(Sigma_emp0));
normMat_emp = stds_emp * stds_emp';
T = size(y_stoch, 1);
Sigma_emp = nan(n,n,numel(lags));
for k = 0:maxLag
    X = y_stoch(1:end-k, :);
    X_lag = y_stoch(1+k:end, :);
    Sigma_emp(:,:,k+1) = (X_lag' * X) / (T - k);
end
Sigma_emp = Sigma_emp ./ normMat_emp;

% Extend to negative lags based on property C_ij(-tau) = C_ji(tau)
lags_full = [-fliplr(lags(2:end)), lags]; % Symmetric lag vector
Cov_th_neg = flip(Cov_th(:,:,2:end), 3);   % flip along 3rd dimension (lags)
Cov_th_neg = permute(Cov_th_neg, [2 1 3]); % swap i<->j for negative lags
Cov_th_full = cat(3, Cov_th_neg, Cov_th);  % concatenate along 3rd dimension

Corr_th_neg = flip(Corr_th(:,:,2:end), 3);
Corr_th_neg = permute(Corr_th_neg, [2 1 3]);
Corr_th_full = cat(3, Corr_th_neg, Corr_th);

Sigma_emp_neg = flip(Sigma_emp(:,:,2:end), 3);
Sigma_emp_neg = permute(Sigma_emp_neg, [2 1 3]);
Sigma_emp_full = cat(3, Sigma_emp_neg, Sigma_emp);

% CROSS-LAG COVARIANCE (CLC) <--------------------------------------
mask = find(triu(ones(n),0)); % Include all pairs

% Isolate only active pairs: not really useful here because the ones with
% dC-Cov null actually do not contribute, so no noise
% mask = find(triu(abs(S) > 0, 1));

% Isolate but keep the diagonal: not useful for analyzing periodic patterns
% due to causal connections. It just adds useless information along y=-x
% and their parallels.
% mask = find(triu(abs(S) > 0 | eye(size(S))));

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags_full, lags_full, crosslagcov(Corr_th_full, mask));
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
% title(sprintf('Theoretical %s', topology));

nexttile, imagesc(lags_full, lags_full, crosslagcov(Sigma_emp_full, mask));
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
% title(sprintf('Empirical %s', topology));

%% SINGLE ROWS
% figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
% clc = nan(length(lags_full),length(lags_full),numel(n));
% for i = 1:n
%     mask = false(size(A)); mask(:, i) = true;
%     clc(:,:,i) = crosslagcov(Sigma_th_full, mask);
%     nexttile(1), imagesc(mask);
%     nexttile(2), imagesc(lags_full, lags_full, clc(:,:,i))
%     colorbar, colormap(magma)
%     drawnow, pause(0.2)
% end
% mask = find(triu(ones(length(lags_full)), 0));
% figure, imagesc(crosslagcov(clc, mask)), colorbar, colormap(magma)

%% FUNCTIONS
function S = buildS(n, topology)
%BUILD S Construct skew-symmetric adjacency matrix S
%   n        : number of nodes
%   topology : 'chain', 'uni_chain', 'ring', 'star'

S = zeros(n);
switch lower(topology)
    case 'chain'  % chain
        for i = 1:n-1
            S(i,i+1) = 1;
            S(i+1,i) = -1;
        end
    case 'ring'
        for i = 1:n
            j = mod(i,n) + 1;
            S(i,j) = 1;
            S(j,i) = -1;
        end
    case 'hub'
        for j = 2:n
            S(1,j) = 1;
            S(j,1) = -1;
        end
    otherwise
        error('Unknown topology: %s', topology)
end
end

function showtop(S)
    G = digraph(S);
    figure; h = plot(G, 'Layout','circle', 'EdgeLabel',G.Edges.Weight);
    h.LineWidth = abs(G.Edges.Weight)/max(abs(G.Edges.Weight));
end

function clc = crosslagcov(matrixvec, mask)
%CLC Construct cross-lag covariance matrix (CLC)
%   matrixvec : series of matrices whose CLC to calculate
%   mask      : mask for the matrices
    vecs = reshape(matrixvec, [], size(matrixvec,3)); % vectorize each matrix
    vecs = vecs(mask(:), :); % apply mask
    clc = corr(vecs);
end