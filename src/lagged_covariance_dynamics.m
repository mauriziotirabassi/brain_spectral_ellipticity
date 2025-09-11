clear; clc; %close all

rng(42); % This affects simulated noise injections
n = 10; I = eye(n);
Sigma_w = 0.1 * I; % Uncorrelated noise

% Topology
topology = 'Ring';
S = buildS(n, topology);
showS(S)

% Initial energy distribution
Sigma = I; % Balanced
% eps = 1e-1; Sigma = eps * I; Sigma(1,1) = 10; % Unbalanced

% Dynamics
A = (-0.5 * Sigma_w + S) / Sigma;
ev = eig(A); fprintf('max real part = %.4g\n', max(real(ev)));

% Define A
% topology = 'Random';
% S = randn(n); S = 0.5 * (S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q' * Q; % Random negative definite part
% A = -Q + S;
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * (A.'));
% showS(S);

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
Sigma_th = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Sigma_th(:,:,k) = Sigma_tau(tau) ./ normMat_th; % Pearson correlation
end

% Empirical covariance
Sigma_emp0 = (y_stoch' * y_stoch) / size(y_stoch,1);
stds_emp = sqrt(diag(Sigma_emp0));
normMat_emp = stds_emp * stds_emp';
Sigma_emp = nan(n,n,numel(lags));
for k = 0:maxLag
    X = y_stoch(1:end-k, :);
    X_lag = y_stoch(1+k:end, :);
    Sigma_emp(:,:,k+1) = (X_lag' * X) / (size(y_stoch, 1) - k) ./ normMat_emp;
end

% Extend to negative lags based on property C_ij(-tau) = C_ji(tau)
lags_full = [-fliplr(lags(2:end)), lags]; % Symmetric lag vector
Sigma_th_full = nan(n,n,numel(lags_full));
Sigma_emp_full = nan(n,n,numel(lags_full));
for i = 1:n
    for j = 1:n
        Cth_pos = squeeze(Sigma_th(i,j,:));
        Cth_neg = squeeze(Sigma_th(j,i,2:end));
        Sigma_th_full(i,j,:) = [flipud(Cth_neg); Cth_pos];

        Cemp_pos = squeeze(Sigma_emp(i,j,:));
        Cemp_neg = squeeze(Sigma_emp(j,i,2:end));
        Sigma_emp_full(i,j,:) = [flipud(Cemp_neg); Cemp_pos];
    end
end

% FCD-STYLE SIMILATIRY ACROSS LAGS <--------------------------------------
% triuIdx = find(triu(ones(n),1)); % Include all pairs 
triuIdx = find(triu(abs(S) > 0, 1)); % Isolate only active pairs
% triuIdx = find(triu(abs(S) > 0 | eye(size(S)))); % Isolate but keep the diagonal

% Theoretical
vecs_th = [];
for k = 1:length(lags_full)
    Ck = Sigma_th_full(:,:,k);
    vecs_th(:,k) = Ck(triuIdx);
end
FCD_th = corr(vecs_th);

% Empirical
vecs_em = [];
for k = 1:length(lags_full)
    Ck = Sigma_emp_full(:,:,k);
    vecs_em(:,k) = Ck(triuIdx);
end
FCD_emp = corr(vecs_em);

% Plot
figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags_full, lags_full, FCD_th);
axis square; colorbar; colormap jet;
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical %s', topology));

nexttile, imagesc(lags_full, lags_full, FCD_emp);
axis square; colorbar; colormap jet;
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Empirical %s', topology));

% 2D FOURIER TRANSFORM
nLags = length(lags_full);
df = 1 / (nLags * tr);
freq_axis = (-floor(nLags/2):ceil(nLags/2)-1) * df; % in Hz

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
FCDfft_th = log(1 + abs(fftshift(fft2(FCD_th))));
nexttile, imagesc(freq_axis, freq_axis, FCDfft_th); colormap(jet); colorbar;
title(sprintf('Theoretical PSD %s', topology));

FCDfft_emp = log(1 + abs(fftshift(fft2(FCD_emp))));
nexttile, imagesc(freq_axis, freq_axis, FCDfft_emp); colormap(jet); colorbar;
title(sprintf('Empirical PSD %s', topology));

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
    case 'star'
        for j = 2:n
            S(1,j) = 1;
            S(j,1) = -1;
        end
    otherwise
        error('Unknown topology: %s', topology)
end
end

function showS(S)
    G = digraph(S);
    figure; h = plot(G, 'Layout','circle', 'EdgeLabel',G.Edges.Weight);
    h.LineWidth = abs(G.Edges.Weight)/max(abs(G.Edges.Weight));
end