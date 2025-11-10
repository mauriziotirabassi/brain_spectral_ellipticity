clear; clc; %close all

rng(42); % This affects simulated noise injections
n = 10; I = eye(n);
Sigma_w = 0.1 * I; % Uncorrelated noise

% Topology
topology = 'Ring';
S = buildS(n, topology);
showtop(S)

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
% mask = find(triu(ones(n),1)); % Include all pairs

% Isolate only active pairs: not really useful here because the ones with
% dC-Cov null actually do not contribute, so no noise
mask = find(triu(abs(S) > 0, 1));

% Isolate but keep the diagonal: not useful for analyzing periodic patterns
% due to causal connections. It just adds useless information along y=-x
% and their parallels.
% mask = find(triu(abs(S) > 0 | eye(size(S))));

% Plot
figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags_full, lags_full, crosslagcov(Sigma_th_full, mask));
axis square; colorbar; colormap(curl);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Theoretical %s', topology));

nexttile, imagesc(lags_full, lags_full, crosslagcov(Sigma_emp_full, mask));
axis square; colorbar; colormap(curl);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title(sprintf('Empirical %s', topology));

%% SINGLE ROWS
figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
clc = nan(length(lags_full),length(lags_full),numel(n));
for i = 1:n
    mask = false(size(A)); mask(:, i) = true;
    clc(:,:,i) = crosslagcov(Sigma_th_full, mask);
    nexttile(1), imagesc(mask);
    nexttile(2), imagesc(lags_full, lags_full, clc(:,:,i))
    colorbar, colormap jet
    drawnow, pause(0.2)
end
mask = find(triu(ones(length(lags_full)), 0));
figure, imagesc(crosslagcov(clc, mask)), colorbar, colormap jet

%%
% % 2D FOURIER TRANSFORM
% nLags = length(lags_full);
% df = 1 / (nLags * tr);
% freq_axis = (-floor(nLags/2):ceil(nLags/2)-1) * df; % in Hz
% 
% figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
% FCDfft_th = log(1 + abs(fftshift(fft2(FCD_th))));
% nexttile, imagesc(freq_axis, freq_axis, FCDfft_th); colormap(jet); colorbar;
% title(sprintf('Theoretical PSD %s', topology));
% 
% FCDfft_emp = log(1 + abs(fftshift(fft2(FCD_emp))));
% nexttile, imagesc(freq_axis, freq_axis, FCDfft_emp); colormap(jet); colorbar;
% title(sprintf('Empirical PSD %s', topology));
% 
% % SINOGRAM
% theta = 0:179;
% figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
% [R_th, xp_th] = radon(FCD_th, theta);
% nexttile, imagesc(theta, xp_th, R_th);
% xlabel('Angle (degrees)'); ylabel('Projection position');
% title(sprintf('Theoretical Sinogram %s', topology));
% 
% [R_emp, xp_emp] = radon(FCD_emp, theta);
% nexttile, imagesc(theta, xp_emp, R_emp);
% xlabel('Angle (degrees)'); ylabel('Projection position');
% title(sprintf('Empirical Sinogram %s', topology));

%% PATTERN FAMILY ISOLATION (HIERARCHICAL CLUSTERING)
D_th = squareform(1 - FCD_th); Z_th = linkage(D_th, 'average');
D_emp = squareform(1 - FCD_emp); Z_emp = linkage(D_emp, 'average');

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, dendrogram(Z_th); nexttile, dendrogram(Z_emp);

k = 4; % number of clusters you want
clusters_emp = cluster(Z_emp, 'maxclust', k);
clusters_th = cluster(Z_th, 'maxclust', k);

%%
X = vecs_th';   % size: nLags x nEdges
[coeff, score, latent, ~, explained] = pca(X);

% coeff:   edges × PCs   (families: which edges define each PC)
% score:   lags × PCs    (how each lag expresses those families)
% latent:  eigenvalues   (variance of each PC)
% explained: % variance explained per PC

% Example: plot first two PCs across lags (family activations)
figure; plot(lags_full, score(:,1:2));
xlabel('Lag index'); ylabel('Activation'); legend('PC1','PC2');
title('Family activations across lags');


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

    vecs = [];
    for k = 1:size(matrixvec, 3)
        Ck = matrixvec(:,:,k);
        vecs(:,k) = Ck(mask);
    end
    clc = corr(vecs);
end