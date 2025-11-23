clear; clc; %close all
rng(42);

n = 2; I = eye(n);
Sigma_w = I;
A = [-.02, -.1; .1, -.2];
Sigma = lyap(A, Sigma_w);
S = 0.5 * (A * Sigma - Sigma * (A.'));

% Topology
% topology = 'Chain'; S = buildS(n, topology); % showtop(S)
% topology = 'Hub'; S2 = buildS(n, topology); % showtop(S)
% S = S + 3 * S2; showtop(S)
% S = zeros(n);

% Dynamics
% Sigma = I;
% A = (-0.5 * Sigma_w + S) / Sigma;

% % Define A
% topology = 'Random';
% S = randn(n); S = 0.5 * (S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q' * Q; % Random negative definite part
% A = -Q + S;
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * (A.'));
% showtop(S);

% TIME-LAGGED AUTO/CROSS-COVARIANCE & CORRELATION FUNCTIONS
maxLag = 10000; % Number of lags to evaluate (excluding zero lag)
delta_tau = 0.01;
lags = (0:maxLag) * delta_tau;

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
Cov_th = nan(n, n, numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Cov_th(:,:,k) = Sigma_tau(tau);
end

% Theoretical correlation
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Corr_th = Cov_th ./ normMat_th;

figure;
if n > 4, npl = 5; else, npl = n; end
for i = 1:npl
    for j = 1:npl
        subplot(npl, npl, (i - 1) * npl + j);
        plot(lags, squeeze(Corr_th(i, j, :)), 'b-');
        title(sprintf('\\Sigma_{%d%d}(\\tau)', i, j));
        xlabel('\tau'); grid on;
    end
end
sgtitle('Auto/Cross-Covariance Functions');

% TODO: Select single node lead or lag profile.

X = reshape(Corr_th, [], size(Corr_th,3));

%% CROSS-LAG COVARIANCE
% Raw Gram across lags (no centering, no scaling)
G_raw = X' * X / (size(X, 1) -1); % m x m
G_raw(tril(true(size(G_raw)), -1)) = NaN;

% Cosine similarity across lags (normalize each column by its L2 norm)
X_cos = X ./ vecnorm(X, 2, 1);
G_cos = X_cos' * X_cos;
% D = pdist2(X', X', 'cosine'); % pairwise cosine distance
% G_cos = 1 - D; % convert distance to similarity
G_cos(tril(true(size(G_cos)), -1)) = NaN;

% Pearson correlation across lags (demean columns then normalize)
X_corr = (X - mean(X, 1)) ./ std(X, 0, 1);
G_corr = (X_corr' * X_corr) / (size(X, 1) - 1);
% G_corr = corr(X);
G_corr(tril(true(size(G_corr)), -1)) = NaN;

figure, tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile, h = imagesc(lags, lags, G_raw); % G_raw
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2'); set(h, 'AlphaData', ~isnan(G_raw));
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Dot Product'))
nexttile, h = imagesc(G_cos); % G_cos
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2'); set(h, 'AlphaData', ~isnan(G_cos));
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Cosine Similarity'))
nexttile, h = imagesc(G_corr); % G_cos
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2'); set(h, 'AlphaData', ~isnan(G_corr));
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Pearson Correlation'))

%% TIME-LAGGED COVARIANCE EIGENFUNCTIONS
% Dot product
[~, ~, V_raw] = svd(X, 'econ');
figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
for r = 1:3
    nexttile(r), plot(lags, V_raw(:,r));
    grid on; xlabel('\tau'); title(sprintf('Raw mode %d', r));
end

% Cosine similarity
[~, ~, V_cos] = svd(X_cos, 'econ');
figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
for r = 1:3
    nexttile(r), plot(lags, V_cos(:,r));
    grid on; xlabel('\tau'); title(sprintf('Cosine mode %d', r));
end

% Pearson correlation
[~, ~, V_corr] = svd(X_corr, 'econ');
figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
for r = 1:3
    nexttile(r); plot(lags, V_corr(:,r));
    grid on; xlabel('\tau'); title(sprintf('Corr mode %d', r));
end

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
