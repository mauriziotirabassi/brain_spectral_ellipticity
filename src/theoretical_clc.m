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
G_raw  = X' * X; % m x m
G_raw(tril(true(size(G_raw)), -1)) = NaN;

% Cosine similarity across lags (normalize each column by its L2 norm)
colL2  = sqrt(sum(X.^2, 1)); % 1 x m
colL2(colL2 == 0) = eps;
G_cos  = (X' * X) ./ (colL2' * colL2); % m x m
G_cos(tril(true(size(G_cos)), -1)) = NaN;

% Pearson correlation across lags (demean columns then normalize)
% (equivalent to corr(X) up to tiny numerical noise)
Xc = X - mean(X, 1); % demean columns
colL2c = sqrt(sum(Xc.^2, 1)); % 1 x m (L2 of demeaned columns)
colL2c(colL2c == 0) = eps;
G_corr = (Xc' * Xc) ./ (colL2c' * colL2c); % m x m
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
[U, Sig, V] = svd(X, 'econ'); % V columns are eigenvectors of X'X

figure, tiledlayout(4, 1, 'TileSpacing','compact','Padding','compact');
for r = 1:4
    time_mode = V(:, r);
    nexttile(r), plot(lags, time_mode); xlabel('\tau'); title(sprintf('Temporal mode %d', r)); grid on;
end

diag(Sig).^2 ./ sum(diag(Sig).^2) % Energy captured

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
