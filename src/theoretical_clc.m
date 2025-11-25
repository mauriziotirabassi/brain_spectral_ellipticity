clear; clc; %close all
rng(42);

n = 2;
% A = diag([-1 -5]);
A = [-1 -1; 1 -1];
Sigma_w = eye(n);
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
maxLag = 450; % Number of lags to evaluate (excluding zero lag)
delta_tau = 0.05;
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

% TODO: Select single node lead or lag profile.

% mask2D = eye(n) > 0;
% mask3D = repmat(mask2D, 1, 1, size(Corr_th,3));
% X = Corr_th(mask3D);
X = squeeze(Corr_th(2, :, :)); % Single-node lag profile
% X = squeeze(Corr_th(:, 1, :)); % Single-node lead profile
X = reshape(X, [], size(Corr_th,3));

% TODO: Avoid division by 0 with eps.

%% CROSS-LAG COVARIANCE (DOT PRODUCT)
G_raw = X' * X; %/ (size(X, 1) - 1);
G_raw(tril(true(size(G_raw)), -1)) = NaN;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, stackedplot(lags, X.');
title('Auto/Cross-Covariance Functions'), xlabel('\tau')
nexttile, h = imagesc(lags, lags, G_raw);
axis square; colormap(magma); colorbar, %clim([-1 1]); % clim(clims)
xlabel('\tau_1'); set(h, 'AlphaData', ~isnan(G_raw)); % ylabel('\tau_2');
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Dot Product'))

%% CROSS-LAG COVARIANCE (COSINE SIMILARITY)
% Cosine similarity across lags (normalize each column by its L2 norm)
X_cos = X ./ vecnorm(X, 2, 1);
G_cos = X_cos' * X_cos; %/ (size(X, 1) - 1);
% D = pdist2(X', X', 'cosine'); % pairwise cosine distance
% G_cos = 1 - D; % convert distance to similarity
G_cos(tril(true(size(G_cos)), -1)) = NaN;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, stackedplot(lags, X_cos.');
title('L2-Normalized Auto/Cross-Covariance Functions'), xlabel('\tau')
nexttile, h = imagesc(lags, lags, G_cos);
axis square; colormap(magma); colorbar, clim([-1 1]); % clim(clims)
xlabel('\tau_1'); set(h, 'AlphaData', ~isnan(G_cos)); ylabel('\tau_2');
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Cosine Similarity'))

[~, ~, V_cos] = svd(X_cos, 'econ');
figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
for r = 1:3
    nexttile(r), plot(lags, V_cos(:,r));
    grid on; xlabel('\tau'); title(sprintf('Cosine mode %d', r));
end

%% CROSS-LAG COVARIANCE (PEARSON CORRELATION)
% Pearson correlation across lags (demean columns then normalize)
X_corr = (X - mean(X, 1)) ./ std(X, 0, 1);
G_corr = (X_corr' * X_corr); %/ (size(X, 1) - 1);
% G_corr = corr(X);
G_corr(tril(true(size(G_corr)), -1)) = NaN;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, stackedplot(lags, X_corr.');
title('Standardized Auto/Cross-Covariance Functions'), xlabel('\tau')
nexttile, h = imagesc(lags, lags, G_corr);
axis square; colormap(magma); colorbar; clim([-1 1]); % clim(clims) 
xlabel('\tau_1'); ylabel('\tau_2'); set(h, 'AlphaData', ~isnan(G_corr));
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Pearson Correlation'))

[~, ~, V_corr] = svd(X_corr, 'econ');
figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
for r = 1:3
    nexttile(r); plot(lags, V_corr(:,r));
    grid on; xlabel('\tau'); title(sprintf('Corr mode %d', r));
end

% clims = [ min([min(G_cos(:)), min(G_corr(:))]), ...
%           max([max(G_cos(:)), max(G_corr(:))]) ];

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