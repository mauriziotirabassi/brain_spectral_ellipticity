clear; clc; %close all
rng(42);

n = 2;
% A = diag([-1 -2]);
% A = [-1 -1; 1 -2.5]; % Bulging
% A = [-1 1; -1 -3];
Sigma_w = eye(n);
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * (A.'));

% Topology
% topology = 'Chain'; S = buildS(n, topology); % showtop(S)
% topology = 'Hub'; S2 = buildS(n, topology); % showtop(S)
% S = S + 3 * S2; showtop(S)
% S = zeros(n);

% Dynamics
% Sigma = I; % scalar
Sigma = diag([1 40]); % general diagonal
% Sigma = [2 0.6 0.3; 0.6 1.5 0.5; 0.3 0.5 1.8]; % anisotropic

% Sigma = diag([0.5, 1]);   S = [0 0.1; -0.1 0]; % overdamped
% Sigma = diag([0.5, 1]);   S = [0 1; -1 0]; % oscillatory

% S = zeros(n); % pure dissipation
omega = 100;
S = omega * [0 1; -1 0]; % 2D
% S = .8 * [0 -1 0; 1 0 0; 0 0 0]; % 3D
% S = 1 * [0 -1  .2 -.3; 1  0 -.4  .5; -.2 .4  0 -.6; .3 -.5 .6  0]; % 4D
% S = 10 * [0 .1 .2 .3 .4; -.1 0 .5 .2 .7; -.2 -.5 0 0 0; -.3 -.2 0 0 .1; -.4 -.7 0 -.1 0]; % 5D
A = (-0.5 * Sigma_w + S) / Sigma;

% 2D KAPPA
mu = trace(A)/2;
Delta = trace(A)^2 - 4*det(A);
gamma = sqrt(abs(Delta))/2;
J = (A - mu*eye(2)) / gamma;
kappa = trace(J'*J);
disp(['Delta = ', num2str(Delta)]);
disp(['Kappa = ', num2str(kappa)]);

% DYNAMIC ANISOTROPY COEFFICIENT (DAC) SPECTRUM
% [~, T] = schur(A, 'real');
% kappa_spec = [];
% i = 1;
% while i <= n
%     if i < n && abs(T(i + 1, i)) ~= 0 % 2x2 oscillatory Schur block
%         b = T(i, i + 1); c = T(i + 1, i);
%         kappa = (b^2 + c^2) / abs(b * c);
%         kappa_spec(end + 1) = kappa; %#ok<SAGROW>
%         i = i + 2;
%     else, i = i + 1;
%     end
% end
% disp('Kappa spectrum:');
% disp(kappa_spec(:));

% % Define A
% topology = 'Random';
% S = randn(n); S = 0.5 * (S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q' * Q; % Random negative definite part
% A = -Q + S;
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * (A.'));
% showtop(S);

% TIME-LAGGED AUTO/CROSS-COVARIANCE & CORRELATION FUNCTIONS
lastLag = .5; numLags = 30000;
lags = linspace(0, lastLag, numLags + 1);

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
% mask3D = repmat(mask2D, 1, 1, size(Corr_th, 3));
% X = Corr_th(mask3D);
% X = squeeze(Corr_th(1, :, :)); % Single-node lag profile
% X = squeeze(Corr_th(:, 1, :)); % Single-node lead profile
X = reshape(Corr_th, [], size(Corr_th, 3));

% TODO: Avoid division by 0 with eps.

% CROSS-LAG COVARIANCE (DOT PRODUCT)
% G_raw = X' * X; %/ (size(X, 1) - 1);
% G_raw(tril(true(size(G_raw)), -1)) = NaN;
% 
% figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
% nexttile, stackedplot(lags, X.');
% title('Auto/Cross-Covariance Functions'), xlabel('\tau')
% nexttile, h = imagesc(lags, lags, G_raw);
% axis square; colormap(magma); colorbar, %clim([-1 1]); % clim(clims)
% xlabel('\tau_1'); set(h, 'AlphaData', ~isnan(G_raw)); % ylabel('\tau_2');
% set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
% title(sprintf('Dot Product'))

% CROSS-LAG COVARIANCE (COSINE SIMILARITY)
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
axis square; colormap(magma); colorbar, %clim([-1 1]); % clim(clims)
xlabel('\tau_1'); set(h, 'AlphaData', ~isnan(G_cos)); ylabel('\tau_2');
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title(sprintf('Cosine Similarity'))

% [U_cos, ~, V_cos] = svd(X_cos, 'econ');
% figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
% for r = 1:3
%     nexttile(r), plot(lags, V_cos(:,r));
%     grid on; xlabel('\tau'); title(sprintf('Cosine mode %d', r));
% end

% KAPPA
tol = 1e-12;
if Delta < -tol       % oscillatory → ellipse
    C = cos(gamma * lags);
    S = sin(gamma * lags);
elseif Delta > tol    % overdamped → hyperbola
    C = cosh(gamma * lags);
    S = sinh(gamma * lags);
else                  % critical → line
    C = ones(size(lags));
    S = tau;
end

% Lag vector
u1 = sqrt(2) * C;
u2 = sqrt(kappa) * S;

% Plot
figure;
plot(u1, u2,'b-','LineWidth',1); hold on
plot(u1(1), u2(1),'go','MarkerFaceColor','g'); % start
plot(u1(end), u2(end),'ro','MarkerFaceColor','r'); % end
axis equal; grid on;
xlabel('$u_1 = \sqrt2\,C(\gamma\tau)$', 'Interpreter', 'latex');
ylabel('$u_2 = \sqrt\kappa\,S(\gamma\tau)$', 'Interpreter', 'latex');
title(['$\omega = ', num2str(omega), ...
       ',\ \Delta = ', num2str(Delta), ...
       ',\ \kappa = ', num2str(kappa), '$'], 'Interpreter', 'latex');

%% CROSS-LAG COVARIANCE (PEARSON CORRELATION)
% Pearson correlation across lags (demean columns then normalize)
X_corr = (X - mean(X, 1)) ./ std(X, 0, 1);
G_corr = (X_corr' * X_corr) / (size(X, 1) - 1);
% G_corr = corr(X);
G_corr(tril(true(size(G_corr)), -1)) = NaN;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, stackedplot(lags, X_corr.');
title('Standardized Auto/Cross-Covariance Functions'), xlabel('\tau')
nexttile, h = imagesc(lags, lags, G_corr);
axis square; colormap(magma); colorbar;% clim([-1 1]); % clim(clims) 
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