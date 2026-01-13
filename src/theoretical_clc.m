clear; clc; %close all
rng(42);

% set(groot, 'defaultTextInterpreter', 'latex');
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'defaultLegendInterpreter', 'latex');

n = 2; Sigma_w = eye(n);

% Sigma = eye(n); % scalar
Sigma = diag([1 10]); % diagonal
% Sigma = [2 0.6 0.3; 0.6 1.5 0.5; 0.3 0.5 1.8]; % full

d = diag(inv(Sigma)); sigma_w_scalar = 1; alpha = (sigma_w_scalar .* d) / 2; 
w_crit = abs(alpha(1) - alpha(2)) / (2 * sqrt(d(1)*d(2)));

omega = w_crit;
S = omega * skewone(n);
% S = omega * buildS(n, 'hub');
% S = omega * [0 -1  .2 -.3; 1  0 -.4  .5; -.2 .4  0 -.6; .3 -.5 .6  0]; % 4D
% S = omega * [0 1 0 0 0; -1 0 0 0 0; 0 0 0 1 0; 0 0 -1 0 0; 0 0 0 0 0]; % 5D

% w1 = 1; w2 = 2;
% w1 = 3; w2 = 3.5;
% w1 = 1; w2 = sqrt(5);
% S = [0  w1 0  0; 
%     -w1 0  0  0; 
%      0  0  0  w2; 
%      0  0 -w2 0];
% Scramble S to make complex topology, CLC will look like the independent
% modes case anyway.
% [Q, ~] = qr(randn(n)); S = Q * S * Q';

A = (-0.5 * Sigma_w + S) / Sigma;

% [~, T] = schur(A, 'real');

S_half = sqrtm(Sigma);
Atilde = inv(S_half) * A * S_half; 
mu = trace(Atilde)/2;
gammaJ = Atilde - mu*eye(n);
Delta = trace(Atilde)^2 - 4 * det(Atilde);
gamma = sqrt(abs(Delta))/2;
J = gammaJ ./ gamma;

% % DYNAMIC ANISOTROPY INDEX (DAI)
kappa = 0.5 * trace(J' * J);
disp(['Delta = ', num2str(Delta)]);
disp(['Kappa = ', num2str(kappa)]);

% DYNAMIC ANISOTROPY INDEX (DAI) SPECTRUM
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

% TIME-LAGGED AUTO/CROSS-COVARIANCE & CORRELATION FUNCTIONS
lastLag = 25; numLags = 500;
lags = linspace(0, lastLag, numLags + 1);
% lags = linspace(-lastLag, lastLag, numLags + 1);

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

X = reshape(Corr_th, [], size(Corr_th, 3)); % Data matrix

% CROSS-LAG COVARIANCE (COSINE SIMILARITY)
X_cos = X ./ vecnorm(X, 2, 1); G_cos = X_cos' * X_cos;
G_cos(tril(true(size(G_cos)), -1)) = NaN;

figure('Color', 'w'); tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile, [r, c] = ind2sub([n, n], 1:n^2);
lbls = "\rho_{" + r + c + "}"; n_limit = min(n^2, 25);
s = stackedplot(lags, X_cos(1:n_limit, :).', 'DisplayLabels', lbls(1:n_limit), ...
    'LineWidth', 1.2, 'GridVisible', 'on');
title('Data Marix $X$'), xlabel('\tau'), s.FontSize = 12;

nexttile, imagesc(lags, lags, G_cos, 'AlphaData', ~isnan(G_cos));
axis square; colormap(magma); colorbar; clim([-1, 1])
title('Cross-Lag Covariance $X^\top X$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('$\tau_k$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\tau_\ell$', 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');

% [U_cos, ~, V_cos] = svd(X_cos, 'econ');
% figure, tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
% for r = 1:3
%     nexttile(r), plot(lags, V_cos(:,r));
%     grid on; xlabel('\tau'); title(sprintf('Cosine mode %d', r));
% end

% showtop(S);

%%

% LAG-VECTOR (SIMILARITY STRUCTURE)
tol = 1e-12;
if Delta < -tol % oscillatory → ellipse
    C = cos(gamma * lags);
    S = sin(gamma * lags);
elseif Delta > tol % overdamped → hyperbola
    C = cosh(gamma * lags);
    S = sinh(gamma * lags);
else % critical → line
    C = ones(size(lags));
    S = tau;
end
u1 = sqrt(2) * C; u2 = sqrt(kappa) * S; % lag vector

figure, plot(u1, u2,'b-','LineWidth', 1); hold on
plot(-u1, u2, 'b--', 'LineWidth', 1);
axis equal; grid on;
xlabel('$u_1 = \sqrt2\,C(\gamma\tau)$', 'Interpreter', 'latex');
ylabel('$u_2 = \sqrt\kappa\,S(\gamma\tau)$', 'Interpreter', 'latex');
% title(['$\omega = ', num2str(omega), ...
%        ',\ \Delta = ', num2str(Delta), ...
%        ',\ \kappa = ', num2str(kappa), '$'], 'Interpreter', 'latex');

function S = skewone(n)
    U = triu(ones(n), 1); S = U - U';
end

function showtop(S)
    G = digraph(S);
    figure; h = plot(G, 'Layout','circle', 'EdgeLabel',G.Edges.Weight);
    h.LineWidth = abs(G.Edges.Weight)/max(abs(G.Edges.Weight));
end

function S = buildS(n, topology)
%BUILD S Construct skew-symmetric adjacency matrix S
%   n        : number of nodes
%   topology : 'chain', 'ring', 'hub'

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