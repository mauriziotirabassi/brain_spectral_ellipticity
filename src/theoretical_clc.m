clear; clc; %close all
rng(42);

n = 2; Sigma_w = eye(n);

% Sigma = I; % scalar
Sigma = diag([1 40]); % diagonal
% Sigma = [2 0.6 0.3; 0.6 1.5 0.5; 0.3 0.5 1.8]; % full

omega = 2;
S = omega * [0 1; -1 0]; % 2D
% S = omega * [0 -1 0; 1 0 0; 0 0 0]; % 3D
% S = omega * [0 -1  .2 -.3; 1  0 -.4  .5; -.2 .4  0 -.6; .3 -.5 .6  0]; % 4D
% S = omega * [0 .1 .2 .3 .4; -.1 0 .5 .2 .7; -.2 -.5 0 0 0; -.3 -.2 0 0 .1; -.4 -.7 0 -.1 0]; % 5D

A = (-0.5 * Sigma_w + S) / Sigma;

% DYNAMIC ANISOTROPY INDEX (DAI)
mu = trace(A) / 2; Delta = trace(A)^2 - 4 * det(A);
gamma = sqrt(abs(Delta)) / 2;
J = (A - mu * eye(2)) / gamma;
kappa = trace(J' * J);
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
lastLag = 50; numLags = 5000;
% lags = linspace(0, lastLag, numLags + 1);
lags = linspace(-lastLag, lastLag, numLags + 1);

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

X = reshape(Corr_th, [], size(Corr_th, 3)); % data matrix

% CROSS-LAG COVARIANCE (COSINE SIMILARITY)
X_cos = X ./ vecnorm(X, 2, 1); G_cos = X_cos' * X_cos;
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
plot(u1(1), u2(1), 'go','MarkerFaceColor', 'g'); % start
plot(u1(end), u2(end), 'ro','MarkerFaceColor', 'r'); % end
axis equal; grid on;
xlabel('$u_1 = \sqrt2\,C(\gamma\tau)$', 'Interpreter', 'latex');
ylabel('$u_2 = \sqrt\kappa\,S(\gamma\tau)$', 'Interpreter', 'latex');
title(['$\omega = ', num2str(omega), ...
       ',\ \Delta = ', num2str(Delta), ...
       ',\ \kappa = ', num2str(kappa), '$'], 'Interpreter', 'latex');
