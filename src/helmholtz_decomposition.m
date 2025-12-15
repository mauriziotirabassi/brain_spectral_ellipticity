clear; clc; close all

rng(42);
n = 2;
Sigma_w = eye(n);
Sigma = diag([1 4]);
S = 1 * [0 1; -1 0];
A = (-0.5 * Sigma_w + S) / Sigma;

ev = eig(A);
fprintf('max real part = %.4g\n', max(real(ev)));

% Grid for plotting vector field
npts = 10;
[x1, x2] = meshgrid(linspace(-1, 1, npts), linspace(-1, 1, npts));
X = [x1(:)'; x2(:)'];

% A\SIGMA HELMHOLTZ DECOMPOSITION
ASigma = A * Sigma;
dA_Cov = -0.5 * Sigma_w;
dC_Cov = 0.5 * (ASigma - ASigma');

figure, tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile, plotvec(ASigma * X, X, 'r'), title('$A\Sigma$', 'Interpreter', 'latex');
nexttile, plotvec(dA_Cov * X, X, 'b'), title('$-\frac{1}{2}\Sigma_w$', 'Interpreter', 'latex');
nexttile, plotvec(dC_Cov * X, X, 'g'), title('$S$', 'Interpreter', 'latex');

% A HELMHOLTZ DECOMPOSITION
A_diss = -0.5 * Sigma_w;
A_rot = dC_Cov / Sigma;

figure, tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile, plotvec(A * X, X, 'r'), title('$A$ ', 'Interpreter', 'latex')
nexttile, plotvec(A_diss * X, X, 'b'), title('Pure Dissipation $-\frac{1}{2}\Sigma_w\Sigma^{-1}$', 'Interpreter', 'latex')
nexttile, plotvec(A_rot * X, X, 'g'), title('Rotational Drift $S\Sigma^{-1}$', 'Interpreter', 'latex');

% A GEOMETRIC DECOMPOSITION
mu = trace(A)/2;
A_iso = mu * eye(n);
J_normalized = A - A_iso;
A_diff = J_normalized - A_rot;

figure, tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile, plotvec(A_iso * X, X, 'r'), title('Isotropic Dissipation $\mu I$ ', 'Interpreter', 'latex')
nexttile, plotvec(A_diff * X, X, 'b'), title('Differential Dissipation', 'Interpreter', 'latex')
nexttile, plotvec(A_rot * X, X, 'g'), title('Rotational Drift $S\Sigma^{-1}$', 'Interpreter', 'latex');

% Sigma RESCALES A IN EIGENSPACE
% This is the same thing as just plotting A vs ASigma
[U,D] = eig(Sigma);
radii = sqrt(diag(D));
A_tilde = U' * A * U;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, plotvec(A * X, X, 'r'), title('$A$', 'Interpreter', 'latex');
hold on

% Covariance ellipse
scale = .4;
theta = linspace(0,2*pi,100);
ellipse = U * (radii .* [cos(theta); sin(theta)]) * scale;
plot(ellipse(1,:), ellipse(2,:), 'g-', 'LineWidth', 1);

% Scaled eigenvectors
for i = 1:2
    vec = U(:,i) * radii(i) * scale;   % scale eigenvector by corresponding sqrt(eigenvalue)
    h = quiver(0, 0, vec(1), vec(2), 0, 'g-', 'LineWidth', 1, 'MaxHeadSize', 0.5);
end
legend(h, '$\Sigma$ Eigvecs', 'interpreter', 'latex')

nexttile, plotvec(A * Sigma * X, X, 'b'), title('$A\Sigma$', 'Interpreter', 'latex');

%% SIGMA^-1 NORMALIZATION (DEPRECATED)
[U,D] = eig(Sigma);
radii = sqrt(diag(D));

% Simple circle
theta = linspace(0, 2*pi, 100);
xc = cos(theta);
yc = sin(theta);

% Ellipsoid
ellipsoid_pts = U * diag(radii) * [xc; yc];
X = ellipsoid_pts(1,:);
Y = ellipsoid_pts(2,:);

figure; hold on; grid on; axis equal
plot(X, Y, 'k-', 'LineWidth', 1);

% Scaled eigenvectors
for i = 1:2
    vec = U(:,i) * radii(i);   % scale eigenvector by corresponding sqrt(eigenvalue)
    quiver(0, 0, vec(1), vec(2), 0, '-', 'LineWidth', 1, 'MaxHeadSize', 0.5);
end

% Vector field points on top of ellipsoid
V = ASigma * ellipsoid_pts;
quiver(X, Y, V(1,:), V(2,:), 0, 'r')

skip = 1;
idx = 1:skip:size(ellipsoid_pts,2);
V_norm = A * ellipsoid_pts;
scale = .05;
quiver(X(idx), Y(idx), scale*V_norm(1,idx), scale*V_norm(2,idx), 0, 'b')

%% SIGMA^-1/2 NORMALIZATION
Sigma_invhalf = U * diag(1./sqrt(diag(D))) * U';

% Whiten points
whitened_pts = Sigma_invhalf * ellipsoid_pts;

% Plot
figure; hold on; axis equal; grid on
plot(whitened_pts(1,:), whitened_pts(2,:), 'k-')

%% FUNCTIONS
function plotvec(vec, X, color)
    quiver(X(1,:), X(2,:), vec(1,:), vec(2,:), color);
    xlabel('x_1'), ylabel('x_2'), axis equal
    xlim([-1 1]); ylim([-1 1]); grid on
end