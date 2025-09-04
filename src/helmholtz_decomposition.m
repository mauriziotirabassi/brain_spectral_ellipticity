clear; clc; close all

rng(42);
n = 2;
I = eye(n);
Sigma_w = 0.1 * I; % Uncorrelated noise

% Define A
A = [-4 1; -1 -.5];
Sigma_w = eye(2)*0.1;
Sigma = lyap(A, Sigma_w); 

% Derive A
% S = randn(n);
% S = 0.5 * (S - S'); % Ensure skew-symmetry
% Sigma = [7 1; 1 1];
% % Sigma = I;
% A = (-1/2 * Sigma_w + S) / Sigma;

% Derive Sigma
% S = randn(n); S = 0.5*(S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q'*Q; % Random negative definite part
% A = -Q + S;            
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * A.');

ev = eig(A);
fprintf('max real part = %.4g\n', max(real(ev)));

% Helmholtz decomposition
ASigma = A*Sigma;
dA_Cov = -0.5*Sigma_w;
dC_Cov = 0.5*(ASigma - ASigma');

% Grid for plotting vector field
npts = 10;
[x1, x2] = meshgrid(linspace(-1,1,npts), linspace(-1,1,npts));
X = [x1(:)'; x2(:)'];

%% PLOT DECOMPOSITION
figure, tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile, plotvec(ASigma * X, X, 'r'), title(sprintf('A\\Sigma'));
nexttile, plotvec(dA_Cov * X, X, 'b'), title('dA-Cov');
nexttile, plotvec(dC_Cov * X, X, 'g'), title('dC-Cov');

%%
figure, plotvec(A * X, X, 'g'), title('dC-Cov');

%% SIGMA^-1 NORMALIZATION
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

skip = 3;
idx = 1:skip:size(ellipsoid_pts,2);
V_norm = A * ellipsoid_pts;
quiver(X(idx), Y(idx), V_norm(1,idx), V_norm(2,idx), 0, 'b')

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