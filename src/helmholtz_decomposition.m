clear; clc; close all

rng(42);
n = 2;
Sigma_w = eye(n);
Sigma = diag([1 4]);
S = 1 * [0 1; -1 0];
A = (-0.5 * Sigma_w + S) / Sigma;

ev = eig(A);
fprintf('max real part = %.4g\n', max(real(ev)));

ASigma = A * Sigma;
dA_Cov = -0.5 * Sigma_w;
dC_Cov = 0.5 * (ASigma - ASigma');

% Grid for plotting vector field
npts = 11;
[x1, x2] = meshgrid(linspace(-1, 1, npts), linspace(-1, 1, npts));
X = [x1(:)'; x2(:)'];

%% A\SIGMA HELMHOLTZ DECOMPOSITION
figure, tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile, plotvec(ASigma * X, X, 'r'), title('$A\Sigma$', 'Interpreter', 'latex');
nexttile, plotvec(dA_Cov * X, X, 'b'), title('$-\frac{1}{2}\Sigma_w$', 'Interpreter', 'latex');
nexttile, plotvec(dC_Cov * X, X, 'g'), title('$S$', 'Interpreter', 'latex');

%% A HELMHOLTZ DECOMPOSITION
A_diss = -0.5 * Sigma_w;
A_rot = dC_Cov / Sigma;
[U,D] = eig(Sigma);
radii = sqrt(diag(D));

figure, tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile, plotvec(A * X, X, 'r'), title('$A$ ', 'Interpreter', 'latex', 'FontSize', 15, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
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

nexttile, plotvec(A_diss * X, X, 'b'), title('Dissipative Flow $-\frac{1}{2}\Sigma_w\Sigma^{-1}$', 'Interpreter', 'latex', 'FontSize', 15, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
ylabel('');
nexttile, plotvec(A_rot * X, X, 'g'), title('Conservative Flow $S\Sigma^{-1}$', 'Interpreter', 'latex', 'FontSize', 15, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
ylabel('');

%% A SCALAR–TRACELESS DECOMPOSITION
S_half = sqrtm(Sigma);
Atilde = inv(S_half) * A * S_half; 

mu      = trace(Atilde)/2;
gammaJ  = Atilde - mu*eye(n);    % = γ J̃
% gamma   = sqrt(abs(det(gammaJ)));
% Jtilde  = gammaJ / gamma;        % normalized deviatoric generator

J_diss = diag(diag(gammaJ));    % Dissipative (Diagonal/Stretching)
J_sol  = gammaJ - J_diss;       % Solenoidal (Off-diagonal/Rotation)

f = figure('Color', 'w', 'Position', [100 100 1000 500]);
t = tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding', 'compact');

% Row 1
nexttile(1); plotvec(A * X, X, 'k');
title('$\tilde A$', 'Interpreter','latex', 'FontSize', 14);
annotation('textbox', [0.3 0.746 0.05 0.05], 'String', '$=$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 20, 'HorizontalAlignment', 'center');

% Plot 2: Deviatoric (Remove x_2 label)
nexttile(2, [1 2]); plotvec(gammaJ * X, X, 'r');
title('Isochoric Flow $\gamma\tilde J$', 'Interpreter','latex', 'FontSize', 14);
ylabel(''); 
annotation('textbox', [0.68 0.746 0 0.05], 'String', '$+$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 20, 'HorizontalAlignment', 'center');

% Plot 3: Average Dissipation (Remove x_2 label)
nexttile(4); plotvec(mu * eye(2) * X, X, 'b');
title('Average Dissipation $\mu I$', 'Interpreter','latex', 'FontSize', 14);
ylabel('');

% Row 2
nexttile(5); axis off;

% Plot 4: Diff. Dissipation (Keep label)
nexttile(6); plotvec(J_diss*X, X, 'm');
title('Differential Dissipation', 'Interpreter','latex', 'FontSize', 12);
annotation('textbox', [0.5 0.25 0.01 0.05], 'String', '$+$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 20, 'HorizontalAlignment', 'center');

% Plot 5: Solenoidal (Remove x_2 label)
nexttile(7); plotvec(J_sol*X, X, 'g');
title('Solenoidal Flow', 'Interpreter','latex', 'FontSize', 12);
ylabel('');

% Brace Overlay
axes('Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
brace_width = 'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo';
text(0.5, 0.45, ['$\overbrace{\phantom{' brace_width '}}^{\phantom{a}}$'], 'Interpreter', 'latex', 'FontSize', 15, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% SOLENOIDAL FLOW WHITNENING
S_white = inv(Sigma^2) * dC_Cov * inv(Sigma^2);

figure, tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile, plotvec(dC_Cov * X, X, 'r'), title('$S$ ', 'Interpreter', 'latex')
nexttile, plotvec(A_rot * X, X, 'r'), title('$S\Sigma^{-1}$ ', 'Interpreter', 'latex')
nexttile, plotvec(S_white * X, X, 'r'), title('$\Sigma^{-1/2}S\Sigma^{-1/2}$ ', 'Interpreter', 'latex')

%% Sigma RESCALES A IN EIGENSPACE ()
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

%% FUNCTIONS
function plotvec(vec, X, color)
    quiver(X(1,:), X(2,:), vec(1,:), vec(2,:), color);
    xlabel('x_1'), ylabel('x_2'), axis equal
    xlim([-1 1]); ylim([-1 1]); grid on
end
