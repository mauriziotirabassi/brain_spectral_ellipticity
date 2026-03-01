clear; clc; close all
rng(42);

% --- Global Knobs ---
fs_title = 18;   % Font size for Titles
fs_label = 50;   % Font size for Axis Labels (x_1, x_2)
fs_annot = 25;   % Font size for Annotations (=, +)
fs_tick  = 14;   % Font size for Axis Ticks
fs_brace = 15;   % Font size specifically for the Brace

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

S_half = sqrtm(Sigma);
Atilde = inv(S_half) * A * S_half; 
mu      = trace(Atilde)/2;
gammaJ  = Atilde - mu*eye(n);    
J_diss = diag(diag(gammaJ));    
J_sol  = gammaJ - J_diss;       

f = figure('Color', 'w', 'Position', [100 100 1100 600]);
t = tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding', 'compact');

% Row 1
% Plot 1: Atilde
nexttile(1); 
plotvec(Atilde * X, X, 'k', fs_label, fs_tick);
title('$\tilde A$', 'Interpreter','latex', 'FontSize', fs_title);
annotation('textbox', [0.3 0.746 0.05 0.05], 'String', '$=$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');

% Plot 2: Isochoric Flow
nexttile(2, [1 2]); 
plotvec(gammaJ * X, X, 'r', fs_label, fs_tick);
title('Isochoric Flow $A_0$', 'Interpreter','latex', 'FontSize', fs_title);
ylabel(''); 
annotation('textbox', [0.68 0.746 0 0.05], 'String', '$+$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');

% Plot 3: Average Dissipation
nexttile(4); 
plotvec(mu * eye(2) * X, X, 'b', fs_label, fs_tick);
title('Average Dissipation $\mu I$', 'Interpreter','latex', 'FontSize', fs_title);
ylabel('');

% Row 2
nexttile(5); axis off;

% Plot 4: Differential Dissipation
nexttile(6); 
plotvec(J_diss*X, X, 'm', fs_label, fs_tick);
title('Differential Dissipation', 'Interpreter','latex', 'FontSize', fs_title);
annotation('textbox', [0.5 0.25 0.01 0.05], 'String', '$+$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');

% Plot 5: Solenoidal Flow
nexttile(7); 
plotvec(J_sol*X, X, 'g', fs_label, fs_tick);
title('Solenoidal Flow', 'Interpreter','latex', 'FontSize', fs_title);
ylabel('');

% Brace Overlay
% Note: The string length below matches the original code provided
brace_width = 'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo';
axes('Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
text(0.5, 0.45, ['$\overbrace{\phantom{' brace_width '}}^{\phantom{a}}$'], 'Interpreter', 'latex', 'FontSize', fs_brace, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

function plotvec(vec, X, color, fs_label, fs_tick)
    quiver(X(1,:), X(2,:), vec(1,:), vec(2,:), color);
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fs_label);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fs_label);
    axis equal
    xlim([-1 1]); ylim([-1 1]); grid on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs_tick);
end