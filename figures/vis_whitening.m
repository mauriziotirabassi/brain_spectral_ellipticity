clear; clc; close all
rng(42);

% --- Global Knobs ---
fs_title = 20;   % Font size for Titles
fs_label = 40;   % Font size for Axis Labels (x_1, x_2, z_1, z_2)
fs_annot = 35;   % Font size for Annotations (=, +)
fs_tick  = 14;   % Font size for Axis Ticks

% --- Annotation Position Knobs ---
eq_x = 0.33; eq_y = 0.466;  % Position for '='
pl_x = 0.637; % Position for '+'

% --- Contour Knobs ---
contour_levels = [0.05, 0.225, 0.5]; % Provide array of values, or a single integer for N slices
contour_style  = 'g-';                  % 'g--' for dashed, 'g-' for solid

n = 2;
Sigma_w = eye(n);
Sigma = [1 1; 1 4];
S = 2 * [0 1; -1 0];

% Grid for plotting vector field
npts = 11;
[x1, x2] = meshgrid(linspace(-1, 1, npts), linspace(-1, 1, npts));
X = [x1(:)'; x2(:)'];

% High-resolution grid for smooth isocontours
[xq, yq] = meshgrid(linspace(-1.2, 1.2, 100), linspace(-1.2, 1.2, 100));

% --- FIGURE 1: Original Frame (x) ---
Sigma_inv = inv(Sigma);
D = -0.5 * Sigma_w * Sigma_inv;
C = S * Sigma_inv;
A = D + C;

f1 = figure('Color', 'w', 'Position', [100 100 1200 400]);
t1 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot 1: A
nexttile(1); 
plotvec(A * X, X, 'k', fs_label, fs_tick, 'x');
hold on;
% Energy landscape U(x) for isocontours
U_x = 0.5 * (Sigma_inv(1,1)*xq.^2 + 2*Sigma_inv(1,2)*xq.*yq + Sigma_inv(2,2)*yq.^2);
contour(xq, yq, U_x, contour_levels, contour_style, 'LineWidth', 1.5);
title('$A$', 'Interpreter','latex', 'FontSize', fs_title);
annotation(f1, 'textbox', [eq_x eq_y 0.05 0.1], 'String', '$=$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');

% Plot 2: Dissipative Component
nexttile(2); 
plotvec(D * X, X, 'b', fs_label, fs_tick, 'x');
title('$-\frac{1}{2}\Sigma_w\Sigma^{-1}$', 'Interpreter','latex', 'FontSize', fs_title);
annotation(f1, 'textbox', [pl_x eq_y 0.05 0.1], 'String', '$+$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');
ylabel([]);

% Plot 3: Conservative Component
nexttile(3); 
plotvec(C * X, X, 'r', fs_label, fs_tick, 'x');
title('$S\Sigma^{-1}$', 'Interpreter','latex', 'FontSize', fs_title);
ylabel([]);

% --- FIGURE 2: Whitened Frame (z) ---
S_half_inv = inv(sqrtm(Sigma));
Atilde = S_half_inv * A * inv(S_half_inv);
Dtilde = -0.5 * Sigma_inv; % Equivalently Dtilde = -0.5 * Sigma_w * Sigma_inv because Sigma_w = I
Stilde = S_half_inv * S * S_half_inv;

f2 = figure('Color', 'w', 'Position', [150 150 1200 400]);
t2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot 1: Atilde
nexttile(1); 
plotvec(Atilde * X, X, 'k', fs_label, fs_tick, 'z');
hold on;
% Isotropic energy landscape U(z) for isocontours
U_z = 0.5 * (xq.^2 + yq.^2);
contour(xq, yq, U_z, contour_levels, contour_style, 'LineWidth', 1.5);
title('$\tilde{A}$', 'Interpreter','latex', 'FontSize', fs_title);
annotation(f2, 'textbox', [eq_x eq_y 0.05 0.1], 'String', '$=$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');

% Plot 2: Whitened Dissipative Component
nexttile(2); 
plotvec(Dtilde * X, X, 'b', fs_label, fs_tick, 'z');
title('$-\frac{1}{2}\Sigma_w\Sigma^{-1}$', 'Interpreter','latex', 'FontSize', fs_title);
annotation(f2, 'textbox', [pl_x eq_y 0.05 0.1], 'String', '$+$', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', fs_annot, 'HorizontalAlignment', 'center');
ylabel([]);

% Plot 3: Whitened Conservative Component
nexttile(3); 
plotvec(Stilde * X, X, 'r', fs_label, fs_tick, 'z');
title('$\Sigma^{-1/2}S\Sigma^{-1/2}$', 'Interpreter','latex', 'FontSize', fs_title);
ylabel([]);

% --- Helper Function ---
function plotvec(vec, X, color, fs_label, fs_tick, var_name)
    quiver(X(1,:), X(2,:), vec(1,:), vec(2,:), color);
    
    % Dynamically set labels to x_1, x_2 or z_1, z_2 based on 'var_name'
    xlabel(['$' var_name '_1$'], 'Interpreter', 'latex', 'FontSize', fs_label);
    ylabel(['$' var_name '_2$'], 'Interpreter', 'latex', 'FontSize', fs_label);
    
    axis equal
    xlim([-1 1]); ylim([-1 1]); grid on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs_tick);
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
end