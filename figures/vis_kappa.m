clear; clc; close all;

% User Configuration

% Manual Alignment Knobs
TOP_OFFSET    = -0.02;  % Positive value shifts the top rows right, negative shifts left
MAIN_X_START  = 0.08;   % Starting X coordinate of the bottom main plot (normalized 0 to 1)
MAIN_WIDTH    = 0.88;   % Width of the bottom main plot (normalized 0 to 1)
CBAR_FONT_SIZE= 20;     % Font size for the Colorbar label
CBAR_GAP      = 0.005;  % Horizontal gap between the last matrix panel and the colorbar

% Panel Layout
PANEL_W      = 0.175;   % Width of individual top and middle panels
GAP_X        = 1e-14;   % Horizontal gap between individual panels
GAP_Y        = 0.02;    % Vertical gap between the top row (matrix) and middle row (trajectory)
Y_MID        = 0.36;    % Starting Y coordinate of the middle row (trajectory plots)

% Label Positioning Tweaks (Normalized Units relative to Axis)
X_LBL_Y_SHIFT = -0.15;  % Vertical shift for the bottom plot X-axis label
TAU_L_SHIFT   = -0.02;  % Horizontal shift for the left Y-axis label (tau_l)
TAU_K_Y_SHIFT = 0.01;   % Vertical shift for the top X-axis label (tau_k)
U1_LBL_SHIFT  = -0.01;  % Vertical shift for the middle plot X-axis label (u1)

% Font Sizes
f_size_lbl      = 19;   % Font size for general labels
f_size_ax       = 14;   % Font size for axis ticks
label_font_size = 20;   % Font size for main axis labels

% Figure Size
FIG_W_PX = 1400;        % Target figure width in pixels for aspect ratio calculation
FIG_H_PX = 850;         % Target figure height in pixels for aspect ratio calculation

% System Parameters
Sigma_w = eye(2);
Sigma   = diag([1, 10]);
d = diag(inv(Sigma));
alpha = (1 .* d) / 2;
w_crit = abs(alpha(1) - alpha(2)) / (2 * sqrt(d(1)*d(2)));
omegas_cases = [0.05, 1.0, 1.3, 1.7, 3.2] * w_crit;

% Colors
cmap = magma(256);
dark_color = cmap(15, :);

% Figure Setup
f = figure('Color', 'w', 'Position', [100 100 1100 600]);

% Calculate Height
FIG_AR  = FIG_W_PX / FIG_H_PX;
PANEL_H = PANEL_W * FIG_AR;

% Main Phase Diagram (Bottom)
ax_main = axes('Position', [MAIN_X_START, 0.10, MAIN_WIDTH, 0.18]);
hold(ax_main, 'on'); grid(ax_main, 'on');

% Calculate Curve for main phase diagram
omega_scan = linspace(0, 3.5*w_crit, 1000);
kappa_scan = nan(size(omega_scan));
for k = 1:length(omega_scan)
    [k_val, ~, ~] = calc_dynamics(omega_scan(k), Sigma, Sigma_w);
    kappa_scan(k) = k_val;
end

% Plot main phase diagram
plot(ax_main, omega_scan, kappa_scan, '-', 'Color', dark_color, 'LineWidth', 2);
xline(ax_main, w_crit, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);

% Critical Label
text(ax_main, w_crit * 1.1, 0.14, '$\tilde{\omega}_{\mathrm{crit}}$', ...
    'Interpreter','latex', 'FontSize', f_size_lbl, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom');

% Axis Formatting
set(ax_main, 'FontSize', f_size_ax);
set(ax_main, 'YScale', 'log');
set(ax_main, 'YMinorGrid', 'on');
set(ax_main, 'YTick', 10.^(0:8));
xlim(ax_main, [0, 3.5*w_crit]);
ylim(ax_main, [1, 1e3]);

% Manual X-Label Positioning
lbl = xlabel(ax_main, '$\tilde{\omega}$ (Hz)', 'Interpreter','latex', 'FontSize', label_font_size);
ylabel(ax_main, '$\kappa$', 'Interpreter','latex', 'FontSize', label_font_size);
set(lbl, 'Units', 'normalized');
pos_lbl = get(lbl, 'Position');
pos_lbl(2) = X_LBL_Y_SHIFT;
set(lbl, 'Position', pos_lbl, 'VerticalAlignment', 'top');

% Panels Generation
n_cols = 5;
total_w = (n_cols * PANEL_W) + ((n_cols - 1) * GAP_X);
start_x = ((1 - total_w) / 2) + TOP_OFFSET;
Y_TOP   = Y_MID + PANEL_H + GAP_Y;

for i = 1:n_cols
    w_val = omegas_cases(i);
    [kappa, u1, u2] = calc_dynamics(w_val, Sigma, Sigma_w);
    
    pos_x = start_x + (i-1)*(PANEL_W + GAP_X);
    
    % A. Matrix Plot (Top) - Cross-Lag Similarity
    ax_mat = axes('Position', [pos_x, Y_TOP, PANEL_W, PANEL_H]);
    
    U = [u1; u2];
    X_cos = U ./ (vecnorm(U, 2, 1) + 1e-9);
    G_cos = X_cos' * X_cos;
    G_cos(tril(true(size(G_cos)), -1)) = NaN;
    
    imagesc(ax_mat, G_cos, 'AlphaData', ~isnan(G_cos));
    colormap(ax_mat, magma);
    clim(ax_mat, [-1 1]);
    
    axis(ax_mat, 'tight');
    axis(ax_mat, 'square');
    set(ax_mat, 'XTick', [], 'YTick', [], 'FontSize', f_size_ax);
    
    if i == 1
        % Top Label
        xt = xlabel(ax_mat, '$\tau_k$', 'Interpreter','latex', 'FontSize', label_font_size);
        set(ax_mat, 'XAxisLocation', 'top');
        
        set(xt, 'Units', 'normalized');
        pos_xt = get(xt, 'Position');
        pos_xt(2) = 1.0 + TAU_K_Y_SHIFT;
        set(xt, 'Position', pos_xt, 'VerticalAlignment', 'bottom');
        
        % Left Label
        yl = ylabel(ax_mat, '$\tau_\ell$', 'Interpreter','latex', 'FontSize', label_font_size);
        set(ax_mat, 'YAxisLocation', 'left');
        
        set(yl, 'Units', 'normalized');
        pos_yl = get(yl, 'Position');
        pos_yl(1) = TAU_L_SHIFT;
        set(yl, 'Position', pos_yl);
    end
    
    % B. Trajectory Plot (Middle) - Lag-Vector Trajectory
    ax_traj = axes('Position', [pos_x, Y_MID, PANEL_W, PANEL_H]);
    hold(ax_traj, 'on'); box(ax_traj, 'on'); grid(ax_traj, 'on');
    
    plot(ax_traj, u1, u2, '-', 'Color', dark_color, 'LineWidth', 2);
    plot(ax_traj, u1, -u2, '--', 'Color', dark_color, 'LineWidth', 2);
    
    limit = 1.1 * max(max(abs(u1)), max(abs(u2)));
    if limit == 0, limit = 1; end
    xlim(ax_traj, [-limit, limit]);
    ylim(ax_traj, [-limit, limit]);
    
    axis(ax_traj, 'square');
    set(ax_traj, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', f_size_ax);
    
    if i == 1
        % u1 label (X)
        xl_traj = xlabel(ax_traj, '$u_1$', 'Interpreter','latex', 'FontSize', label_font_size);
        set(xl_traj, 'Units', 'normalized');
        p_xl = get(xl_traj, 'Position');
        p_xl(2) = U1_LBL_SHIFT;
        set(xl_traj, 'Position', p_xl, 'VerticalAlignment', 'top');
        
        % u2 label (Y)
        yl_traj = ylabel(ax_traj, '$u_2$', 'Interpreter','latex', 'FontSize', label_font_size);
        set(yl_traj, 'Units', 'normalized');
        p_yl = get(yl_traj, 'Position');
        p_yl(1) = TAU_L_SHIFT;
        set(yl_traj, 'Position', p_yl);
    end
    
    % C. Arrow Logic
    main_pos = get(ax_main, 'Position');
    
    x_rel = (w_val - 0) / (3.5*w_crit - 0);
    x_arrow_start = main_pos(1) + x_rel * main_pos(3);
    
    y_lims = get(ax_main, 'YLim');
    y_log_min = log10(y_lims(1));
    y_log_max = log10(y_lims(2));
    k_clamped = min(max(kappa, y_lims(1)), y_lims(2));
    y_val_log = log10(k_clamped);
    y_rel = (y_val_log - y_log_min) / (y_log_max - y_log_min);
    y_arrow_start = main_pos(2) + y_rel * main_pos(4);
    
    target_x = pos_x + PANEL_W/2;
    target_y = Y_MID - 0.05;
    
    annotation('arrow', [x_arrow_start, target_x], [y_arrow_start, target_y], ...
        'Color', dark_color, 'LineWidth', 1, 'HeadStyle', 'vback2');
end

% Colorbar
c = colorbar(ax_mat);
c.Label.String = 'Cross-Lag Similarity';
c.Label.Interpreter = 'latex';
c.Label.FontSize = CBAR_FONT_SIZE;
c.TickLabelInterpreter = 'tex';
c.FontSize = f_size_ax;
c.Position = [start_x + total_w + CBAR_GAP, Y_TOP, 0.012, PANEL_H];

% Dynamics Function
function [kappa, u1, u2] = calc_dynamics(omega, Sigma, Sigma_w)
    S = omega * [0 1; -1 0];
    A = (-0.5 * Sigma_w + S) / Sigma;
    S_half = sqrtm(Sigma);
    Atilde = inv(S_half) * A * S_half;
    
    mu = trace(Atilde)/2;
    Delta = trace(Atilde)^2 - 4*det(Atilde);
    
    lags = linspace(0, 45, 2000);
    
    % Critical Regime
    if abs(Delta) < 1e-5
        h = (Atilde(1,1) - Atilde(2,2)) / 2;
        gamma_val = (Atilde(1,2) + Atilde(2,1)) / 2;
        delta = sqrt(h^2 + gamma_val^2);
        
        u1 = ones(size(lags));
        u2 = sqrt(2) * delta * lags;
        kappa = 1e5;
    % Oscillatory and Overdamped Regimes
    else
        nu = sqrt(abs(Delta))/2;
        A0 = Atilde - mu*eye(2);
        J = A0 / nu;
        kappa = trace(J'*J) / 2;
        
        if Delta > 0
            C = cosh(nu * lags);
            S_func = sinh(nu * lags);
        else
            C = cos(nu * lags);
            S_func = sin(nu * lags);
        end
        
        u1 = C;
        u2 = sqrt(kappa) * S_func;
    end
end