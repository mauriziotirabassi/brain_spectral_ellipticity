clear; clc; close all;

% --- 1. System Parameters ---
Sigma_w = eye(2);
Sigma   = diag([1, 10]); 
d = diag(inv(Sigma));    
alpha = (1 .* d) / 2; 
w_crit = abs(alpha(1) - alpha(2)) / (2 * sqrt(d(1)*d(2)));

% 5 Test Cases
omegas_cases = [0.05, 1.0, 1.3, 1.7, 3.2] * w_crit;

% Colors
cmap = magma(256);
dark_color = cmap(15, :); 

% --- 2. Figure Setup ---
f = figure('Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9], 'Color', 'w');

% --- 3. Main Phase Diagram (Bottom) ---
ax_main = axes('Position', [0.05, 0.08, 0.90, 0.22]); 
hold(ax_main, 'on'); grid(ax_main, 'on');

% Calculate Curve
omega_scan = linspace(0, 3.5*w_crit, 1000);
kappa_scan = nan(size(omega_scan));
for k = 1:length(omega_scan)
    [k_val, ~, ~] = calc_dynamics(omega_scan(k), Sigma, Sigma_w);
    kappa_scan(k) = k_val;
end

% Plot
plot(ax_main, omega_scan, kappa_scan, '-', 'Color', dark_color, 'LineWidth', 2);
xline(ax_main, w_crit, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);

% RESTORED: Critical Label
text(ax_main, w_crit, 1.5, '$\tilde\omega_\text{crit}$', ...
    'Interpreter','latex', 'FontSize', 14, 'VerticalAlignment','bottom');

% Axis Settings
set(ax_main, 'YScale', 'log');
xlim(ax_main, [0, 3.5*w_crit]);
ylim(ax_main, [1, 1e5]); 

xlabel(ax_main, '$\tilde\omega$', 'Interpreter','latex', 'FontSize', 16);
ylabel(ax_main, '$\kappa$', 'Interpreter','latex', 'FontSize', 16);

% --- 4. Panel Layout Logic ---
n_cols = 5;
margin_x = 0.02; % Minimized margins
spacing_x = 0.01; % Tight spacing
total_w = 1 - 2*margin_x;
p_width = (total_w - (n_cols-1)*spacing_x) / n_cols;

% Calculate Height for Square Aspect Ratio
fig_ar = f.Position(3) / f.Position(4); 
p_height = p_width * fig_ar; 

y_traj = 0.38;               % Shifted up to give space for arrows
y_mat  = y_traj + p_height + 0.04; 

x_starts = margin_x + (0:n_cols-1)*(p_width + spacing_x);

for i = 1:n_cols
    w_val = omegas_cases(i);
    [kappa, u1, u2] = calc_dynamics(w_val, Sigma, Sigma_w);
    
    % --- A. Matrix Plot (Top) ---
    ax_mat = axes('Position', [x_starts(i), y_mat, p_width, p_height]);
    
    U = [u1; u2];
    X_cos = U ./ (vecnorm(U, 2, 1) + 1e-9); 
    G_cos = X_cos' * X_cos;
    G_cos(tril(true(size(G_cos)), -1)) = NaN; 
    
    imagesc(ax_mat, G_cos, 'AlphaData', ~isnan(G_cos));
    colormap(ax_mat, magma);
    caxis(ax_mat, [-1 1]);
    axis(ax_mat, 'square'); 
    set(ax_mat, 'XTick', [], 'YTick', []); 
    
    % Labels ONLY on First Panel
    if i == 1
        xlabel(ax_mat, '$\tau_k$', 'Interpreter','latex', 'FontSize', 14);
        set(ax_mat, 'XAxisLocation', 'top');
        ylabel(ax_mat, '$\tau_\ell$', 'Interpreter','latex', 'FontSize', 14);
        set(ax_mat, 'YAxisLocation', 'right'); 
    end
    
    % --- B. Trajectory Plot (Middle) ---
    ax_traj = axes('Position', [x_starts(i), y_traj, p_width, p_height]);
    hold(ax_traj, 'on'); box(ax_traj, 'on'); grid(ax_traj, 'on');
    
    plot(ax_traj, u1, u2, '-', 'Color', dark_color, 'LineWidth', 1.5);
    plot(ax_traj, -u1, u2, '--', 'Color', dark_color, 'LineWidth', 1);
    
    % Force Square Limits
    limit = 1.1 * max(max(abs(u1)), max(abs(u2)));
    if limit == 0, limit = 1; end
    xlim(ax_traj, [-limit, limit]);
    ylim(ax_traj, [-limit, limit]);
    axis(ax_traj, 'square'); 
    
    % Label Logic: Remove numbers, keep titles on first panel
    set(ax_traj, 'XTickLabel', [], 'YTickLabel', []); 
    if i == 1
        xlabel(ax_traj, '$u_1$', 'Interpreter','latex', 'FontSize', 14);
        ylabel(ax_traj, '$u_2$', 'Interpreter','latex', 'FontSize', 14);
    end
    
    % --- C. Arrow Logic ---
    main_pos = get(ax_main, 'Position');
    
    % X Coordinate
    x_rel = (w_val - 0) / (3.5*w_crit - 0);
    x_arrow_start = main_pos(1) + x_rel * main_pos(3);
    
    % Y Coordinate (Log Scale Clamped)
    y_lims = get(ax_main, 'YLim');
    y_log_min = log10(y_lims(1)); 
    y_log_max = log10(y_lims(2));
    k_clamped = min(max(kappa, y_lims(1)), y_lims(2));
    y_val_log = log10(k_clamped);
    y_rel = (y_val_log - y_log_min) / (y_log_max - y_log_min);
    y_arrow_start = main_pos(2) + y_rel * main_pos(4);
    
    % Target Point: Stop arrow nicely below the plot (margin preserved)
    target_x = x_starts(i) + p_width/2;
    
    % Adjusted target to preserve margin and not hit labels
    % The labels are usually ~0.04 units high. We stop the arrow below them.
    target_y = y_traj - 0.06; 
    
    annotation('arrow', [x_arrow_start, target_x], [y_arrow_start, target_y], ...
        'Color', dark_color, 'LineWidth', 1, 'HeadStyle', 'vback2');
end

% --- 5. Colorbar ---
c = colorbar(ax_mat); % Attach to last axis
c.Label.String = 'Cross-Lag Similarity';
% Position manually to right of last matrix
c.Position = [x_starts(end) + p_width + 0.005, y_mat, 0.012, p_height];

% --- Helper Function ---
function [kappa, u1, u2] = calc_dynamics(omega, Sigma, Sigma_w)
    S = omega * [0 1; -1 0];
    A = (-0.5 * Sigma_w + S) / Sigma; 
    S_half = sqrtm(Sigma);             
    Atilde = inv(S_half) * A * S_half;
    
    mu = trace(Atilde)/2;                  
    Delta = trace(Atilde)^2 - 4*det(Atilde);
    
    d = diag(inv(Sigma));
    lags = linspace(0, 45, 300);
    
    if abs(Delta) < 1e-5 % Critical
        u1 = sqrt(2) * ones(size(lags));
        limit_slope = omega * (d(1) + d(2)); 
        u2 = limit_slope * lags;
        kappa = 1e5; 
    else
        gamma = sqrt(abs(Delta))/2;
        J = (Atilde - mu*eye(2)) / gamma;
        kappa = trace(J'*J); 
        
        if Delta > 0 % Overdamped
            C = cosh(gamma * lags);
            S = sinh(gamma * lags);
        else % Oscillatory
            C = cos(gamma * lags);
            S = sin(gamma * lags);
        end
        u1 = sqrt(2) * C; 
        u2 = sqrt(kappa) * S;
    end
end