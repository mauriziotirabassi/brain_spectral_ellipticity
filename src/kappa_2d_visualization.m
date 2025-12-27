clear; clc; close all;

% --- 1. System Parameters ---
Sigma_w = eye(2);
Sigma   = diag([1, 100]); 
d = diag(inv(Sigma));    
sigma_w_scalar = 1;      
alpha = (sigma_w_scalar .* d) / 2; 

% Critical Frequency Calculation
w_crit = abs(alpha(1) - alpha(2)) / (2 * sqrt(d(1)*d(2)));

% Define the 3 Test Cases
omegas_cases = [0.2 * w_crit, w_crit, 2.5 * w_crit];
titles = {'Overdamped', 'Critical', 'Oscillatory'};
colors = {[0.85, 0.33, 0.1], [0.47, 0.67, 0.19], [0, 0.45, 0.74]}; 

% --- 2. Main Canvas Setup ---
f = figure('Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9], 'Color', 'w');

% Define Main Axis Position (Bottom Panel)
% [left bottom width height]
ax_main = axes('Position', [0.08, 0.1, 0.84, 0.30]); 
hold(ax_main, 'on'); grid(ax_main, 'on');

% --- 3. Plot Main Phase Diagram (Kappa vs Omega) ---
numPts = 1000;
omega_scan = linspace(0, 3*w_crit, numPts);
kappa_scan = nan(size(omega_scan));

for k = 1:numPts
    w_val = omega_scan(k);
    [kappa_val, ~, ~] = calc_dynamics(w_val, Sigma, Sigma_w);
    kappa_scan(k) = kappa_val;
end

plot(ax_main, omega_scan, kappa_scan, 'k-', 'LineWidth', 2);
set(ax_main, 'YScale', 'log');
xlim(ax_main, [0, 3*w_crit]);
ylim(ax_main, [1, max(kappa_scan)*2]);
xlabel(ax_main, 'Frequency \omega', 'FontSize', 12);
ylabel(ax_main, 'Dynamic Anisotropy \kappa(\omega)', 'FontSize', 12);
% title(ax_main, 'Phase Diagram: System Regimes', 'FontSize', 14);

% Add Critical Line
xline(ax_main, w_crit, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
text(ax_main, w_crit, 1.5, ' \omega_{crit}', 'FontSize', 12, 'VerticalAlignment','bottom');

% --- 4. Loop to Create the "Pop-up" Panels ---
panel_y = 0.5; 
panel_w = 0.1290; 
panel_h = 0.2070;
gap = 0.005;      

% Centers adjusted: Start at 0.135 (prevents neg left), end at 0.77 (saves right margin)
group_x_centers = [0.19, 0.48, 0.77]; 

% Store handles to link colorbar later
ax_matrices = [];
last_matrix_pos = []; 

for i = 1:3
    w_case = omegas_cases(i);
    col = colors{i};
    
    % A. Calculate Dynamics
    [kappa, u1, u2] = calc_dynamics(w_case, Sigma, Sigma_w);
    
    % B. Calculate Similarity Matrix
    U = [u1; u2];
    X_cos = U ./ (vecnorm(U, 2, 1) + 1e-9); 
    G_cos = X_cos' * X_cos;
    
    % Mask Lower Triangle
    G_cos(tril(true(size(G_cos)), -1)) = NaN;
    lags = linspace(0, 25, length(u1));

    % --- DRAW ARROWS ---
    ax_pos = get(ax_main, 'Position');
    x_limits = get(ax_main, 'XLim');
    y_limits = get(ax_main, 'YLim');
    
    x_norm = ax_pos(1) + ax_pos(3) * (w_case - x_limits(1)) / (x_limits(2) - x_limits(1));
    log_y = log10(kappa);
    log_ylim = log10(y_limits);
    y_norm = ax_pos(2) + ax_pos(4) * (log_y - log_ylim(1)) / (log_ylim(2) - log_ylim(1));
    
    target_x = group_x_centers(i); 
    target_y = panel_y - 0.03; 
    
    annotation('arrow', [x_norm, target_x], [y_norm, target_y], ...
        'Color', col, 'LineWidth', 1.5, 'HeadStyle', 'vback2');
    
    % --- PLOT 1: GEOMETRY ---
    pos_geom = [group_x_centers(i) - panel_w - gap/2, panel_y, panel_w, panel_h];
    ax_geom = axes('Position', pos_geom);
    box(ax_geom, 'on'); hold(ax_geom, 'on'); grid(ax_geom, 'on');
    
    plot(ax_geom, u1, u2, 'Color', col, 'LineWidth', 2);
    plot(ax_geom, -u1, u2, '--', 'Color', col, 'LineWidth', 1.5);
    
    axis(ax_geom, 'equal'); 
    xlabel(ax_geom, 'u_1', 'Interpreter', 'tex', 'FontSize', 9);
    
    % FIX 1: Only show u2 label on the first plot
    if i == 1
        ylabel(ax_geom, 'u_2', 'Interpreter', 'tex', 'FontSize', 9);
    end
    
    set(ax_geom, 'FontSize', 8);
    title(ax_geom, titles{i}, 'Color', col, 'FontWeight', 'bold', 'FontSize', 14);

    % --- PLOT 2: MATRIX ---
    pos_mat = [group_x_centers(i) + gap/2, panel_y, panel_w, panel_h];
    ax_mat = axes('Position', pos_mat);
    ax_matrices = [ax_matrices, ax_mat]; 
    last_matrix_pos = pos_mat; 
    
    h = imagesc(ax_mat, lags, lags, G_cos);
    set(h, 'AlphaData', ~isnan(G_cos)); 
    
    axis(ax_mat, 'square');
    colormap(ax_mat, magma);
    clim(ax_mat, [-1 1]);
    
    set(ax_mat, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
    xlabel(ax_mat, '\tau_1', 'FontSize', 9);
    
    % FIX 2: Only show tau2 label on the last plot (Oscillatory)
    if i == 3
        ylabel(ax_mat, '\tau_2', 'FontSize', 9);
    end

    % --- DRAW BOUNDING BOX (With Safety Clamps) ---
    % Calculate exact coordinates
    rect_x = pos_geom(1) - 0.012; 
    rect_y = pos_geom(2) - 0.02;
    rect_w = (pos_mat(1) + pos_mat(3)) - pos_geom(1) + 0.024;
    rect_h = pos_geom(4) + 0.06;
    
    % CRITICAL FIX: Ensure values stay within [0, 1]
    rect_x = max(0, rect_x);
    if (rect_x + rect_w) > 1, rect_w = 1 - rect_x; end
    
    annotation('rectangle', [rect_x, rect_y, rect_w, rect_h], ...
        'Color', [0.8 0.8 0.8], 'LineStyle', ':');
end

% --- 5. Global Colorbar ---
c = colorbar(ax_matrices(end)); 
c.Label.String = 'Cosine Similarity';
c.Label.FontSize = 12;

% Position Calculation:
% Place it further right.
cb_x_start = last_matrix_pos(1) + last_matrix_pos(3) + 0.055; 
cb_width = 0.012; 

% Ensure colorbar doesn't go off screen (safety check)
if cb_x_start + cb_width > 0.98
   cb_x_start = 0.965; % Hard cap if screen aspect ratio is weird
end

c.Position = [cb_x_start - 0.02, panel_y, cb_width, panel_h]; 

% --- Helper Function ---
function [kappa, u1, u2] = calc_dynamics(omega, Sigma, Sigma_w)
    S = omega * [0 1; -1 0];
    A = (-0.5 * Sigma_w + S) / Sigma; 
    S_half = sqrtm(Sigma);             
    Atilde = inv(S_half) * A * S_half;
    
    mu = trace(Atilde)/2;                  
    Delta = trace(Atilde)^2 - 4*det(Atilde);
    
    d = diag(inv(Sigma));
    lags = linspace(0, 25, 300);

    if abs(Delta) < 1e-5 % Critical
        u1 = sqrt(2) * ones(size(lags));
        limit_slope = omega * (d(1) + d(2)); 
        u2 = limit_slope * lags;
        kappa = 1e4; 
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