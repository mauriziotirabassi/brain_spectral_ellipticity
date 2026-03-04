clear; clc; close all;

% =========================================================================
% --- USER CONFIGURATION ---
% =========================================================================

% 1. Simulation Parameters
% Anisotropy levels (delta) to compare
deltas = [0.2, 0.5, 1.0, 1.5]; 

% Frequency range (omega) for the scan
omega_max = 3.0;
omega_step = 0.001;
omega_scan = 0:omega_step:omega_max;

% 2. Visualization
FIG_W_PX = 1100;
FIG_H_PX = 600;
f_size_lbl = 18;
f_size_ax  = 14;
line_width = 2.0;

% Colors (Magma colormap for distinct levels)
cmap = magma(length(deltas) + 2); 
colors = cmap(1:length(deltas), :);

% =========================================================================
% --- CALCULATION & PLOTTING ---
% =========================================================================

figure('Color', 'w', 'Position', [100 100 FIG_W_PX FIG_H_PX]);
ax = axes('Position', [0.12, 0.15, 0.8, 0.75]);
hold(ax, 'on'); grid(ax, 'on');

% Loop through each anisotropy level
for i = 1:length(deltas)
    d_val = deltas(i);
    
    % Preallocate
    kappa_vals = zeros(size(omega_scan));
    
    % Calculate Kappa based on the user-provided snippet logic
    % Formula: kappa = (omega^2 + delta^2) / abs(omega^2 - delta^2)
    for k = 1:length(omega_scan)
        w_val = omega_scan(k);
        
        num = w_val^2 + d_val^2;
        den = abs(w_val^2 - d_val^2);
        
        % Handle numerical singularity
        if den < 1e-9 * num
            kappa_vals(k) = NaN; % Break line at singularity
        else
            kappa_vals(k) = num / den;
        end
    end
    
    % Plot Curve
    plot(ax, omega_scan, kappa_vals, '-', ...
        'Color', colors(i,:), ...
        'LineWidth', line_width, ...
        'DisplayName', sprintf('$\\delta = %.1f$', d_val));
    
    % Mark the singularity (Exceptional Point)
    xline(ax, d_val, '--', ...
        'Color', [colors(i,:) 0.4], ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off');
end

% =========================================================================
% --- FORMATTING ---
% =========================================================================

set(ax, 'YScale', 'log');
set(ax, 'FontSize', f_size_ax);
xlim(ax, [0, omega_max]);
ylim(ax, [1, 1e4]); % Clamp y-axis to keep plot readable

xlabel(ax, 'Frequency Parameter $\omega$', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
ylabel(ax, 'Non-Normality Factor $\kappa$', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
title(ax, '$\kappa$ Divergence at Exceptional Points ($\omega = \delta$)', 'Interpreter', 'latex', 'FontSize', f_size_lbl);

% Legend
lgd = legend(ax, 'Location', 'northeast');
set(lgd, 'Interpreter', 'latex', 'FontSize', f_size_ax);

box(ax, 'on');