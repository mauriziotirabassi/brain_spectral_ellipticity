clear; clc;
rng(42);

% --- Control Panel ---
% Parameters 
alpha1  = 1;     
alpha2  = 1;     
gamma   = 0;     
omega   = 1.5;   

% Visualization Settings
tick_font_size    = 15;   
label_font_size   = 23;    
sigma_font_adder  = 1;    % Knob: Manual adder for Sigma labels
tau_label_y_knob  = -1; 
ab_y_pos_knob     = 1;  
sigma_x_pos_knob  = -0.08; % Knob: Horizontal position of Sigma labels (Normalized)

color_tau_l       = [1, 0, 0];    
color_tau_k       = [0, 0, 1];    
color_grid        = [0.8 0.8 0.8]; 

% Case Definition
A = [-alpha1, gamma - omega; gamma + omega, -alpha2];

% Computation
lastLag  = 15;
numLags  = 500;
lags     = linspace(0, lastLag, numLags + 1);

% Get the darkest color from the magma colormap
magma_map = magma(256);
dark_magma = magma_map(1, :);

% Select Lags
idx_l = 100; 
idx_k = 250; 
val_tau_l = lags(idx_l);
val_tau_k = lags(idx_k);

% Compute Covariance Evolution
Sigma_tau = @(tau) expm(A * tau);
Corr = nan(2, 2, numel(lags));
for k = 1:numel(lags)
    Corr(:,:,k) = Sigma_tau(lags(k));
end

X = reshape(Corr, [], size(Corr, 3));
tlc_labels = ["$\tilde{\Sigma}_{11}$", "$\tilde{\Sigma}_{21}$", "$\tilde{\Sigma}_{12}$", "$\tilde{\Sigma}_{22}$"];

X_cos = X ./ vecnorm(X, 2, 1);
G_cos = X_cos' * X_cos;
G_cos(tril(true(size(G_cos)), -1)) = NaN;

% --- Plotting ---
figure('Color', 'w', 'Position', [100 100 1100 600]);

t_outer = tiledlayout(1, 5, 'TileSpacing', 'loose', 'Padding', 'compact');

t_left = tiledlayout(t_outer, 4, 1);
t_left.Layout.Tile = 1; 
t_left.Layout.TileSpan = [1, 2]; 
t_left.TileSpacing = 'none';
t_left.Padding = 'none';

% --- LEFT COLUMN: Stacked Plots ---
for i = 1:4
    ax = nexttile(t_left);
    hold on;
    
    % Manual Grid
    x_grid = 0:5:15;
    for xg = x_grid
        xline(xg, '-', 'Color', color_grid, 'LineWidth', 0.5);
    end
    yline(0, '-', 'Color', color_grid, 'LineWidth', 0.5); 
    
    xline(val_tau_l, '--', 'Color', color_tau_l, 'LineWidth', 1.5);
    xline(val_tau_k, '--', 'Color', color_tau_k, 'LineWidth', 1.5);
    
    plot(lags, X_cos(i, :), 'LineWidth', 1.2, 'Color', dark_magma);
    
    % --- MANUAL SIGMA LABELS (Bypassing ylabel threshold) ---
    text(sigma_x_pos_knob, 0.5, tlc_labels(i), 'Units', 'normalized', ...
         'Interpreter', 'latex', 'FontSize', label_font_size + sigma_font_adder, ...
         'HorizontalAlignment', 'right', 'Clipping', 'off');
       
    xlim([0, lastLag]);
    y_min_padded = min(X_cos(i, :)) - 0.25; 
    y_max_padded = max(X_cos(i, :)) + 0.25;
    ylim([y_min_padded, y_max_padded]);
    box off;
    set(gca, 'FontSize', tick_font_size);
    
    if i < 4
        set(gca, 'XColor', 'none'); 
    else
        xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', label_font_size);
        set(gca, 'XColor', 'k'); 
        
        text(val_tau_l, tau_label_y_knob, '$\tau_\ell$', 'Interpreter', 'latex', 'Color', color_tau_l, ...
             'FontSize', label_font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        text(val_tau_k, tau_label_y_knob, '$\tau_k$', 'Interpreter', 'latex', 'Color', color_tau_k, ...
             'FontSize', label_font_size, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
    
    if i == 1
        text(-0.15, ab_y_pos_knob, '\textbf{A}', 'Units', 'normalized', 'Interpreter', 'latex', ...
             'FontSize', label_font_size+4, 'Clipping', 'off');
    end
end

% --- RIGHT COLUMN: CLS Matrix ---
ax_right = nexttile(t_outer, 3, [1, 3]);

imagesc(lags, lags, G_cos, 'AlphaData', ~isnan(G_cos));
axis square;
colormap(magma);
c = colorbar;
c.Label.Interpreter = 'latex';
clim([-1, 1]);

hold on; 

gap = 0.25; 
line([val_tau_k, val_tau_k], [0, val_tau_l - gap], 'Color', color_tau_k, 'LineStyle', '--', 'LineWidth', 1.5);
line([val_tau_k, val_tau_k], [val_tau_l + gap, lastLag], 'Color', color_tau_k, 'LineStyle', '--', 'LineWidth', 1.5);
line([0, val_tau_k - gap], [val_tau_l, val_tau_l], 'Color', color_tau_l, 'LineStyle', '--', 'LineWidth', 1.5);
line([val_tau_k + gap, lastLag], [val_tau_l, val_tau_l], 'Color', color_tau_l, 'LineStyle', '--', 'LineWidth', 1.5);

plot(val_tau_k, val_tau_l, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'MarkerSize', 10);

set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'YDir', 'reverse');
set(gca, 'FontSize', tick_font_size);

box on; 
set(gca, 'Layer', 'top'); 

xlabel('$\tau_k$', 'Interpreter', 'latex', 'FontSize', label_font_size); 
ylabel('$\tau_\ell$', 'Interpreter', 'latex', 'FontSize', label_font_size);

text(-0.15, ab_y_pos_knob, '\textbf{B}', 'Units', 'normalized', 'Interpreter', 'latex', ...
     'FontSize', label_font_size+4, 'Clipping', 'off');