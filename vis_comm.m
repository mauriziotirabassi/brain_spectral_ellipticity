clear; clc;
rng(42);

% --- Control Panel ---
% Parameters 
alpha1  = 1;     % Dissipative rate 1
alpha2  = 1.2;   % Dissipative rate 2
gamma   = 0;     % Shear
omega   = 1.5;   % Rotational frequency

% Font Size Control
tick_font_size  = 13;   
label_font_size = 19;   

% Case Definition (Using alpha2 for general case)
A = [-alpha1, gamma - omega; gamma + omega, -alpha2];

% Computation
lastLag  = 15;
numLags  = 500;
lags     = linspace(0, lastLag, numLags + 1);

% Compute Covariance Evolution
Sigma_tau = @(tau) expm(A * tau);
Corr = nan(2, 2, numel(lags));
for k = 1:numel(lags)
    Corr(:,:,k) = Sigma_tau(lags(k));
end

% Reshape (Order: 11, 21, 12, 22)
X = reshape(Corr, [], size(Corr, 3));
tlc_labels = ["$\tilde{\Sigma}_{11}$", "$\tilde{\Sigma}_{21}$", "$\tilde{\Sigma}_{12}$", "$\tilde{\Sigma}_{22}$"];

% Compute CLS
X_cos = X ./ vecnorm(X, 2, 1);
G_cos = X_cos' * X_cos;
G_cos(tril(true(size(G_cos)), -1)) = NaN;

% --- Plotting ---
figure('Color', 'w', 'Position', [100 100 1100 600]);

% Main Layout: 4 rows, 2 columns
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Left Column: Stacked Plots
for i = 1:4
    % Select odd tiles (1, 3, 5, 7)
    nexttile(2*i - 1);
    
    plot(lags, X_cos(i, :), 'k-', 'LineWidth', 1.2);
    grid on;
    
    % Y-Label (Function Name)
    ylabel(tlc_labels(i), 'Interpreter', 'latex', 'FontSize', label_font_size, ...
           'Rotation', 0, 'HorizontalAlignment', 'right');
       
    % Axis formatting
    set(gca, 'FontSize', tick_font_size);
    xlim([0, lastLag]);
    ylim([-1, 1]);
    
    % Handle X-Labels: Only show on the bottom plot
    if i < 4
        xticklabels({});
    else
        xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', label_font_size);
    end
end

% 2. Right Column: CLS Matrix
% Span all 4 rows in the 2nd column
nexttile(2, [4, 1]);

imagesc(lags, lags, G_cos, 'AlphaData', ~isnan(G_cos));
axis square;
colormap(magma);
c = colorbar;
clim([-1, 1]);

% Labels
xlabel('$\tau_k$', 'Interpreter', 'latex', 'FontSize', label_font_size + 4); % Slightly larger for emphasis
ylabel('$\tau_\ell$', 'Interpreter', 'latex', 'FontSize', label_font_size + 4);

% Axis formatting
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'YDir', 'reverse');
set(gca, 'FontSize', tick_font_size);