clear; clc;
rng(42);

% --- Control Panel ---
% Parameters 
alpha1  = 1;     % Dissipative rate 1
alpha2  = 1.2;   % Dissipative rate 2
gamma   = 1;     % Shear
omega   = 1.5;   % Rotational frequency

% Font Size Control
tick_font_size  = 13;   % Size for axis ticks (numbers)
label_font_size = 19;   % Size for labels (Sigma, tau, etc.)

% Label Positioning Knobs (Normalized Figure Coordinates 0 to 1)
% Y-Labels (Function Names)
y_label_x     = -0.08;   % Horizontal position (Right-Left displacement)
y_label_start = 0.88;    % Vertical position of the top label
y_label_gap   = 0.255;   % Vertical distance between labels
% X-Label (Lag tau)
x_label_x     = 0.5;     % Horizontal position
x_label_y     = -0.035;  % Vertical position

% Case 1: Trivial Commutativity (omega = 0)
A = [-alpha1, gamma - omega; gamma + omega, -alpha2];
lastLag  = 15;
numLags  = 500;
lags     = linspace(0, lastLag, numLags + 1);

% Screen size constraint (<14 inches)
figure('Color', 'w', 'Position', [100 100 1100 600]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Labels for Stacked Plot (LaTeX)
tlc_labels = ["$\tilde{\Sigma}_{11}$", "$\tilde{\Sigma}_{21}$", "$\tilde{\Sigma}_{12}$", "$\tilde{\Sigma}_{22}$"];
    
% Compute Covariance Evolution
Sigma_tau = @(tau) expm(A * tau);
Corr = nan(2, 2, numel(lags));
for k = 1:numel(lags)
    Corr(:,:,k) = Sigma_tau(lags(k));
end

% Reshape (Order: 11, 21, 12, 22)
X = reshape(Corr, [], size(Corr, 3));

% Compute CLS
% Normalize columns (vectors at each time point)
X_cos = X ./ vecnorm(X, 2, 1);

% Gram Matrix
G_cos = X_cos' * X_cos;
% Mask lower triangle
G_cos(tril(true(size(G_cos)), -1)) = NaN;

% Plot Left: TLC Functions (Stacked Plot)
nexttile
s = stackedplot(lags, X_cos.', 'LineWidth', 1.2, 'GridVisible', 'on');
s.DisplayLabels = repmat("", 1, 4); % Set empty strings to hide native labels
s.XLabel = '';                      % Remove default xlabel
s.FontSize = tick_font_size;        % Set tick font size

% Plot Right: CLS Matrix
nexttile
imagesc(lags, lags, G_cos, 'AlphaData', ~isnan(G_cos));
axis square;
colormap(magma);
colorbar;
clim([-1, 1]);
xlabel('$\tau_k$', 'Interpreter', 'latex', 'FontSize', 50);
ylabel('$\tau_\ell$', 'Interpreter', 'latex', 'FontSize', 50);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'YDir', 'reverse');
set(gca, 'FontSize', tick_font_size); % Set tick font size

% --- Manual Label Overlays ---
% Create a transparent axes over the entire figure for absolute positioning
% We use the tiledlayout 't' as parent to ensure alignment
overlay_ax = axes(t, 'Visible', 'off', 'Position', [0 0 1 1]);
xlim(overlay_ax, [0 1]);
ylim(overlay_ax, [0 1]);

% Loop to place Function Name labels
for i = 1:4
    % Calculate vertical position based on start and gap
    y_pos = y_label_start - (i-1) * y_label_gap;
    
    text(overlay_ax, y_label_x, y_pos, tlc_labels(i), ...
        'Interpreter', 'latex', 'FontSize', label_font_size, ...
        'HorizontalAlignment', 'center'); 
end

% Place Single X-Label for the stacked plot
text(overlay_ax, x_label_x, x_label_y, '$\tau$', ...
    'Interpreter', 'latex', 'FontSize', label_font_size, ...
    'HorizontalAlignment', 'center');