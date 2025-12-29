clear; clc; close all;

% --- GLOBAL PARAMETERS ---
n = 4; 
Sigma = eye(n); 
Sigma_w = eye(n);

% Time settings
lastLag = 25; 
numLags = 600; 
lags = linspace(0, lastLag, numLags + 1);

% Define the three spectral regimes
regimes = {
    'Harmonic',  [1, 2],       'Harmonic (1:2)';
    'Beat',      [3, 3.5],       'Modulation (Beats)';
    'Dissonant', [1, sqrt(5)],   'Dissonant (Irrational)'
};

% --- VISUALIZATION SETUP ---
figure('Color', 'w', 'Position', [100, 100, 1400, 450]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:3
    % 1. Construct S Matrix
    w1 = regimes{i, 2}(1);
    w2 = regimes{i, 2}(2);
    S = [0  w1 0  0; 
        -w1 0  0  0; 
         0  0  0  w2; 
         0  0 -w2 0];
    
    % 2. System Dynamics
    A = (-0.5 * Sigma_w + S) / Sigma;
    
    % 3. Compute Covariance & Similarity
    Cov_th = nan(n, n, numel(lags));
    for k = 1:numel(lags)
        Cov_th(:,:,k) = expm(A * lags(k)) * Sigma;
    end
    stds_th = sqrt(diag(Sigma));
    Corr_th = Cov_th ./ (stds_th * stds_th.');
    X = reshape(Corr_th, [], size(Corr_th, 3)); 
    X_cos = X ./ vecnorm(X, 2, 1); 
    G_cos = X_cos' * X_cos;
    
    % Mask lower triangle
    G_cos(tril(true(size(G_cos)), -1)) = NaN;
    
    % 4. Plotting
    nexttile;
    imagesc(lags, lags, G_cos, 'AlphaData', ~isnan(G_cos));
    axis square; 
    colormap(magma); 
    clim([-1, 1]);
    
    % Titles and X-Labels (Common)
    title(['{' regimes{i, 3} '}'], 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('$\tau_k$', 'Interpreter', 'latex', 'FontSize', 15);
    set(gca, 'XAxisLocation', 'top', 'YDir', 'reverse');
    
    % --- AXIS HANDLING ---
    if i < 3
        % Plots 1 & 2: Hide Y-axis numbers completely
        set(gca, 'YTickLabel', []);
    else
        % Plot 3: Move axis to RIGHT and show label
        set(gca, 'YAxisLocation', 'right');
        % Standard rotation (90) reads bottom-to-top.
        ylabel('$\tau_\ell$', 'Interpreter', 'latex', 'FontSize', 15);
    end
end

% Shared Elements
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Cosine Similarity';
cb.Label.FontSize = 12;