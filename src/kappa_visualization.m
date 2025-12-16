clear; clc; close all;

% --- Parameters ---
Sigma_w = eye(2);
Sigma   = diag([1, 40]); % Anisotropic covariance
d = diag(inv(Sigma));    % Damping coefficients (inverse variance)
sigma_w_scalar = 1;      

% Calculate Theoretical Parameters
alpha = (sigma_w_scalar .* d) / 2; 

% 1. Calculate Critical Frequency Threshold
% |w| = |alpha1 - alpha2| / (2 * sqrt(d1*d2))
w_crit = abs(alpha(1) - alpha(2)) / (2 * sqrt(d(1)*d(2)));

% Define the 3 Test Cases
omegas = [0.1 * w_crit,   ... % Case 1: Overdamped
          w_crit,         ... % Case 2: Critical
          2.5 * w_crit];      % Case 3: Oscillatory 

titles = {'Overdamped \Delta > 0', 'Critical \Delta \approx 0', 'Oscillatory \Delta < 0'};

% --- Visualization Setup ---
f = figure('Color','w');
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(t, 'Evolution of Dynamics: From Hyperbolic Loss to Elliptical Memory', 'FontSize', 10, 'FontWeight', 'bold');

% Common Lag settings
numLags = 500;
lastLag = 15; 
% lags = linspace(0, lastLag, numLags);
lags = linspace(-lastLag, lastLag, numLags);

for i = 1:3
    % 1. Setup System
    omega = omegas(i);
    S = omega * [0 1; -1 0];
    A = (-0.5 * Sigma_w + S) / Sigma; 
    
    % 2. Calc Theoretical Constants
    mu = trace(A)/2;
    Delta = trace(A)^2 - 4*det(A);
    
    % --- ROBUST CALCULATION BLOCK ---
    if abs(Delta) < 1e-5 % Critical Regime (Handle Limit explicitly)
        
        % In the limit Delta->0, the hyperbolic/trig functions become linear
        % u1 -> sqrt(2)
        % u2 -> Linear growth slope determined by rotation vs damping
        u1 = sqrt(2) * ones(size(lags));
        
        % Derived limit slope for u2 when Delta=0: omega*(d1+d2)
        limit_slope = omega * (d(1) + d(2)); 
        u2 = limit_slope * lags;
        
        kappa = Inf; % Technically infinite at singularity
        regime_color = [0.47, 0.67, 0.19]; % Green
        
    else % Non-Critical Cases (Standard Formulas)
        
        gamma = sqrt(abs(Delta))/2;
        J = (A - mu*eye(2)) / gamma;
        kappa = trace(J'*J);
        
        if Delta > 0 % Overdamped
            C = cosh(gamma * lags);
            S = sinh(gamma * lags);
            regime_color = [0.85, 0.33, 0.1]; % Burnt Orange
        else % Oscillatory
            C = cos(gamma * lags);
            S = sin(gamma * lags);
            regime_color = [0, 0.45, 0.74];   % Blue
        end
        
        u1 = sqrt(2) * C; 
        u2 = sqrt(kappa) * S;
    end
    
    % Compute Cosine Similarity Matrix
    U = [u1; u2]; 
    U_norm = U ./ vecnorm(U, 2, 1);
    SimMat = U_norm' * U_norm;
    SimMat(tril(true(size(SimMat)), -1)) = NaN;
    
    % --- PLOT TOP ROW: GEOMETRY (Lag Vector) ---
    nexttile(i);
    hold on; grid on; axis equal;
    
    % Visual guide lines for Overdamped
    if Delta > 1e-5 
        x_asym = linspace(0, max(u1), 100);
        slope = sqrt(kappa/2);
        plot(x_asym, slope*x_asym, '--k', 'Color', [0.7 0.7 0.7]); hold on
        plot(-x_asym, slope*x_asym, '--k', 'Color', [0.7 0.7 0.7]);
        plot(x_asym, -slope*x_asym, '--k', 'Color', [0.7 0.7 0.7]);
        plot(-x_asym, -slope*x_asym, '--k', 'Color', [0.7 0.7 0.7]);
        % max_u2 = max(u2); if max_u2 == 0, max_u2 = 1; end
        % ylim([0, max_u2*1.1]);
    end
    
    plot(u1, u2, 'LineWidth', 2, 'Color', regime_color); hold on
    plot(-u1, u2, '--', 'LineWidth', 2, 'Color', regime_color);
    
    xlabel('$u_1$', 'Interpreter', 'latex'); 
    if i==1, ylabel('$u_2$', 'Interpreter', 'latex'); end
    title(titles{i});
    
    % --- PLOT BOTTOM ROW: SIMILARITY MATRIX ---
    nexttile(i + 3);
    h = imagesc(lags, lags, SimMat);
    set(h, 'AlphaData', ~isnan(SimMat));
    set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
    colormap(gca, magma); 
    axis square;
    clim([-1 1]);
    
    if i==1
        ylabel('Lag \tau_2');
        xlabel('Lag \tau_1');
    end
    
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Cosine Similarity';

function cmap = custom_colormap()
    top = [178, 24, 43]/255;   % Red
    mid = [247, 247, 247]/255; % White
    bot = [33, 102, 172]/255;  % Blue
    len = 256;
    cmap = [linspace(bot(1),mid(1),len/2)', linspace(bot(2),mid(2),len/2)', linspace(bot(3),mid(3),len/2)'; 
            linspace(mid(1),top(1),len/2)', linspace(mid(2),top(2),len/2)', linspace(mid(3),top(3),len/2)'];
end