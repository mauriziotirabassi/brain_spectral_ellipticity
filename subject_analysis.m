clear all; clc

iSub = 5;
genFontSize = 12; % Global font size for axes, labels, and titles

% Data setup
% "Macro_Area", "Functional_Network", "Division" (cortical/subcortical), "Abbreviation" (node)
clusterCol = 'Macro_Area';
dataDir = fullfile(pwd, 'data');
outDir = fullfile(dataDir, 'regressed_001_01_sim62131');

% Load data
d = dir(fullfile(outDir, '*.mat'));
subj = load(fullfile(outDir, d(iSub).name)); 
T = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');

% Spectral decomposition
results = get_kappa_spectrum(subj.A, eye(size(subj.A,1)) * subj.output.eff_conn.NoiseVar);
freqs = results.Oscillatory.frequencies;
kappas = results.Oscillatory.kappas;
W_osc = results.Oscillatory.node_weights;
W_real = results.Real.node_weights;

% Sort by frequency
[sorted_freqs, idx] = sort(freqs);
sorted_logK = log(kappas(idx));
sorted_W_osc = W_osc(:, idx);

% Calculate global energy and participation metrics
global_E = sum(sorted_W_osc, 1);
dots_size = (global_E / sum(W_osc(:))) * 2500; 
dots_size(dots_size < 10) = 10; 
dots_PR = (global_E.^2) ./ (size(subj.A,1) * sum(sorted_W_osc.^2, 1));

% Cluster loop
u_clusters = unique(string(T.(clusterCol)));
n_clust = length(u_clusters);
O_C = zeros(n_clust, 1); 
PR_C = zeros(n_clust, 1);
H_mat_Ell = zeros(n_clust, length(freqs)); 

for i = 1:n_clust
    node_idx = string(T.(clusterCol)) == u_clusters(i);
    
    % Oscillatory fraction
    E_osc = sum(sum(sorted_W_osc(node_idx, :)));
    E_real = sum(sum(W_real(node_idx, :)));
    O_C(i) = E_osc / (E_osc + E_real);
    
    % Modal diversity
    E_modes = sum(sorted_W_osc(node_idx, :), 1);
    PR_C(i) = sum(E_modes)^2 / (length(E_modes) * sum(E_modes.^2));
    
    % Regional ellipticity
    norm_E = E_modes / sum(E_modes);
    H_mat_Ell(i, :) = norm_E .* sorted_logK';
end

% Visualization
figure('Color', 'w', 'Position', [100 100 1100 600]); 
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Tile 1: Dynamical landscape
ax1 = nexttile([1 2]);
scatter(sorted_freqs, sorted_logK, dots_size, dots_PR, 'filled', 'MarkerEdgeColor', 'k');
ylabel('$\log \kappa$', 'Interpreter', 'latex', 'FontSize', genFontSize); 
title('\textbf{Dynamical Landscape}', 'Interpreter', 'latex', 'FontSize', genFontSize); 
grid on; box on;
c1 = colorbar; 
set(c1.Label, 'String', 'Spatial Reach', 'Interpreter', 'latex', 'FontSize', genFontSize);
colormap(ax1, magma);
clim(ax1, [0 1]);
xlim([0 max(sorted_freqs)*1.05]); 

% Tile 2: Ellipticity contribution
ax2 = nexttile([2 2]);
imagesc(H_mat_Ell); 
set(gca, 'YDir', 'normal');
yticks(1:n_clust); yticklabels(u_clusters);
title('\textbf{Ellipticity Contribution}', 'Interpreter', 'latex', 'FontSize', genFontSize);
c2 = colorbar;
set(c2.Label, 'String', '$\eta \log \kappa$', 'Interpreter', 'latex', 'FontSize', genFontSize);
colormap(ax2, magma);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', genFontSize); 

% Custom x-ticks for frequency
[~, unique_idx] = unique(round(sorted_freqs, 2), 'stable');
xticks(unique_idx); 
xticklabels(string(round(sorted_freqs(unique_idx), 2)));
xtickangle(45);

% Tile 3: Regional composition
ax3 = nexttile([3 1]);
b = barh(O_C);
b.FaceColor = 'flat'; b.CData = PR_C;
yticks(1:n_clust); yticklabels(u_clusters); 
title('\textbf{Regional Composition}', 'Interpreter', 'latex', 'FontSize', genFontSize); 
xlabel('Oscillatory Fraction ($O_{\mathcal{C}}$)', 'Interpreter', 'latex', 'FontSize', genFontSize); 
grid on; box on;
colormap(ax3, magma);
c3 = colorbar; 
set(c3.Label, 'String', 'Spectral Participation', 'Interpreter', 'latex', 'FontSize', genFontSize);
clim(ax3, [0 1]);

% Apply global font size to all axes
set([ax1, ax2, ax3], 'FontSize', genFontSize);