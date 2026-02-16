clear all; clc
iSub = 5;

% --- 1. SETUP ---
clusterCol = 'Macro_Area'; % Macro_Area or Functional_Network
dataDir = fullfile(pwd, 'data');
outDir = fullfile(dataDir, 'regressed_001_01_sim62131');

% Load Data
d = dir(fullfile(outDir, '*.mat'));
subj = load(fullfile(outDir, d(iSub).name)); 
T = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');

% Spectral Decomposition
results = get_kappa_spectrum(subj.A, eye(size(subj.A,1)) * subj.output.eff_conn.NoiseVar);
freqs = results.Oscillatory.frequencies;
kappas = results.Oscillatory.kappas;
W_osc = results.Oscillatory.node_weights;
W_real = results.Real.node_weights;

% Sort by Frequency
[sorted_freqs, idx] = sort(freqs);
sorted_logK = log(kappas(idx));
sorted_W_osc = W_osc(:, idx);

% Global Metrics (Scatter Sizing)
% Size = (Mode Energy / Total Brain Energy) * Scaling Factor
global_E = sum(sorted_W_osc, 1);
dots_size = (global_E / sum(W_osc(:))) * 2500; 
dots_size(dots_size < 10) = 10; % Minimum visibility

% Color = participation ratio of oscillatory mode across regions
dots_PR = (global_E.^2) ./ (size(subj.A,1) * sum(sorted_W_osc.^2, 1));

% --- 2. LOOP ---
u_clusters = unique(string(T.(clusterCol)));
n_clust = length(u_clusters);
O_C = zeros(n_clust, 1); PR_C = zeros(n_clust, 1);

% Initialize Heatmap
H_mat_Ell = zeros(n_clust, length(freqs)); % Weighted by log(kappa)

for i = 1:n_clust
    node_idx = string(T.(clusterCol)) == u_clusters(i);
    
    % 1. Oscillatory Fraction
    E_osc = sum(sum(sorted_W_osc(node_idx, :)));
    E_real = sum(sum(W_real(node_idx, :)));
    O_C(i) = E_osc / (E_osc + E_real);
    
    % 2. Modal Diversity (PR)
    E_modes = sum(sorted_W_osc(node_idx, :), 1);
    PR_C(i) = sum(E_modes)^2 / (length(E_modes) * sum(E_modes.^2));
    
    % Normalized energy distribution for this cluster
    norm_E = E_modes / sum(E_modes);

    % 3. Regional Ellipticity
    H_mat_Ell(i, :) = norm_E .* sorted_logK';
end

%% --- 3. VISUALIZATION ---
figure('Color','w', 'Position', [100 100 1100 600]); 
t = tiledlayout(3, 3, 'TileSpacing','compact', 'Padding','compact');

% === Tile 1: Dynamical Landscape (SCATTER) ===
% X-Axis: True Frequency (Continuous)
ax1 = nexttile([1 2]);
scatter(sorted_freqs, sorted_logK, dots_size, dots_PR, 'filled', 'MarkerEdgeColor','k');
ylabel('Ellipticity ($\log \kappa$)', 'Interpreter','latex'); 
title('\textbf{Dynamical Landscape}', 'Interpreter','latex'); 
grid on; box on;
c1 = colorbar; set(c1.Label, 'String', 'Spatial Reach ($\mathrm{PR}$)', 'Interpreter', 'latex');
colormap(ax1, magma);
clim(ax1, [0 1]);
xlim([0 max(sorted_freqs)*1.05]); % Slight padding to not cram the dot(s) at the end

% === Tile 2: Criticality Heatmap (High Shear) ===
% X-Axis: Mode Index (Discrete)
ax2 = nexttile([2 2]);
imagesc(H_mat_Ell); 
set(gca, 'YDir', 'normal');
yticks(1:n_clust); yticklabels(u_clusters);
title('\textbf{Oscillatory Criticality} ($\eta \cdot \log \kappa$)', 'Interpreter','latex');
c2 = colorbar;
colormap(ax2, magma);
xlabel('Frequency (Hz)', 'Interpreter','latex'); 

% Shared X-Axis Ticks (Frequency) applied to bottom plot
[~, unique_idx] = unique(round(sorted_freqs, 2), 'stable');
xticks(unique_idx); 
xticklabels(string(round(sorted_freqs(unique_idx), 2)));
xtickangle(45);

% === Tile 3: Composition (BAR CHART) ===
ax3 = nexttile([3 1]);
b = barh(O_C);
b.FaceColor = 'flat'; b.CData = PR_C;
yticks(1:n_clust); yticklabels(u_clusters); 
title('\textbf{Oscillatory Fraction} ($O_{\mathcal{C}}$)', 'Interpreter','latex'); 
grid on; box on;
colormap(ax3, magma);
c3 = colorbar; set(c3.Label, 'String', 'Modal Diversity ($\mathrm{PR}$)', 'Interpreter', 'latex');
clim(ax3, [0 1]);
