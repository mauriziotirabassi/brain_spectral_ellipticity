clear all; clc

% --- Parameters ---
clusterCol = 'Network'; 
dataDir = fullfile(pwd,'data');
outDir  = fullfile(dataDir,'regressed_001_01_sim62131');

% --- Load Data ---
d = dir(fullfile(outDir,'*.mat')); 
files = {d.name};
iSub = 1; 
subj = load(fullfile(outDir,files{iSub}));

A = subj.A; 
n = size(A, 1); 
Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar; 
T_table = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');

% --- Run Spectral Decomposition ---
% Assumes get_kappa_spectrum is in the path
results = get_kappa_spectrum(A, Sigma_w);

% 1. Extract Oscillatory Data
W_osc_mat = results.Oscillatory.node_weights;   
kappas    = results.Oscillatory.kappas;         
freqs     = results.Oscillatory.frequencies;    
log_kappas = log(kappas);

% 2. Extract Real Data
W_real_mat = results.Real.node_weights;         

% --- 3. Calculate Nodal Oscillatory Fraction (O_i) ---
W_osc_total = sum(W_osc_mat, 2); 
W_real_total = sum(W_real_mat, 2);
W_total_reconstructed = W_osc_total + W_real_total;
% O_i = W_osc_total ./ W_total_reconstructed;

% --- Participation Ratio (Mode Spatial Extent) ---
PR_modes = (sum(W_osc_mat, 1).^2) ./ (n * sum(W_osc_mat.^2, 1)); 

% --- Cluster Setup ---
node_labels = string(T_table.(clusterCol)); 
u_clusters = unique(node_labels); 
n_clust = length(u_clusters);

% --- Cluster Statistics Loop ---
mu_C    = zeros(n_clust, 1);
sigma_C = zeros(n_clust, 1);
O_C     = zeros(n_clust, 1); 
PR_C    = zeros(n_clust, 1); 
E_C_modes = zeros(n_clust, length(kappas)); 

for i = 1:n_clust
    idx = node_labels == u_clusters(i);
    
    % Region oscillatory fraction
    cluster_osc_energy = sum(W_osc_total(idx));
    cluster_total_energy = sum(W_total_reconstructed(idx));
    O_C(i) = cluster_osc_energy / cluster_total_energy;
    
    mode_energies = sum(W_osc_mat(idx, :), 1);
    E_C_modes(i, :) = mode_energies; 
    
    % Region modal diversity (oscillatory fraction)
    num = sum(mode_energies)^2;
    den = length(mode_energies) * sum(mode_energies.^2);
    PR_C(i) = num / den;
    
    eta_C = mode_energies / sum(mode_energies); % Region spectral profile
    mu_C(i)    = sum(eta_C .* log_kappas'); % Region (weighted) mean non-normality
    sigma_C(i) = sqrt(sum(eta_C .* (log_kappas' - mu_C(i)).^2)); % Region (weighted) non-normality std
end

% --- Text Output ---
res_table = table(u_clusters, O_C, PR_C, mu_C, sigma_C, ...
    'VariableNames', {'Cluster', 'Osc_Frac', 'Modal_Diversity', 'Mean_LogK', 'Std_LogK'});
disp('Cluster Statistics:');
disp(res_table);

% --- Visualization Preparation ---
[sorted_freqs, sort_idx] = sort(freqs, 'ascend');
sorted_logK = log_kappas(sort_idx);
sorted_PR_modes = PR_modes(sort_idx);

Rho_C_mat = E_C_modes ./ sum(E_C_modes, 1); 
sorted_Rho_C = Rho_C_mat(:, sort_idx);

% Eta_C_mat = E_C_modes ./ sum(E_C_modes, 2); 
E_C_weighted = E_C_modes .* reshape(log_kappas, 1, []);
Eta_C_mat = E_C_weighted ./ sum(E_C_weighted, 2); % Weighted by \log\kappa
sorted_Eta_C = Eta_C_mat(:, sort_idx);

%% --- Visualization ---
figure('Color','w', 'Position', [100 100 1100 600]);
t = tiledlayout(3, 3, 'TileSpacing','compact', 'Padding', 'compact');
colormap(magma); 

% Prepare shared X-axis labels
n_modes = length(sorted_logK);
num_ticks = 25; 
tick_indices = round(linspace(1, n_modes, num_ticks));
tick_labels = round(sorted_logK(tick_indices), 1); 
xLabelHeat = 'Low Frequency $\leftarrow \rightarrow$ High Frequency (Values: $\log \kappa$)';

% --- Updated TILE 1: Dynamical Landscape (Energy Scaled) ---
E_global_modes = sum(W_osc_mat, 1); 
marker_sizes = 50 + (500 * (E_global_modes(sort_idx) / max(E_global_modes)));
nexttile([1 2]);
scatter(sorted_freqs, sorted_logK, marker_sizes, sorted_PR_modes, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
xlabel('Frequency (Hz)', 'Interpreter', 'latex'); 
ylabel('Integration Capacity ($\log \kappa$)', 'Interpreter', 'latex');
title('\textbf{Dynamical Landscape}: Size = Mode Energy, Color = Spatial Reach ($PR$)', 'Interpreter', 'latex');
grid on; box on; 
c1 = colorbar;
c1.Label.String = 'Spatial Reach ($PR$)';
c1.Label.Interpreter = 'latex';
clim([0 1]);

% --- TILE 2: Spectral Profile (Heatmap Eta) ---
nexttile([1 2]);
imagesc(sorted_Eta_C);
set(gca, 'YDir', 'normal'); 
c4 = colorbar;
c4.Label.String = '$\eta$';
c4.Label.Interpreter = 'latex';
clim([0 max(sorted_Eta_C(:))]); 
title('Region Spectral Profile: How much does a mode drive a region?', 'Interpreter', 'latex');
xticks(tick_indices);
xticklabels(string(tick_labels));
xlabel(xLabelHeat, 'Interpreter', 'latex');
yticks(1:n_clust);
yticklabels(u_clusters);

% --- TILE 3: Mode Topography (Heatmap Rho) ---
nexttile([1 2]);
imagesc(sorted_Rho_C);
set(gca, 'YDir', 'normal'); 
c2 = colorbar;
c2.Label.String = '$\rho$';
c2.Label.Interpreter = 'latex';
clim([0 max(sorted_Rho_C(:))]); 
title('Mode Spatial Topography: How much does a region contribute to a mode?', 'Interpreter', 'latex');
xticks(tick_indices);
xticklabels(string(tick_labels));
xlabel(xLabelHeat, 'Interpreter', 'latex');
yticks(1:n_clust);
yticklabels(u_clusters);

% --- TILE 4: Oscillatory Composition (Bar Chart) ---
nexttile([3 1]);
b = barh(O_C);
b.FaceColor = 'flat';      
b.CData = PR_C;            
yticks(1:n_clust);
yticklabels(u_clusters);
xlabel('Oscillatory Fraction ($O_{\mathcal{C}}$)', 'Interpreter', 'latex');
title('Region Oscillatory Fraction: How much is a region oscillatory?', 'Interpreter', 'latex');
xlim([0 1]); 
grid on;
c3 = colorbar;
c3.Label.String = 'Region Modal Diversity (Participation Ratio)';
c3.Label.Interpreter = 'latex';
clim([0 1]);
