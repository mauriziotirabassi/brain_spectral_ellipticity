clear all; clc

% --- Parameters ---
clusterCol = 'Network'; 
dominance_thresh = 0.05; % Cluster must own >5% of mode energy to count
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
O_i = W_osc_total ./ W_total_reconstructed;

% --- Participation Ratio (PR) Calculation (Mode-Centric) ---
PR_modes = (sum(W_osc_mat, 1).^2) ./ (n * sum(W_osc_mat.^2, 1)); 

% --- Cluster Setup ---
node_labels = string(T_table.(clusterCol)); 
u_clusters = unique(node_labels); 
n_clust = length(u_clusters);

% --- Pre-Calculate Total Mode Energy (Global Denominator) ---
Total_Mode_Energy_All = sum(W_osc_mat, 1); 

% --- Cluster Statistics Loop ---
mu_C    = zeros(n_clust, 1);
sigma_C = zeros(n_clust, 1);
min_K   = zeros(n_clust, 1); 
max_K   = zeros(n_clust, 1); 
O_C     = zeros(n_clust, 1); 
PR_C    = zeros(n_clust, 1); 
E_C_modes = zeros(n_clust, length(kappas)); 

for i = 1:n_clust
    idx = node_labels == u_clusters(i);
    
    % Oscillatory Fraction
    cluster_osc_energy = sum(W_osc_total(idx));
    cluster_total_energy = sum(W_total_reconstructed(idx));
    O_C(i) = cluster_osc_energy / cluster_total_energy;
    
    % Cluster Energy per Mode
    mode_energies = sum(W_osc_mat(idx, :), 1);
    E_C_modes(i, :) = mode_energies; 
    
    % Cluster PR (Broadband Index)
    num = sum(mode_energies)^2;
    den = length(mode_energies) * sum(mode_energies.^2);
    PR_C(i) = num / den;
    
    % Weighted Statistics (Mean/Std based on Internal Energy)
    p_C = mode_energies / sum(mode_energies); 
    mu_C(i)    = sum(p_C .* log_kappas');
    sigma_C(i) = sqrt(sum(p_C .* (log_kappas' - mu_C(i)).^2));
    
    % --- NEW: Generator Range (Based on Spatial Dominance) ---
    % Calculate Participation Coefficient P_C for this cluster
    P_C_cluster = mode_energies ./ Total_Mode_Energy_All;
    
    % Find modes where this cluster is a dominant driver (>10%)
    dominant_indices = P_C_cluster > dominance_thresh;
    
    if any(dominant_indices)
        dominant_kappas = log_kappas(dominant_indices);
        min_K(i) = min(dominant_kappas);
        max_K(i) = max(dominant_kappas);
    else
        min_K(i) = NaN; 
        max_K(i) = NaN;
    end
end

% Create Results Table
res_table = table(u_clusters, O_C, PR_C, mu_C, sigma_C, min_K, max_K, ...
    'VariableNames', {'Cluster', 'Osc_Frac', 'PR_Broadband', 'Mean_LogK', 'Std_LogK', 'Min_LogK', 'Max_LogK'});
disp('Cluster Statistics:');
disp(res_table);

% --- Visualization Preparation ---
Total_Mode_Energy = sum(E_C_modes, 1);
P_C_mat = E_C_modes ./ Total_Mode_Energy; 

[sorted_logK, sort_k_idx] = sort(log_kappas, 'ascend');
sorted_freqs = freqs(sort_k_idx);
sorted_P_C    = P_C_mat(:, sort_k_idx); 
sorted_PR_modes = PR_modes(sort_k_idx);

%% --- Visualization ---
figure('Color','w', 'Position', [0 100 1800 350]);
t = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding', 'compact');
colormap(magma); 

% Tile 1: Scatter (Log Kappa vs Omega)
nexttile;
scatter(sorted_freqs, sorted_logK, 70, sorted_PR_modes, 'filled');
xlabel('Mode Frequency $\nu$', 'Interpreter', 'latex');
ylabel('Log Non-Normality ($\log \kappa$)', 'Interpreter', 'latex');
title('Intrinsic Non-Normality Spectrum');
grid on; box on; 
c1 = colorbar;
c1.Label.String = 'Mode Spatial PR';
clim([0 1]); 

% Tile 2: Heatmap (Cluster Participation)
nexttile;
imagesc(sorted_P_C);
set(gca, 'YDir', 'normal'); 
c2 = colorbar;
c2.Label.String = 'Participation P_C';
clim([0 0.5]); 
title('Cluster Spectral Profile');

% -- Update X-Axis to show Log Kappa values --
n_modes = length(sorted_logK);
num_ticks = 10; 
tick_indices = round(linspace(1, n_modes, num_ticks));
tick_labels = round(sorted_logK(tick_indices), 1); 

xticks(tick_indices);
xticklabels(string(tick_labels));
xlabel('Log Non-Normality ($\log \kappa$)', 'Interpreter', 'latex');
yticks(1:n_clust);
yticklabels(u_clusters);

% Tile 3: Bar Chart (Oscillatory Fraction)
nexttile;
b = barh(O_C);
b.FaceColor = 'flat';      
b.CData = PR_C;            
yticks(1:n_clust);
yticklabels(u_clusters);
xlabel('Oscillatory Fraction ($O_\mathcal{C}$)', 'Interpreter', 'latex');
title('Cluster Oscillatory Fraction');
xlim([0 1]); 
grid on;

c3 = colorbar;
c3.Label.String = 'Cluster Broadband PR';
clim([0 1]);
