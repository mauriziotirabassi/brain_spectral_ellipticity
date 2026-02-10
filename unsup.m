clear all; clc; close all;

% --- Parameters ---
clusterCol = 'Network'; 
dominance_thresh = 0.05; 
dataDir = fullfile(pwd,'data');
outDir  = fullfile(dataDir,'regressed_001_01_sim62131');

% --- Load File List ---
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};
num_subs = length(files);

% Load Region Names (Assuming consistent across subjects)
T_table = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
node_labels = string(T_table.(clusterCol)); 
u_clusters = unique(node_labels); 
n_clust = length(u_clusters);

% --- Initialize Feature Matrix ---
% 3 metrics per cluster: Osc_Frac, Mean_LogK, PR_Broadband
num_features = n_clust * 3; 
Subject_Features = zeros(num_subs, num_features);
Feature_Labels = strings(1, num_features);

fprintf('Processing %d subjects...\n', num_subs);

%% 1. Feature Extraction Loop
for i = 1:num_subs
    fprintf('  Subject %d: %s\n', i, files{i});
    
    % Load Subject
    subj = load(fullfile(outDir,files{i}));
    A = subj.A; 
    n = size(A, 1); 
    Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar; 
    
    % Run Spectral Decomposition
    results = get_kappa_spectrum(A, Sigma_w);
    
    % Extract Weights & Params
    W_osc_mat = results.Oscillatory.node_weights;   
    kappas    = results.Oscillatory.kappas;
    log_kappas = log(kappas);
    W_real_mat = results.Real.node_weights;         
    
    % Node-level stats
    W_osc_total = sum(W_osc_mat, 2); 
    W_total_reconstructed = W_osc_total + sum(W_real_mat, 2);
    
    % Cluster Loop
    feat_idx = 1;
    for c = 1:n_clust
        idx = node_labels == u_clusters(c);
        
        % A. Oscillatory Fraction (O_C)
        cluster_osc = sum(W_osc_total(idx));
        cluster_tot = sum(W_total_reconstructed(idx));
        O_C = cluster_osc / cluster_tot;
        
        % B. Mean Log Kappa (Mu_K) - Weighted by Mode Energy in Cluster
        mode_energies = sum(W_osc_mat(idx, :), 1);
        if sum(mode_energies) > 0
            p_C = mode_energies / sum(mode_energies);
            mu_K = sum(p_C .* log_kappas');
        else
            mu_K = 0; % Fallback
        end
        
        % C. Participation Ratio (PR)
        num = sum(mode_energies)^2;
        den = length(mode_energies) * sum(mode_energies.^2);
        if den > 0, PR_C = num / den; else, PR_C = 0; end
        
        % Store Features
        Subject_Features(i, feat_idx)   = O_C;
        Subject_Features(i, feat_idx+1) = mu_K;
        Subject_Features(i, feat_idx+2) = PR_C;
        
        % Create Labels (only once)
        if i == 1
            Feature_Labels(feat_idx)   = u_clusters(c) + "_OscFrac";
            Feature_Labels(feat_idx+1) = u_clusters(c) + "_LogK";
            Feature_Labels(feat_idx+2) = u_clusters(c) + "_PR";
        end
        feat_idx = feat_idx + 3;
    end
end

%% 2. PCA & Clustering Analysis

% A. Standardization (Z-Score)
% Essential because LogK and OscFrac have different scales
X = zscore(Subject_Features);

% B. PCA
[coeff, score, latent, tsquared, explained] = pca(X);

% C. Hierarchical Clustering (on PCA scores or raw standardized data)
% Using Ward's linkage for compact clusters
Z_link = linkage(score(:, 1:5), 'ward'); % Use top 5 PCs for clustering

%% 3. Visualization

figure('Color','w', 'Position', [100 100 1600 900]);
t = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');

% --- Tile 1: Scree Plot (Variance Explained) ---
nexttile;
pareto(explained);
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot: Dimensionality of Phenotypes');
grid on;

% --- Tile 2: Dendrogram (The Family Tree) ---
nexttile;
[H, T, outperm] = dendrogram(Z_link, 0, 'Labels', files, 'Orientation', 'left');
title('Subject Clustering (Dynamical Similarity)');
xlabel('Distance (Ward)');
set(gca, 'TickLabelInterpreter', 'none'); 

% --- Tile 3: PCA Biplot (PC1 vs PC2) ---
nexttile([1, 2]); % Span 2 columns
% Find top features contributing to PC1 and PC2 for cleaner plot
[~, sort_idx] = sort(sum(coeff(:,1:2).^2, 2), 'descend');
top_n = 15; % Show only top 15 most important features arrows
top_features = sort_idx(1:top_n);

biplot(coeff(top_features, 1:2), 'Scores', score(:, 1:2), ...
    'VarLabels', Feature_Labels(top_features), ...
    'MarkerSize', 20);
title(sprintf('PCA Biplot (PC1: %.1f%%, PC2: %.1f%%)', explained(1), explained(2)));
xlabel(['PC 1 (' num2str(explained(1), '%.1f') '%)']);
ylabel(['PC 2 (' num2str(explained(2), '%.1f') '%)']);
grid on;

%%
function props = calculate_kappa_properties_fixed(B)
    % Extract 2x2 invariants
    off_diff = B(2,1) - B(1,2);
    h        = (B(1,1) - B(2,2)) / 2;
    gamma    = (B(2,1) + B(1,2)) / 2;
    
    omega = abs(off_diff) / 2;
    delta_sq = h^2 + gamma^2;
    omega_sq = omega^2;
    
    num = omega_sq + delta_sq;
    den = abs(omega_sq - delta_sq);
    
    if den < 1e-9 * num, kappa = Inf; else, kappa = num / den; end
    
    props.kappa     = kappa;
    props.frequency = sqrt(den); % CORRECTED: Square root of discriminant
end