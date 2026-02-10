clear all; clc; close all;

% --- Parameters ---
clusterCol = 'Network'; 
dominance_thresh = 0.05; % Threshold for Min/Max K calculation
dataDir = fullfile(pwd,'data');
outDir  = fullfile(dataDir,'regressed_001_01_sim62131');

% --- Load File List ---
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};
num_subs = length(files);
sub_labels = string(1:num_subs)'; % ID: 1, 2, 3...

% Load Region Names
T_table = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
node_labels = string(T_table.(clusterCol)); 
u_clusters = unique(node_labels); 
n_clust = length(u_clusters);

% --- Define Feature Set (6 per Cluster) ---
metric_names = ["OscFrac", "PR", "MeanLogK", "StdLogK", "MinLogK", "MaxLogK"];
num_metrics = length(metric_names);
num_features = n_clust * num_metrics; 

Subject_Features = zeros(num_subs, num_features);
Feature_Labels = strings(1, num_features);

fprintf('Processing %d subjects for 3D Fingerprinting...\n', num_subs);

%% 1. Feature Extraction Loop
for i = 1:num_subs
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
    
    W_osc_total = sum(W_osc_mat, 2); 
    W_total_reconstructed = W_osc_total + sum(W_real_mat, 2);
    Total_Mode_Energy_All = sum(W_osc_mat, 1);
    
    % Cluster Loop
    feat_idx = 1;
    for c = 1:n_clust
        idx = node_labels == u_clusters(c);
        
        % A. Oscillatory Fraction (O_C)
        cluster_osc = sum(W_osc_total(idx));
        cluster_tot = sum(W_total_reconstructed(idx));
        O_C = cluster_osc / cluster_tot;
        
        % B. Mode Statistics (Weighted by Mode Energy)
        mode_energies = sum(W_osc_mat(idx, :), 1);
        
        if sum(mode_energies) > 0
            p_C = mode_energies / sum(mode_energies);
            mu_K = sum(p_C .* log_kappas');
            sigma_K = sqrt(sum(p_C .* (log_kappas' - mu_K).^2));
            
            % PR Broadband
            num = sum(mode_energies)^2;
            den = length(mode_energies) * sum(mode_energies.^2);
            PR_C = num / den;
        else
            mu_K = 0; sigma_K = 0; PR_C = 0;
        end
        
        % C. Min/Max K (Based on Dominance)
        P_C_cluster = mode_energies ./ Total_Mode_Energy_All;
        dominant_indices = P_C_cluster > dominance_thresh;
        
        if any(dominant_indices)
            dom_k = log_kappas(dominant_indices);
            min_K = min(dom_k);
            max_K = max(dom_k);
        else
            min_K = NaN; max_K = NaN; % Handle later
        end
        
        % Store 6 Features
        Subject_Features(i, feat_idx)   = O_C;
        Subject_Features(i, feat_idx+1) = PR_C;
        Subject_Features(i, feat_idx+2) = mu_K;
        Subject_Features(i, feat_idx+3) = sigma_K;
        Subject_Features(i, feat_idx+4) = min_K;
        Subject_Features(i, feat_idx+5) = max_K;
        
        % Create Labels (Once)
        if i == 1
            for m = 1:num_metrics
                Feature_Labels(feat_idx+m-1) = u_clusters(c) + "_" + metric_names(m);
            end
        end
        feat_idx = feat_idx + num_metrics;
    end
end

% --- Data Cleaning (NaN Imputation) ---
% Replace NaNs in Min/Max K (occuring if no dominant mode) with 0
Subject_Features(isnan(Subject_Features)) = 0;

%% 2. Display Feature Table Head
disp('--- Feature Table Preview (First 5 Subjects, First 6 Cols) ---');
T_display = array2table(Subject_Features(1:5, 1:6), ...
    'VariableNames', Feature_Labels(1:6), ...
    'RowNames', sub_labels(1:5));
disp(T_display);

%% 3. 3D PCA & Clustering Analysis
% A. Standardization (Z-Score)
X = zscore(Subject_Features);

% B. PCA (Extract Top 3 Components)
[coeff, score, latent, tsquared, explained] = pca(X);

% C. Hierarchical Clustering (Ward's Linkage on 3D PCA scores)
Z_link = linkage(score(:, 1:3), 'ward'); 
clusters = cluster(Z_link, 'MaxClust', 3); % Cut tree into 3 groups for coloring

%% 4. Visualization
figure('Color','w', 'Position', [50 50 1600 900]);
t = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');

% --- Tile 1: Scree Plot ---
nexttile;
pareto(explained);
xlabel('Principal Component'); ylabel('Variance Explained (%)');
title('Scree Plot'); grid on;

% --- Tile 2: Dendrogram ---
nexttile;
[H, T, outperm] = dendrogram(Z_link, 0, 'Labels', sub_labels, 'Orientation', 'left');
title('Subject Clustering (Ward Linkage)'); xlabel('Distance');

% --- Tile 3: 3D Biplot (Spanning Bottom Row) ---
nexttile([1, 2]); 

% Find top features contributing to PC1, PC2, PC3 combined
% Magnitude of feature vector in 3D space
feature_magnitudes = sum(coeff(:, 1:3).^2, 2);
[~, sort_idx] = sort(feature_magnitudes, 'descend');
top_n = 12; % Show top 12 most important arrows
top_features = sort_idx(1:top_n);

% 3D Biplot
h = biplot(coeff(top_features, 1:3), 'Scores', score(:, 1:3), ...
    'VarLabels', Feature_Labels(top_features), ...
    'ObsLabels', sub_labels, ...
    'MarkerSize', 20);

% Customizing Appearance
title(sprintf('3D Dynamical State Space (PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%)', ...
    explained(1), explained(2), explained(3)));
xlabel(['PC 1 (' num2str(explained(1), '%.1f') '%)']);
ylabel(['PC 2 (' num2str(explained(2), '%.1f') '%)']);
zlabel(['PC 3 (' num2str(explained(3), '%.1f') '%)']);
grid on; view(45, 30); % Set isometric view