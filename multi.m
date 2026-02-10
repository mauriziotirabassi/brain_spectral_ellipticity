clear all; clc

% --- Parameters ---
clusterCol = 'Macro'; 
var_thresh = 0.5;
dataDir = fullfile(pwd,'data');
outDir  = fullfile(dataDir,'regressed_001_01_sim62131');

% --- Load File List ---
d = dir(fullfile(outDir,'*.mat')); 
files = {d.name};
n_subs = length(files);

% --- Load Cluster Names (Once) ---
T = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
node_labels = string(T.(clusterCol)); 
u_clusters = unique(node_labels);
n_clust = length(u_clusters);

fprintf('Found %d subjects. Starting analysis...\n', n_subs);

% --- Main Subject Loop ---
for iSub = 1:n_subs
    fprintf('\n========================================\n');
    fprintf('Processing Subject %d: %s\n', iSub, files{iSub});
    fprintf('========================================\n');
    
    try
        % --- Load Data ---
        subj = load(fullfile(outDir,files{iSub}));
        A = subj.A; 
        n = size(A, 1); 
        Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar; 
        
        % --- 1. Calculate Steady-State Covariance (Lyapunov) ---
        Sigma = lyap(A, Sigma_w);
        node_variances = diag(Sigma); 

        % --- 2. Run Spectral Decomposition ---
        % Using evalc to suppress output from get_kappa_spectrum if it is chatty
        [~, results] = evalc('get_kappa_spectrum(A, Sigma_w)'); 
        
        W = results.node_weights;      % Node Energy (Nodes x Modes)
        kappas = results.kappas;       % Intrinsic non-normality
        log_kappas = log(kappas);      % Log scale
        
        % --- 3. Cluster Statistics Calculation ---
        mu_C    = zeros(n_clust, 1);
        sigma_C = zeros(n_clust, 1);
        min_K   = zeros(n_clust, 1); 
        max_K   = zeros(n_clust, 1); 
        O_C     = zeros(n_clust, 1); 

        for i = 1:n_clust
            idx = node_labels == u_clusters(i);
            
            % --- Cluster Oscillatory Fraction (O_C) ---
            % Numerator: Sum of all oscillatory energy in the cluster
            cluster_osc_energy = sum(sum(W(idx, :)));
            % Denominator: Sum of total variance (Sigma) in the cluster
            cluster_total_var = sum(node_variances(idx));
            
            O_C(i) = cluster_osc_energy / cluster_total_var;
            
            % --- Spectral Statistics (Weighted by Oscillatory Energy) ---
            mode_energies = sum(W(idx, :), 1);
            p_C = mode_energies / sum(mode_energies); 
            
            % Weighted Statistics
            mu_C(i)    = sum(p_C .* log_kappas');
            sigma_C(i) = sqrt(sum(p_C .* (log_kappas' - mu_C(i)).^2));
            
            % Complex Metrics via Functions
            [min_K(i), max_K(i)] = calculate_robust_range(log_kappas, p_C, var_thresh);
        end
        
        % --- Output Table ---
        res_table = table(u_clusters, O_C, mu_C, sigma_C, min_K, max_K, ...
            'VariableNames', {'Cluster', 'Osc_Frac', 'Mean_LogK', 'Std_LogK', 'Min_LogK', 'Max_LogK'});
        disp(res_table);
        
    catch ME
        fprintf('Error on Subject %d: %s\n', iSub, ME.message);
    end
end

% --- Helper Functions ---
function [min_k, max_k] = calculate_robust_range(log_kappas, p_C, thresh)
    [sorted_p, sort_idx] = sort(p_C, 'descend');
    cutoff_idx = find(cumsum(sorted_p) >= thresh, 1, 'first');
    if isempty(cutoff_idx), cutoff_idx = length(p_C); end
    active_kappas = log_kappas(sort_idx(1:cutoff_idx));
    min_k = min(active_kappas);
    max_k = max(active_kappas);
end