function T_final = create_dataset_v2(dataDir, clusterCol)
% CREATE_DATASET_V2 Generates a dataset with regional entropy, weighted sum, 
% and weighted variance for ellipticity distributions.
    
    % --- 1. Setup ---
    outDir = fullfile(dataDir, 'regressed_001_01_sim62131');
    T_table = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
    
    node_labels = string(T_table.(clusterCol));
    u_clusters = unique(node_labels);
    n_clust = length(u_clusters);
    
    d = dir(fullfile(outDir, '*.mat'));
    files = {d.name};

    n_subj = length(files);
    
    rows_cell = cell(n_subj, 1); 
    fprintf('Processing %d subjects...\n', n_subj);
    
    % --- 2. Main Loop over subjects ---
    for i = 1:n_subj
        % Load subject data
        subj = load(fullfile(outDir, files{i}));
        A = subj.A;
        Sigma_w = eye(size(A,1)) * subj.output.eff_conn.NoiseVar;
        
        % Spectral decomposition
        res = get_kappa_spectrum(A, Sigma_w);
        W_osc = res.Oscillatory.node_weights;  % nodes x modes
        kappas = res.Oscillatory.kappas;       % mode-specific kappa
        
        % Sort indices by frequency
        [~, freq_idx] = sort(res.Oscillatory.frequencies); 
        W_osc = W_osc(:, freq_idx);
        kappas = kappas(freq_idx);
        
        % Pre-calculate log(kappa) vector
        log_k = log(kappas'); 
        
        % Calculate max entropy (log K) for normalization
        K = length(kappas);
        max_H = log(K);
        
        % Initialize row struct for this subject
        row = struct();
        row.SubjectID = i;
        
        % --- Loop over clusters ---
        for c = 1:n_clust
            raw_name = u_clusters(c);
            clean_name = matlab.lang.makeValidName(raw_name);
            node_idx = node_labels == raw_name;
            
            % Extract energy distribution (eta)
            E_modes = sum(W_osc(node_idx, :), 1); 
            eta = E_modes / (sum(E_modes) + eps); 
            
            % --- Ellipticity Metrics ---
            
            % 1. Weighted Sum (Mean Ellipticity)
            % Contribution of each mode to the region's ellipticity
            ellip_contrib = eta .* log_k;   
            ellip_w_sum = sum(ellip_contrib);
            
            % 2. Weighted Variance (Dynamical Diversity)
            % Variance of log(kappa) weighted by energy (eta)
            ellip_var = sum(eta .* (log_k - ellip_w_sum).^2);
            
            % Store results
            row.(sprintf('%s_SumW', clean_name)) = ellip_w_sum;
            row.(sprintf('%s_VarW', clean_name)) = ellip_var;
        end
        
        rows_cell{i} = row;
    end
    
    % --- 3. Convert to table ---
    T_final = struct2table([rows_cell{:}], 'AsArray', true);
    fprintf('Done. Dataset: %d subjects x %d features.\n', n_subj, width(T_final));
end