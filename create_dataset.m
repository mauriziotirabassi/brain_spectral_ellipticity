function T_final = create_dataset(dataDir, clusterCol)
% CREATE_DATASET Generates a feature matrix.

    % --- 1. Setup ---
    outDir = fullfile(dataDir, 'regressed_001_01_sim62131');
    T_table = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
    
    node_labels = string(T_table.(clusterCol));
    u_clusters = unique(node_labels);
    
    d = dir(fullfile(outDir, '*.mat'));
    files = {d.name};
    n_subj = length(files);
    
    rows_cell = cell(n_subj, 1); 
    
    fprintf('Processing %d subjects...\n', n_subj);

    % --- 2. Main Loop ---
    for i = 1:n_subj
        % Load Data
        subj = load(fullfile(outDir, files{i}));
        A = subj.A;
        Sigma_w = eye(size(A,1)) * subj.output.eff_conn.NoiseVar;
        
        % Spectral Decomposition
        res = get_kappa_spectrum(A, Sigma_w);
        W_osc = res.Oscillatory.node_weights; 
        log_K = log(res.Oscillatory.kappas); 
        
        % Pre-calculations
        W_osc_tot = sum(W_osc, 2);
        W_tot_recon = W_osc_tot + sum(res.Real.node_weights, 2);
        
        % Global Stats
        E_modes_glob = sum(W_osc, 1);      
        p_modes_glob = E_modes_glob / sum(E_modes_glob); 
        
        % Topological Matrices
        Rho = W_osc ./ sum(W_osc, 1); 
        Eta = W_osc ./ sum(W_osc, 2); 

        % --- 3. Build Subject Row ---
        row = struct();
        row.SubjectID = i;
        
        % Global Metrics
        sorted_K = sort(log_K, 'descend');
        row.Global_MaxK1 = sorted_K(1);
        row.Global_MaxK2 = sorted_K(2);
        row.Global_MeanK = mean(log_K);
        row.Global_MeanK_W = p_modes_glob * log_K; 
        % TODO: row.Global_Std = sqrt();
        row.Global_StdK_W = sqrt(p_modes_glob * (log_K - row.Global_MeanK_W).^2);

        % Cluster Metrics
        for c = 1:length(u_clusters)
            raw_name = u_clusters(c);
            clean_name = matlab.lang.makeValidName(raw_name); 
            idx = node_labels == raw_name; 
            
            % --- A. Energy & Diversity ---
            E_osc_C = sum(W_osc_tot(idx));
            E_tot_C = sum(W_tot_recon(idx));
            
            % Raw Oscillatory Energy (The Volume)
            row.(sprintf('%s_EnergyRaw', clean_name)) = E_osc_C;
            % TODO: Add weighted oscillatory energy where each oscillatory
            % mode enegry is weighte dby its \kappa
            
            % Standard Osc Fraction
            row.(sprintf('%s_OscFrac', clean_name)) = E_osc_C / E_tot_C;
            
            % Diversity (PR)
            E_modes_C = sum(W_osc(idx, :), 1); 
            row.(sprintf('%s_PR', clean_name)) = sum(E_modes_C)^2 / (length(E_modes_C) * sum(E_modes_C.^2));
            
            % --- B. Integrator Profile ---
            p_modes_C = E_modes_C / sum(E_modes_C);
            mu_C = p_modes_C * log_K;
            
            % Instability Load (The Pressure - Total Energy weighted by K)
            row.(sprintf('%s_EnergyW', clean_name)) = E_modes_C * log_K; 
            
            % Instability-Weighted Frequency (The Speed of Chaos)
            % (Sum(E * logK * f) / Sum(E * logK))
            weighted_energy_dist = E_modes_C' .* log_K;
            row.(sprintf('%s_FreqW', clean_name)) = sum(weighted_energy_dist .* freqs) / sum(weighted_energy_dist);
            
            row.(sprintf('%s_MeanK', clean_name)) = mu_C;
            row.(sprintf('%s_StdK', clean_name)) = sqrt(p_modes_C * (log_K - mu_C).^2);
            
            % --- C. Topology ---
            Rho_C_vec = sum(Rho(idx, :), 1); 
            row.(sprintf('%s_Hubness', clean_name)) = sum(Rho_C_vec);
            row.(sprintf('%s_Hubness_W', clean_name)) = Rho_C_vec * log_K;
            
            Eta_C_vec = sum(Eta(idx, :), 1);
            row.(sprintf('%s_Dispersion', clean_name)) = gini_coeff(Eta_C_vec);
            row.(sprintf('%s_Dispersion_W', clean_name)) = gini_coeff(Eta_C_vec' .* log_K);
        end
        
        rows_cell{i} = row;
    end
    
    % --- 4. Finalize ---
    T_final = struct2table([rows_cell{:}], 'AsArray', true);
    fprintf('Done. Dataset: %d subjects x %d features.\n', n_subj, width(T_final));
end

function G = gini_coeff(x)
    x = abs(x(:));
    if sum(x) == 0, G = 0; return; end
    n = length(x);
    x = sort(x);
    G = (2 * sum((1:n)' .* x) / (n * sum(x))) - (n + 1) / n;
end