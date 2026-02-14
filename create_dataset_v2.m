function T_final = create_dataset_v2(dataDir, clusterCol)
% CREATE_DATASET_V2 Generates a dataset with regional entropy for criticality 
% and coherence distributions. This is the conservative version focusing 
% exclusively on dynamical complexity.

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

        % Initialize row struct for this subject
        row = struct();
        row.SubjectID = i;

        % --- Loop over clusters ---
        for c = 1:n_clust
            raw_name = u_clusters(c);
            clean_name = matlab.lang.makeValidName(raw_name);
            node_idx = node_labels == raw_name;

            % Extract energy distribution over modes for this cluster
            E_modes = sum(W_osc(node_idx, :), 1); 
            eta = E_modes / (sum(E_modes) + eps); 

            % 1. Criticality Distribution (instability-weighted)
            crit_values = eta .* log(kappas');   
            crit_norm = crit_values / (sum(crit_values) + eps); 
            row.(sprintf('%s_Entropy_Crit', clean_name)) = -sum(crit_norm .* log(crit_norm + eps));

            % 2. Coherence Distribution (stability-weighted)
            coh_values = eta ./ (kappas' + eps);        
            coh_norm = coh_values / (sum(coh_values) + eps); 
            row.(sprintf('%s_Entropy_Coh', clean_name)) = -sum(coh_norm .* log(coh_norm + eps));
        end
        
        rows_cell{i} = row;
    end

    % --- 3. Convert to table ---
    T_final = struct2table([rows_cell{:}], 'AsArray', true);
    fprintf('Done. Dataset: %d subjects x %d features.\n', n_subj, width(T_final));
end