function T_final = create_dataset(dataDir, clusterCol)
% CREATE_DATASET Generates a dataset with weighted sum and weighted variance 
% for ellipticity distributions.

    % Setup
    outDir = fullfile(dataDir, 'regressed_001_01_sim62131');
    T_names = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
    
    node_labels = string(T_names.(clusterCol));
    u_clusters = unique(node_labels);
    n_clust = length(u_clusters);
    
    fileList = dir(fullfile(outDir, '*.mat'));
    files = {fileList.name};
    n_subj = length(files);
    
    rows_cell = cell(n_subj, 1); 
    fprintf('Processing %d subjects...\n', n_subj);
    
    % Process subjects
    for i = 1:n_subj
        % Load subject data
        subj = load(fullfile(outDir, files{i}));
        A = subj.A;
        Sigma_w = eye(size(A,1)) * subj.output.eff_conn.NoiseVar;
        
        % Spectral decomposition
        results = get_kappa_spectrum(A, Sigma_w);
        W_osc = results.Oscillatory.node_weights; 
        kappas = results.Oscillatory.kappas;      
        
        % Sort by frequency
        [~, freq_idx] = sort(results.Oscillatory.frequencies); 
        W_osc = W_osc(:, freq_idx);
        kappas = kappas(freq_idx);
        log_k = log(kappas'); 
        
        % Initialize subject container
        subjStruct = struct();
        subjStruct.SubjectID = i;
        
        % Process clusters
        for c = 1:n_clust
            clusterName = u_clusters(c);
            fieldName = matlab.lang.makeValidName(clusterName);
            node_idx = node_labels == clusterName;
            
            % Normalized energy distribution (eta)
            E_modes = sum(W_osc(node_idx, :), 1); 
            eta = E_modes / (sum(E_modes) + eps); 
            
            % Weighted mean (Average Ellipticity)
            ellip_mean = sum(eta .* log_k);
            
            % Weighted variance (Ellipticity Dispersion)
            ellip_var = sum(eta .* (log_k - ellip_mean).^2);
            
            % Store results
            subjStruct.(sprintf('%s_SumW', fieldName)) = ellip_mean;
            subjStruct.(sprintf('%s_VarW', fieldName)) = ellip_var;
        end
        
        rows_cell{i} = subjStruct;
    end
    
    % Final table generation
    T_final = struct2table([rows_cell{:}], 'AsArray', true);
    fprintf('Done. Dataset: %d subjects x %d features.\n', n_subj, width(T_final));
end