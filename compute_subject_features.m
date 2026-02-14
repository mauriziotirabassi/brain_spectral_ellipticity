function feature_matrix = compute_subject_features(clusterCol, dataDir, outDir)
% COMPUTE_SUBJECT_FEATURES  Constructs the raw feature matrix for all subjects
%
% INPUTS:
%   clusterCol : string, name of column in names.xlsx to use for clusters 
%                e.g., 'Region', 'Macro', 'Network'
%   dataDir    : string, path to folder containing names.xlsx
%   outDir     : string, path to folder containing subject .mat files
%
% OUTPUT:
%   feature_matrix : (nSubjects x nFeatures) raw feature matrix
%
% USAGE:
%   feature_matrix = compute_subject_features('Region', 'data', 'data/regressed_001_01_sim62131');

    %% --- Load subject files ---
    d = dir(fullfile(outDir,'*.mat')); 
    files = {d.name};
    nSubj = length(files);

    % --- Load cluster definitions ---
    T_table = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
    node_labels = string(T_table.(clusterCol));
    u_clusters = unique(node_labels);
    n_clust = length(u_clusters);

    % --- Initialize feature storage ---
    % O_C + PR_C per cluster + PR_eta mean+std + PR_rho mean+std
    feature_matrix = zeros(nSubj, n_clust*2 + 4);

    for iSub = 1:nSubj
        fprintf('Processing subject %d/%d\n', iSub, nSubj);

        % --- Load subject ---
        subj = load(fullfile(outDir,files{iSub}));
        A = subj.A;
        n = size(A,1);
        Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar;

        % --- Run spectral decomposition ---
        results = get_kappa_spectrum(A, Sigma_w);

        % Oscillatory matrices
        W_osc_mat = results.Oscillatory.node_weights;   % N x M_osc
        W_real_mat = results.Real.node_weights;         % N x M_real
        kappas = results.Oscillatory.kappas;

        W_osc_total = sum(W_osc_mat, 2); 
        W_real_total = sum(W_real_mat, 2);
        W_total_reconstructed = W_osc_total + W_real_total;

        %% --- 1. Regional Oscillatory Fraction O_C and PR_C ---
        O_C = zeros(n_clust,1);
        PR_C = zeros(n_clust,1);
        E_C_modes = zeros(n_clust, length(kappas));

        for i = 1:n_clust
            idx = node_labels == u_clusters(i);

            cluster_osc_energy = sum(W_osc_total(idx));
            cluster_total_energy = sum(W_total_reconstructed(idx));
            O_C(i) = cluster_osc_energy / cluster_total_energy;

            mode_energies = sum(W_osc_mat(idx, :), 1);
            E_C_modes(i, :) = mode_energies;

            % Participation ratio per cluster
            num = sum(mode_energies)^2;
            den = length(mode_energies) * sum(mode_energies.^2);
            PR_C(i) = num / den;
        end

        %% --- 2. η_C: column-wise participation ratio ---
        eta_mat = E_C_modes ./ sum(E_C_modes,1); % clusters x modes
        pr_eta_cols = sum(eta_mat).^2 ./ sum(eta_mat.^2);  % vector: PR per mode
        PR_eta_mean = mean(pr_eta_cols);
        PR_eta_std = std(pr_eta_cols);

        %% --- 3. ρ_C: row-wise participation ratio ---
        rho_mat = E_C_modes ./ sum(E_C_modes,2); % clusters x modes
        pr_rho_rows = sum(rho_mat,2).^2 ./ sum(rho_mat.^2,2);  % vector: PR per cluster
        PR_rho_mean = mean(pr_rho_rows);
        PR_rho_std = std(pr_rho_rows);

        %% --- 4. Assemble feature vector ---
        feat_vec = [O_C(:); PR_C(:); PR_eta_mean; PR_eta_std; PR_rho_mean; PR_rho_std];
        feature_matrix(iSub,:) = feat_vec';
    end

end
