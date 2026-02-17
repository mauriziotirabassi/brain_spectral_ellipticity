clear; clc;

% Configuration
dataDir = fullfile(pwd, 'data');
level = 'Macro_Area'; 

% --- 1. DATA GENERATION ---
fprintf('--- LOADING DATA ---\n');
T_raw = create_dataset_v2(dataDir, level);

% Extract just the numeric data for QC
featureNames = T_raw.Properties.VariableNames;
isFeature = ~strcmp(featureNames, 'SubjectID');
X_raw = table2array(T_raw(:, isFeature));
feat_names = featureNames(isFeature);

% Isolate Variance Metrics for QC
idx_Var = contains(feat_names, '_VarW');
X_Var_Raw = X_raw(:, idx_Var);

% --- 2. QUALITY CONTROL (Outlier Removal) ---
fprintf('\n--- QUALITY CONTROL ---\n');
global_div_score = mean(X_Var_Raw, 2);

% Calculate Z-scores
mu_global = mean(global_div_score);
sigma_global = std(global_div_score);
z_scores = (global_div_score - mu_global) / sigma_global;

% Threshold (e.g., > 2.5 Standard Deviations)
threshold = 2.5; 
is_outlier = abs(z_scores) > threshold;

if any(is_outlier)
    fprintf('Detected %d outlier(s) with Global Variance Z-score > %.1f:\n', sum(is_outlier), threshold);
    disp(T_raw.SubjectID(is_outlier));
    
    % Apply Filter
    T = T_raw(~is_outlier, :);
    fprintf('Removed outliers. Final N = %d.\n', height(T));
else
    fprintf('No outliers detected. Final N = %d.\n', height(T_raw));
    T = T_raw;
end

% --- 3. PRE-PROCESSING (On Clean Data) ---
% Re-extract X from the clean table T
X = table2array(T(:, isFeature));

% Clean NaNs/Infs
X(isnan(X)) = 0; X(isinf(X)) = 0;

% Re-Isolate Groups
X_Var = X(:, idx_Var); names_Var = feat_names(idx_Var);
idx_Sum = contains(feat_names, '_SumW');
X_Sum = X(:, idx_Sum); names_Sum = feat_names(idx_Sum);

% --- 4. STATISTICAL PROFILING (CV) ---
mu_X = mean(X, 1);
sigma_X = std(X, 0, 1);
cv_X = abs(sigma_X ./ (mu_X + eps));

[sorted_cv, idx_cv] = sort(cv_X, 'descend');
sorted_names = feat_names(idx_cv);

% --- 5. PCA on DIVERSITY (Variance Metrics) ---
X_Var_z = zscore(X_Var);
[coeff, score, ~, ~, explained] = pca(X_Var_z);

% Top Drivers of PC1
[sorted_loadings, idx_load] = sort(abs(coeff(:,1)), 'descend');
top_drivers = names_Var(idx_load);

% --- 6. INDEPENDENCE CHECK ---
global_stab = mean(X_Sum, 2);
global_div  = mean(X_Var, 2);
[r_coup, p_coup] = corr(global_stab, global_div);

% --- 7. REPORTING ---
fprintf('\n--- RESULTS: %s (N=%d) ---\n', upper(level), height(T));

fprintf('\nA. Invariants (Low CV)\n');
fprintf('%-35s | %-10s | %-20s\n', 'Region', 'CV', 'Mean +/- SD');
fprintf('%s\n', repmat('-',1,70));
for k = length(sorted_cv):-1:max(1, length(sorted_cv)-9)
    idx = idx_cv(k);
    fprintf('%-35s | %.4f     | %.4f +/- %.4f\n', ...
        sorted_names{k}, sorted_cv(k), mu_X(idx), sigma_X(idx));
end

fprintf('\nB. Discriminants (High CV)\n');
fprintf('%-35s | %-10s | %-20s\n', 'Region', 'CV', 'Mean +/- SD');
fprintf('%s\n', repmat('-',1,70));
for k = 1:min(10, length(sorted_cv))
    idx = idx_cv(k);
    fprintf('%-35s | %.4f     | %.4f +/- %.4f\n', ...
        sorted_names{k}, sorted_cv(k), mu_X(idx), sigma_X(idx));
end

fprintf('\nC. Latent Structure\n');
fprintf('PC1 Variance: %.2f%%\n', explained(1));
fprintf('PC2 Variance: %.2f%%\n', explained(2));
fprintf('Top Drivers (PC1):\n');
for k = 1:5
    fprintf('   %d. %-30s (Loading: %.3f)\n', k, top_drivers{k}, sorted_loadings(k));
end

fprintf('\nD. Independence Check\n');
fprintf('Correlation (Sum vs Var): r = %.4f (p = %.4f)\n', r_coup, p_coup);

% --- 8. VISUALIZATION ---
figure('Color','w', 'Position', [100 100 1200 400]);

% CV Profile
subplot(1,3,1); 
bar(sorted_cv); title('CV Profile'); 
ylabel('CV'); xlim([0 length(sorted_cv)]); grid on; axis square;

% PCA
subplot(1,3,2);
scatter(score(:,1), score(:,2), 50, 'filled', 'MarkerFaceColor', 'k');
text(score(:,1), score(:,2), string(1:size(score,1)), 'FontSize', 8, 'VerticalAlignment','bottom');
title(sprintf('PCA (PC1: %.1f%%)', explained(1)));
xlabel('PC1'); ylabel('PC2'); grid on; axis square;

% Coupling
subplot(1,3,3);
scatter(global_stab, global_div, 50, 'filled', 'MarkerFaceColor', 'k');
xlabel('Mean Tension (\mu)'); ylabel('Diversity (\sigma^2)');
title('Metric Independence'); grid on; axis square;