clear; clc;

% Configuration
dataDir = fullfile(pwd, 'data');
level = 'Macro_Area'; 

% 1. Data Loading
T = create_dataset_v2(dataDir, level);

% Pre-processing
featureNames = T.Properties.VariableNames;
isFeature = ~strcmp(featureNames, 'SubjectID');
X = table2array(T(:, isFeature));
feat_names = featureNames(isFeature);

% Isolate Metric Groups (Weighted Variance vs Weighted Sum)
idx_Var = contains(feat_names, '_VarW');
idx_Sum = contains(feat_names, '_SumW');

X_Var = X(:, idx_Var); names_Var = feat_names(idx_Var);
X_Sum = X(:, idx_Sum); names_Sum = feat_names(idx_Sum);

% 2. Coefficient of Variation (CV) Profiling
mu_X = mean(X, 1);
sigma_X = std(X, 0, 1);
cv_X = abs(sigma_X ./ (mu_X + eps));

[sorted_cv, idx_cv] = sort(cv_X, 'descend');
sorted_names = feat_names(idx_cv);

% 3. PCA on Dynamical Diversity (Variance metrics)
% Analyzing latent structure of subject variability
X_Var_z = zscore(X_Var);
[coeff, score, ~, ~, explained] = pca(X_Var_z);

% Top drivers for PC1
[sorted_loadings, idx_load] = sort(abs(coeff(:,1)), 'descend');
top_drivers = names_Var(idx_load);

% 4. Global Coupling Analysis
% Correlating average tension (Sum) with dynamical diversity (Var)
global_stab = mean(X_Sum, 2);
global_div  = mean(X_Var, 2);
[r_coup, p_coup] = corr(global_stab, global_div);

%% --- Reporting ---
fprintf('\n--- DYNAMICAL ARCHITECTURE ANALYSIS: %s ---\n', upper(level));

fprintf('\nA. Low Variability Features (Invariants)\n');
fprintf('%-35s | %-10s | %-10s\n', 'Region', 'CV', 'Mean');
fprintf('%s\n', repmat('-',1,60));
for k = length(sorted_cv):-1:max(1, length(sorted_cv)-9)
    fprintf('%-35s | %.4f     | %.4f\n', sorted_names{k}, sorted_cv(k), mu_X(idx_cv(k)));
end

fprintf('\nB. High Variability Features (Discriminants)\n');
fprintf('%-35s | %-10s | %-10s\n', 'Region', 'CV', 'Mean');
fprintf('%s\n', repmat('-',1,60));
for k = 1:min(10, length(sorted_cv))
    fprintf('%-35s | %.4f     | %.4f\n', sorted_names{k}, sorted_cv(k), mu_X(idx_cv(k)));
end

fprintf('\nC. Latent Structure (PCA on Variance)\n');
fprintf('PC1 Explained Variance: %.2f%%\n', explained(1));
fprintf('PC2 Explained Variance: %.2f%%\n', explained(2));
fprintf('Top Drivers (PC1):\n');
for k = 1:5
    fprintf('   %d. %-30s (Loading: %.3f)\n', k, top_drivers{k}, sorted_loadings(k));
end

fprintf('\nD. Global Independence Check\n');
fprintf('Correlation (Mean Tension vs Diversity): r = %.4f (p = %.4f)\n', r_coup, p_coup);

%% --- Visualization ---
figure('Color','w', 'Position', [100 100 1200 400]);

% CV Profile
subplot(1,3,1); 
bar(sorted_cv); 
title('CV Profile'); 
ylabel('Coefficient of Variation'); 
xlim([0 length(sorted_cv)]); 
grid on; axis square;

% PCA Space
subplot(1,3,2);
scatter(score(:,1), score(:,2), 50, 'filled', 'MarkerFaceColor', 'k');
title(sprintf('PCA latent space (PC1: %.1f%%)', explained(1)));
xlabel('PC1 (Global Diversity)'); 
ylabel('PC2'); 
grid on; axis square;
text(score(:,1), score(:,2), string(1:size(score,1)), ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 8);

% Orthogonality Check
subplot(1,3,3);
scatter(global_stab, global_div, 50, 'filled', 'MarkerFaceColor', 'k');
xlabel('Global Mean Tension (\mu)'); 
ylabel('Global Diversity (\sigma^2)');
title('Metric Coupling'); 
grid on; axis square;
text(score(:,1), score(:,2), string(1:size(score,1)), ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 8);
