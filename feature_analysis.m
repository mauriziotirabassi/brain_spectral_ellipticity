clear; clc;

% Configuration
dataDir = fullfile(pwd, 'data');
level = 'Macro_Area'; 

% --- 1. DATA GENERATION ---
fprintf('--- LOADING DATA ---\n');
T_raw = create_dataset_v2(dataDir, level);

% Create persistent index column
T_raw.OriginalIndex = (1:height(T_raw))'; 

% Extract numeric data
featureNames = T_raw.Properties.VariableNames;
isFeature = ~strcmp(featureNames, 'SubjectID') & ~strcmp(featureNames, 'OriginalIndex');
X_raw = table2array(T_raw(:, isFeature));
feat_names = featureNames(isFeature);

% Isolate Variance Metrics for QC
idx_Var = contains(feat_names, '_VarW');
X_Var_Raw = X_raw(:, idx_Var);

% --- 2. QUALITY CONTROL (Outlier Removal) ---
fprintf('\n--- QUALITY CONTROL ---\n');
global_div_score = mean(X_Var_Raw, 2);
mu_global = mean(global_div_score);
sigma_global = std(global_div_score);
z_scores = (global_div_score - mu_global) / sigma_global;

threshold = 2.5; 
is_outlier = abs(z_scores) > threshold;

if any(is_outlier)
    fprintf('Detected %d outlier(s) > %.1f SD:\n', sum(is_outlier), threshold);
    disp(T_raw(is_outlier, {'OriginalIndex', 'SubjectID'}));
    T = T_raw(~is_outlier, :);
    fprintf('Removed outliers. Final N = %d.\n', height(T));
else
    fprintf('No outliers detected. Final N = %d.\n', height(T_raw));
    T = T_raw;
end

% --- 3. PRE-PROCESSING ---
X = table2array(T(:, isFeature));
X(isnan(X)) = 0; X(isinf(X)) = 0;

% Groups
idx_Var = contains(feat_names, '_VarW');
X_Var = X(:, idx_Var); names_Var = feat_names(idx_Var);
idx_Sum = contains(feat_names, '_SumW');
X_Sum = X(:, idx_Sum); names_Sum = feat_names(idx_Sum);

% --- 4. STATISTICAL PROFILING (CV) ---
mu_X = mean(X, 1);
sigma_X = std(X, 0, 1);
cv_X = abs(sigma_X ./ (mu_X + eps));

[sorted_cv, idx_cv] = sort(cv_X, 'descend');
sorted_names = feat_names(idx_cv);

% --- 5. PCA on DIVERSITY ---
X_Var_z = zscore(X_Var);
[coeff, score, ~, ~, explained] = pca(X_Var_z);
[sorted_loadings, idx_load] = sort(abs(coeff(:,1)), 'descend');
top_drivers = names_Var(idx_load);

% --- 6. INDEPENDENCE CHECK ---
global_stab = mean(X_Sum, 2);
global_div  = mean(X_Var, 2);
[r_coup, p_coup] = corr(global_stab, global_div);

% --- 6b. HYPOTHESIS TESTING (One-Tailed Wilcoxon) ---
% Extract paired CVs
regions = unique(erase(feat_names, {'_SumW', '_VarW'}));
cv_avg = []; cv_disp = [];

for i = 1:length(regions)
    idx_s = find(strcmp(feat_names, [regions{i} '_SumW']));
    idx_v = find(strcmp(feat_names, [regions{i} '_VarW']));
    if ~isempty(idx_s) && ~isempty(idx_v)
        cv_avg(end+1) = cv_X(idx_s);
        cv_disp(end+1) = cv_X(idx_v);
    end
end

% Test if Discriminants (Var) > Invariants (Sum)
[p_wilc, h_wilc, stats_wilc] = signrank(cv_disp, cv_avg, 'tail', 'right');

% Handle statistic name (zval vs signedrank): can change with clustering
if isfield(stats_wilc, 'zval')
    stat_val = stats_wilc.zval; stat_name = 'Z';
else
    stat_val = stats_wilc.signedrank; stat_name = 'W';
end

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

fprintf('\nC. Statistical Validation (CV_Disp > CV_Avg)\n');
fprintf('Test: One-tailed Wilcoxon Signed-Rank\n');
fprintf('Statistic: %s = %.2f\n', stat_name, stat_val);
fprintf('P-value:   %.5e\n', p_wilc);
fprintf('Result:    %s\n', string(h_wilc == 1));

fprintf('\nD. Latent Structure\n');
fprintf('PC1 Variance: %.2f%%\n', explained(1));
clean_names = erase(names_Var, {'_VarW', '_Var'});
[sorted_regions, idx_alpha] = sort(clean_names);
sorted_coeffs = coeff(idx_alpha, 1);
fprintf('Region Loadings (PC1 - Alphabetical):\n');
for k = 1:length(sorted_regions)
    fprintf('   %d. %-30s (Loading: %.3f)\n', k, sorted_regions{k}, sorted_coeffs(k));
end

fprintf('\nE. Independence Check\n');
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
text(score(:,1), score(:,2), string(T.OriginalIndex), ...
    'FontSize', 8, 'VerticalAlignment','bottom');
title(sprintf('PCA (PC1: %.1f%%)', explained(1)));
xlabel('PC1'); ylabel('PC2'); grid on; axis square;

% Coupling
subplot(1,3,3);
scatter(global_stab, global_div, 50, 'filled', 'MarkerFaceColor', 'k');
xlabel('Mean Tension (\mu)'); ylabel('Diversity (\sigma^2)');
title('Metric Independence'); grid on; axis square;

%% --- 9. NEW FIGURE: LOADINGS BARCHART (STYLED) ---
% Graphics Settings
f_size_lbl = 23;   
f_size_ax = 19;    

% Set figure width to 12 inches to fit screen
figure('Units', 'inches', 'Position', [1 1 12 6], 'Color', 'w');

cmap = magma(256);
marker_color = cmap(200, :); 

% Create Bar Chart
b = bar(sorted_coeffs, 'FaceColor', marker_color, 'EdgeColor', 'k', 'FaceAlpha', 0.6);

% Axis Configuration
set(gca, 'FontSize', f_size_ax);
set(gca, 'XDir', 'reverse'); % Flip X Axis direction

% Labels & Ticks
ylabel('PC1 Loading', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
xticks(1:length(sorted_regions));
xticklabels(sorted_regions);
xtickangle(45); % Rotate 45 degrees

% Final Touches
grid on; box on;
title(''); % No title as requested
xlim([0 length(sorted_regions)+1]);