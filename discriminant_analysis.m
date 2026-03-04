clear; clc;

% Configuration
dataDir = fullfile(pwd, 'data');
level = 'Macro_Area'; 
axisFontSize = 19;
labelFontSize = 23;
cmap = magma(256);
marker_color = cmap(200, :); 

% Data loading
fprintf('Processing dataset...\n');
T_raw = create_dataset(dataDir, level);
T_raw.OriginalIndex = (1:height(T_raw))'; 

% Feature extraction
featureNames = T_raw.Properties.VariableNames;
isFeature = ~strcmp(featureNames, 'SubjectID') & ~strcmp(featureNames, 'OriginalIndex');
X_raw = table2array(T_raw(:, isFeature));
feat_names = featureNames(isFeature);

% Quality control
idx_Var = contains(feat_names, '_VarW');
X_Var_Raw = X_raw(:, idx_Var);
global_dispersion = mean(X_Var_Raw, 2);
z_scores = (global_dispersion - mean(global_dispersion)) / std(global_dispersion);
threshold = 2.5; 
is_outlier = abs(z_scores) > threshold;

if any(is_outlier)
    fprintf('Removed %d outlier(s) based on global dispersion.\n', sum(is_outlier));
    T = T_raw(~is_outlier, :);
else
    T = T_raw;
end

% Final data preparation
X = table2array(T(:, isFeature));
X(isnan(X)) = 0; 

idx_Var = contains(feat_names, '_VarW');
X_Sigma = X(:, idx_Var); 
names_Sigma = feat_names(idx_Var);

idx_Mu = contains(feat_names, '_SumW');
X_Mu = X(:, idx_Mu); 
names_Mu = feat_names(idx_Mu);

% Statistical profiling
mu_X = mean(X, 1);
sigma_X = std(X, 0, 1);
cv_X = abs(sigma_X ./ (mu_X + eps));

% PCA on ellipticity dispersion
X_Sigma_z = zscore(X_Sigma);
[coeff, score, ~, ~, explained] = pca(X_Sigma_z);

% Enforce positive loading convention
if mean(coeff(:,1)) < 0
    coeff(:,1) = -coeff(:,1);
    score(:,1) = -score(:,1);
end

% Sort regions by loading weight
clean_names = erase(names_Sigma, {'_VarW', '_Var'});
[sorted_regions, idx_alpha] = sort(clean_names);
sorted_coeffs = coeff(idx_alpha, 1);

% Hypothesis testing
regions_list = unique(clean_names);
cv_invariants = zeros(length(regions_list), 1);
cv_discriminants = zeros(length(regions_list), 1);

for i = 1:length(regions_list)
    id_mu = find(strcmp(feat_names, [regions_list{i} '_SumW']));
    id_sig = find(strcmp(feat_names, [regions_list{i} '_VarW']));
    cv_invariants(i) = cv_X(id_mu);
    cv_discriminants(i) = cv_X(id_sig);
end

[p_val, h_val, stats] = signrank(cv_discriminants, cv_invariants, 'tail', 'right');
stat_val = stats.signedrank;

% Reporting
fprintf('\n--- ANALYSIS REPORT (N=%d) ---\n', height(T));

[sorted_cv_inv, idx_inv] = sort(cv_invariants, 'ascend');
fprintf('\nA. Architectural Invariants (Low CV)\n');
fprintf('%-20s | %-8s | %-20s\n', 'Region', 'CV', 'Mean +/- SD');
for k = 1:length(regions_list)
    reg_idx = idx_inv(k);
    feat_idx = find(strcmp(feat_names, [regions_list{reg_idx} '_SumW']));
    fprintf('%-20s | %.3f    | %.3f +/- %.3f\n', ...
        regions_list{reg_idx}, sorted_cv_inv(k), mu_X(feat_idx), sigma_X(feat_idx));
end

[sorted_cv_disc, idx_disc] = sort(cv_discriminants, 'descend');
fprintf('\nB. Individual Discriminants (High CV)\n');
fprintf('%-20s | %-8s | %-20s\n', 'Region', 'CV', 'Mean +/- SD');
for k = 1:length(regions_list)
    reg_idx = idx_disc(k);
    feat_idx = find(strcmp(feat_names, [regions_list{reg_idx} '_VarW']));
    fprintf('%-20s | %.3f    | %.3f +/- %.3f\n', ...
        regions_list{reg_idx}, sorted_cv_disc(k), mu_X(feat_idx), sigma_X(feat_idx));
end

fprintf('\nC. Statistical Validation\n');
fprintf('Wilcoxon Signed-Rank: W = %.1f, p = %.5e\n', stat_val, p_val);
fprintf('Metric Independence: r = %.3f\n', corr(mean(X_Mu,2), mean(X_Sigma,2)));
fprintf('PC1 Variance Explained: %.2f%%\n', explained(1));

%% Visualization
% Figure 1: PC1 Loadings
figure('Color', 'w', 'Position', [100 100 1100 600]);
bar(sorted_coeffs, 'FaceColor', marker_color, 'EdgeColor', 'k', 'FaceAlpha', 0.6);

set(gca, 'FontSize', axisFontSize);
ylabel('PC1 Loading', 'Interpreter', 'latex', 'FontSize', labelFontSize);
xticks(1:length(sorted_regions));
xticklabels(sorted_regions);
xtickangle(45);
xlim([0 length(sorted_regions)+1]);
grid on; box on;

% Figure 2: Metrics Overview
figure('Color', 'w', 'Position', [90 90 1000 500]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Scatter: Stability vs Diversity
nexttile;
scatter(mean(X_Mu,2), mean(X_Sigma,2), 80, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.6);
xlabel('Global Average Ellipticity ($\mu$)', 'Interpreter', 'latex', 'FontSize', 12); 
ylabel('Global Ellipticity Dispersion ($\sigma^2$)', 'Interpreter', 'latex', 'FontSize', 12);
title('Metric Independence', 'Interpreter', 'latex', 'FontSize', 14); 
grid on; axis square;

% PCA Score Plot
nexttile;
scatter(score(:,1), score(:,2), 80, 'filled', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'k');
text(score(:,1), score(:,2), string(T.OriginalIndex), 'FontSize', 9, 'VerticalAlignment','bottom');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)), 'FontSize', 12); 
ylabel(sprintf('PC2 (%.1f%%)', explained(2)), 'FontSize', 12);
title('Subject Embedding', 'Interpreter', 'latex', 'FontSize', 14); 
grid on; axis square;