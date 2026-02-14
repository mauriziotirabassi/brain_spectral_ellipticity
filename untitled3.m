% ANALYZE_FEATURES
% Statistical profiling of spectral features to identify Invariants vs. Discriminators.

clear all; clc; % close all;

% --- 1. Settings & Data Generation ---
dataDir = fullfile(pwd, 'data');
clusterCol = 'Network'; % Adjust if your column name is different

disp('--- Generating Dataset ---');
T_subject_params = create_dataset(dataDir, clusterCol);

% Extract Numeric Matrix (exclude SubjectID)
% Assumes first column is ID, rest are features
feature_names = T_subject_params.Properties.VariableNames(2:end);
X = table2array(T_subject_params(:, 2:end)); 
[n_subj, n_features] = size(X);

% --- 2. Statistical Profiling (The Filters) ---
mu = mean(X, 1);
sigma = std(X, 0, 1);

% Metric 1: Coefficient of Variation (CV)
% We use absolute value of mean to handle potentially negative metrics (though rare here)
CV = sigma ./ abs(mu); 

% Sort features by CV
[CV_sorted, sort_idx] = sort(CV, 'ascend');
sorted_names = feature_names(sort_idx);

% Identify Invariants (Bottom 10%) and Discriminators (Top 10%)
n_top = round(0.10 * n_features);
if n_top < 5, n_top = 5; end % Ensure at least 5

invariants = sorted_names(1:n_top);
discriminators = sorted_names(end-n_top+1:end);

% --- 3. Redundancy Analysis (Correlation) ---
% Compute Correlation Matrix
R = corr(X);

% Find Highly Correlated Pairs (|r| > 0.95)
% We only look at the upper triangle to avoid duplicates
mask = triu(true(n_features), 1); 
[row, col] = find(abs(R) > 0.95 & mask);
redundant_pairs = [feature_names(row)', feature_names(col)'];

% --- 4. Dimensionality Analysis (PCA) ---
% Standardize Data (Z-score) before PCA
X_z = zscore(X);
[coeff, score, latent, tsquared, explained, mu_pca] = pca(X_z);

% Cumulative Variance
cum_var = cumsum(explained);
n_pcs_90 = find(cum_var >= 90, 1);

% Generalized Variance (Determinant of Covariance)
gen_var = det(cov(X_z));

%% --- 5. Visualization ---
figure('Color', 'w', 'Position', [100 100 1400 500]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact');

% Plot 1: Stability Spectrum (CV)
nexttile;
bar(CV_sorted, 'FaceColor', [0.2 0.4 0.6]);
yline(0.1, 'r--', 'Threshold (10%)', 'LineWidth', 1.5);
xlabel('Feature Rank');
ylabel('Coefficient of Variation (CV)');
title('Stability Spectrum: Invariants vs. Discriminators');
grid on;
xlim([0 n_features+1]);

% Plot 2: Redundancy Map (Correlation)
nexttile;
imagesc(abs(R));
colormap(gca, 'parula');
c = colorbar;
c.Label.String = '|Correlation|';
clim([0 1]);
xlabel('Features'); ylabel('Features');
title('Redundancy Map');
axis square;

% Plot 3: Subject Space (PCA)
nexttile;
scatter(score(:,1), score(:,2), 80, 'filled', 'MarkerFaceColor', [0.8 0.2 0.2]);
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title(sprintf('Subject Space (90%% Var = %d PCs)', n_pcs_90));
grid on; axis square;
% Add labels for subjects (optional, if IDs are short)
text(score(:,1)+0.1, score(:,2), string(T_subject_params.SubjectID), 'FontSize', 8, 'Interpreter', 'none');

% --- 6. Output Report ---
disp(' ');
disp('========================================');
disp('       FEATURE EVALUATION REPORT        ');
disp('========================================');
fprintf('Total Features: %d\n', n_features);
fprintf('Effective Dimensionality (90%% Var): %d PCs\n', n_pcs_90);
fprintf('Generalized Variance (Determinant): %.2e\n', gen_var);

disp(' ');
disp('--- TOP 5 INVARIANTS (Most Stable) ---');
disp(table(CV_sorted(1:5)', sorted_names(1:5)', 'VariableNames', {'CV', 'Feature'}));

disp(' ');
disp('--- TOP 5 DISCRIMINATORS (Most Variable) ---');
disp(table(CV_sorted(end-4:end)', sorted_names(end-4:end)', 'VariableNames', {'CV', 'Feature'}));

disp(' ');
fprintf('Found %d redundant feature pairs (|r| > 0.95).\n', length(row));
if ~isempty(redundant_pairs)
    disp('Sample Redundancies:');
    disp(redundant_pairs(1:min(5, length(row)), :));
end
disp('========================================');
