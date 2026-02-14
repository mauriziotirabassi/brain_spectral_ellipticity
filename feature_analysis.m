clear all; clc

% --- Setup Parameters ---
dataDir = fullfile(pwd, 'data');
clusterCol = 'Network'; 

% --- 1. Generate Dataset ---
% Assumes create_dataset_v2 is available in the environment
T_final = create_dataset_v2(dataDir, clusterCol);

% --- 2. Feature Analysis ---
featureNames = T_final.Properties.VariableNames;
isFeature = ~strcmp(featureNames, 'SubjectID');
X = table2array(T_final(:, isFeature));
names = featureNames(isFeature);

% Handle potential NaNs or Infs
X(isnan(X)) = 0;
X(isinf(X)) = 0;

% A. Coefficient of Variation (CV = std/mean)
mu = mean(X, 1);
sigma = std(X, 0, 1);
cv_vals = abs(sigma ./ (mu + eps));
[sorted_cv, cv_idx] = sort(cv_vals, 'descend');
sorted_names = names(cv_idx);

% B. Correlation Analysis (Redundancy)
R = corr(X);
R_tri = triu(abs(R), 1); 
[r_vals, r_idx] = sort(R_tri(:), 'descend');
[row_idx, col_idx] = ind2sub(size(R), r_idx);

%% --- 3. Reporting ---
dashLine = repmat('-', 1, 50);
equalLine = repmat('=', 1, 50);

fprintf('\n%s\n', equalLine);
fprintf('FEATURE DISCRIMINATION (Top 10 High CV)\n');
fprintf('%s\n', dashLine);
for j = 1:min(10, length(names))
    fprintf('%-35s : %.4f\n', names{cv_idx(j)}, sorted_cv(j));
end

fprintf('\n%s\n', equalLine);
fprintf('FEATURE INVARIANCE (Top 10 Low CV)\n');
fprintf('%s\n', dashLine);
rev_cv_idx = flip(cv_idx);
sorted_cv_rev = flip(sorted_cv);
for j = 1:min(10, length(names))
    fprintf('%-35s : %.4f\n', names{rev_cv_idx(j)}, sorted_cv_rev(j));
end

fprintf('\n%s\n', equalLine);
fprintf('FEATURE REDUNDANCY (Top 10 Similar Pairs)\n');
fprintf('%s\n', dashLine);
for j = 1:min(10, length(r_vals))
    if r_vals(j) == 0, break; end
    fprintf('%s <-> %s : r = %.4f\n', names{row_idx(j)}, names{col_idx(j)}, r_vals(j));
end

% --- 4. Visualization ---
figure('Color','w', 'Position', [100 100 1100 600]);

% Subplot 1: Sorted CV Bar Chart
subplot(1, 2, 1);
bar(sorted_cv, 'FaceColor', [0.2 0.4 0.6]);
title('Feature Discrimination (Sorted CV)');
ylabel('Coefficient of Variation (std/mean)');
xticks(1:length(sorted_names));
xticklabels(sorted_names);
xtickangle(90);
set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 8);
grid on;

% Subplot 2: Correlation Matrix
subplot(1, 2, 2);
imagesc(abs(R));
title('Feature Redundancy (|Correlation Matrix|)');
xticks(1:length(names));
yticks(1:length(names));
xticklabels(names);
yticklabels(names);
xtickangle(90);
set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 8);
axis square;
colorbar;
colormap(magma);
