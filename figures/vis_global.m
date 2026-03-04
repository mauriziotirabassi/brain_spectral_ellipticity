clear all; clc
genFontSize = 20; % Global font size

% Configuration
dataDir = fullfile(pwd, 'data');
outDir = fullfile(dataDir, 'regressed_001_01_sim62131');

% Data processing
d = dir(fullfile(outDir, '*.mat'));
files = {d.name};
num_subs = length(files);
C_kappas = cell(num_subs, 1);
C_freqs = cell(num_subs, 1);
C_subIDs = cell(num_subs, 1);

fprintf('Processing %d subjects...\n', num_subs);

for i = 1:num_subs
    % Load data
    subj = load(fullfile(outDir, files{i}));
    A = subj.A;
    n = size(A, 1);
    
    % Spectral decomposition
    Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar; 
    results = get_kappa_spectrum(A, Sigma_w);
    kappas = results.Oscillatory.kappas;
    freqs = results.Oscillatory.frequencies;
    
    % Collect results
    C_kappas{i} = kappas;
    C_freqs{i} = freqs;
    C_subIDs{i} = repmat(i, length(kappas), 1);
end

% Flatten data structures
all_kappas = vertcat(C_kappas{:});
all_freqs = vertcat(C_freqs{:});
all_subIDs = vertcat(C_subIDs{:});

% Visualization
figure('Color', 'w', 'Position', [100 100 1100 600]);

% Scatter plot colored by subject ID
scatter(all_freqs, all_kappas, 75, all_subIDs, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);

colormap(magma);
set(gca, 'YScale', 'log');
set(gca, 'FontSize', genFontSize);

xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', genFontSize);
ylabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', genFontSize);

grid on; box on; 
axis tight;
ylim([min(all_kappas)*0.9 max(all_kappas)*1.2]); 
xlim([0 .39]);