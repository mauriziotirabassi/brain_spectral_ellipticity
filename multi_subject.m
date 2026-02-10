clear all; clc

% --- Setup ---
dataDir = fullfile(pwd,'data');
outDir  = fullfile(dataDir,'regressed_001_01_sim62131');
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};
num_subs = length(files);

% --- Preallocation ---
C_kappas = cell(num_subs, 1);
C_freqs  = cell(num_subs, 1);
C_subIDs = cell(num_subs, 1);

fprintf('Processing %d subjects...\n', num_subs);

% --- Processing Loop ---
for i = 1:num_subs
    % Load Data
    subj = load(fullfile(outDir, files{i}));
    A = subj.A;
    n = size(A, 1);
    Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar; 

    % Compute Spectrum
    results = get_kappa_spectrum(A, Sigma_w);
    kappas  = results.kappas;
    freqs   = results.frequencies;
    
    % Collect Data
    C_kappas{i} = kappas;
    C_freqs{i}  = freqs;
    C_subIDs{i} = repmat(i, length(kappas), 1);
end

% --- Flatten Data ---
all_kappas = vertcat(C_kappas{:});
all_freqs  = vertcat(C_freqs{:});
all_subIDs = vertcat(C_subIDs{:});

%% --- Visualization ---
figure('Color','w', 'Position', [50 100 1200 500]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Tile 1: Global Scatter (Log Kappa vs Frequency)
nexttile;
scatter(all_freqs, log(all_kappas), 15, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Log Non-Normality ($\log \kappa$)', 'Interpreter', 'latex');
title(sprintf('Global Spectrum (All Modes, N=%d)', num_subs));
grid on; axis square;

% Tile 2: Boxplot of Kappa Distribution per Subject
nexttile;
boxplot(log(all_kappas), all_subIDs);
xlabel('Subject ID');
ylabel('Log Non-Normality ($\log \kappa$)', 'Interpreter', 'latex');
title('Inter-Subject Variability');
grid on; 

sgtitle('Multi-Subject Spectral Analysis');