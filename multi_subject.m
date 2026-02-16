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
    kappas  = results.Oscillatory.kappas;
    freqs   = results.Oscillatory.frequencies;
    
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
scatter(all_freqs, log(all_kappas), 50, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
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

%% --- Visualization: Global Frequency-Stability Constraints ---
figure('Color','w', 'Position', [100 100 600 500]);

% 1. Scatter Plot with Transparency
% Dark grey [0.2 0.2 0.2] with alpha 0.3 to reveal density
s = scatter(all_freqs, all_kappas, 75, [0.2 0.2 0.2], ...
    'filled', 'MarkerFaceAlpha', 0.35);

% 2. Axis Configuration
set(gca, 'YScale', 'log');           % Logarithmic Y-axis
set(gca, 'XScale', 'linear');        % Linear X-axis
axis tight;                          % Snap limits to data
ylim([1 max(all_kappas)*1.5]);       % Ensure Y starts at 1 (theoretical min)

% 3. Aesthetics & Grid
grid on;                             % Enable grid
box on;                              % Enable bounding box
set(gca, 'Layer', 'top');            % Places ticks and grid lines *above* the patches
set(gca, 'TickDir', 'out');          % Ticks point OUTWARD (prevents data overlap)
set(gca, 'LineWidth', 1.2);          % Publication-quality line weight
set(gca, 'FontSize', 12);            % Readable font size (approx 10-12pt)
set(gca, 'FontName', 'Helvetica');   % Standard sans-serif font

% 4. Labels (Using LaTeX)
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Non-Normality Factor ($\kappa$)', 'Interpreter', 'latex', 'FontSize', 14);

% NO TITLE: The title will be provided in the LaTeX figure caption.