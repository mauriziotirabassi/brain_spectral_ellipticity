clear all; clc

% --- 1. CONFIGURATION ---
dataDir = fullfile(pwd,'data');
outDir  = fullfile(dataDir,'regressed_001_01_sim62131');
%%
% Graphics Settings (Matching Script 1)
f_size_lbl = 23;   
f_size_ax = 19;    
label_font_size = 24; 

% Prepare Color (Fixed value from Magma palette)
figure('Visible', 'off'); 
cmap = magma(256);
% Pick a vibrant color from the upper-mid range (e.g., orange/red)
marker_color = cmap(200, :); 
close;

%% --- 2. DATA PROCESSING ---
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};
num_subs = length(files);

C_kappas = cell(num_subs, 1);
C_freqs  = cell(num_subs, 1);
C_subIDs = cell(num_subs, 1);

fprintf('Processing %d subjects...\n', num_subs);

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

% Flatten Data
all_kappas = vertcat(C_kappas{:});
all_freqs  = vertcat(C_freqs{:});
all_subIDs = vertcat(C_subIDs{:});

%% --- 3. VISUALIZATION ---
%
% === Figure 1: Tiled Layout ===
% figure('Color','w', 'Position', [50 100 1200 600]);
% t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % Tile 1: Global Scatter (Kappa vs Frequency)
% ax1 = nexttile;
% % Style: Magma Color, Filled, Black Edges (matching Script 1 style)
% scatter(all_freqs, all_kappas, 50, marker_color, 'filled', ...
%     'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6); 
% 
% set(gca, 'YScale', 'log'); % Logarithmic Y-axis as requested
% xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
% ylabel('Non-Normality ($\kappa$)', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
% title(sprintf('Global Spectrum (N=%d)', num_subs), 'FontSize', f_size_ax);
% grid on; box on; axis square;
% set(gca, 'FontSize', f_size_ax); 
% 
% % Add Panel Label A
% text(-0.15, 1.05, '\textbf{A}', 'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', label_font_size);
% 
% % Tile 2: Boxplot (using Log Kappa for distribution visibility)
% ax2 = nexttile;
% boxplot(log(all_kappas), all_subIDs);
% xlabel('Subject ID', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
% ylabel('$\log \kappa$', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
% title('Inter-Subject Variability', 'FontSize', f_size_ax);
% grid on; box on; axis square;
% set(gca, 'FontSize', f_size_ax);
% 
% % Add Panel Label B
% text(-0.15, 1.05, '\textbf{B}', 'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', label_font_size);

%% === Figure 2: Standalone Global Scatter (Exact Copy of Tile 1) ===
figure('Color','w', 'Position', [100 100 800 600]);

% Scatter: Same Magma color, black edges, filled, transparency
scatter(all_freqs, all_kappas, 200, marker_color, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);

% Axis Configuration
set(gca, 'YScale', 'log');           % Logarithmic Y-axis
set(gca, 'FontSize', f_size_ax);     % Match font size
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', f_size_lbl);
ylabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', f_size_lbl);

grid on; box on; 
axis tight;
ylim([min(all_kappas)*0.9 max(all_kappas)*1.2]); % Adjust limits for log scale
xlim([0 .39])