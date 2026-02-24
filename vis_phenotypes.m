clear all; clc
% --- 1. CONFIGURATION ---
subjects_to_plot = [19, 10]; 
clusterCol = 'Macro_Area'; 
dataDir = fullfile(pwd, 'data');
outDir = fullfile(dataDir, 'regressed_001_01_sim62131');
T_names = readtable(fullfile(dataDir, 'names.xlsx'), 'VariableNamingRule', 'preserve');
u_clusters = unique(string(T_names.(clusterCol)));
n_clust = length(u_clusters);
% Graphics Settings
f_size_lbl = 18;  
f_size_ax = 18;   
label_font_size = 24; 
% Figure Size (Single Subject)
fig_w = 1100;
fig_h = 700; 
% --- 2. MAIN LOOP (ONE FIGURE PER SUBJECT) ---
for s_idx = 1:length(subjects_to_plot)
    iSub = subjects_to_plot(s_idx);
    
    % --- LOAD & PROCESS ---
    d = dir(fullfile(outDir, '*.mat'));
    subj = load(fullfile(outDir, d(iSub).name)); 
    results = get_kappa_spectrum(subj.A, eye(size(subj.A,1)) * subj.output.eff_conn.NoiseVar);
    
    freqs = results.Oscillatory.frequencies;
    kappas = results.Oscillatory.kappas;
    W_osc = results.Oscillatory.node_weights;
    W_real = results.Real.node_weights;
    [sorted_freqs, idx] = sort(freqs);
    sorted_logK = log(kappas(idx));
    sorted_W_osc = W_osc(:, idx);
    
    % Global Metrics
    global_E = sum(sorted_W_osc, 1);
    dots_size = (global_E / sum(W_osc(:))) * 3500; 
    dots_size(dots_size < 15) = 15; 
    dots_PR = (global_E.^2) ./ (size(subj.A,1) * sum(sorted_W_osc.^2, 1));
    
    % Regional Metrics
    O_C = zeros(n_clust, 1); PR_C = zeros(n_clust, 1);
    H_mat_Ell = zeros(n_clust, length(freqs));
    for i = 1:n_clust
        node_idx = string(T_names.(clusterCol)) == u_clusters(i);
        E_osc_reg = sum(sum(sorted_W_osc(node_idx, :)));
        E_real_reg = sum(sum(W_real(node_idx, :)));
        O_C(i) = E_osc_reg / (E_osc_reg + E_real_reg);
        E_modes = sum(sorted_W_osc(node_idx, :), 1);
        PR_C(i) = sum(E_modes)^2 / (length(E_modes) * sum(E_modes.^2));
        norm_E = E_modes / sum(E_modes);
        H_mat_Ell(i, :) = norm_E .* sorted_logK';
    end
    
    % --- 3. VISUALIZATION ---
    figure('Color','w', 'Position', [50+(s_idx*20), 50, fig_w, fig_h], 'Name', ['Subject ' num2str(iSub)]);
    t = tiledlayout(3, 3, 'TileSpacing','compact', 'Padding','compact');
    
    % === PANEL A: SCATTER (Row 1, Cols 1-2) ===
    ax1 = nexttile(1, [1 2]);
    scatter(sorted_freqs, sorted_logK, dots_size, dots_PR, 'filled', 'MarkerEdgeColor','k');
    ylabel('$\log \kappa$', 'Interpreter','latex', 'FontSize', f_size_lbl); 
    grid on; box on;
    colormap(ax1, magma);
    c1 = colorbar; 
    set(c1.Label, 'String', 'Spatial Reach', 'Interpreter', 'latex', 'FontSize', f_size_ax);
    clim(ax1, [0 1]);
    xlim([0 max(sorted_freqs)*1.05]);
    ylim([0 3.3]); 
    set(gca, 'FontSize', f_size_ax);
    
    % Label A
    if s_idx == 1
        text(-0.08, 1.15, '\textbf{A}', 'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', label_font_size);
    end
    % === PANEL C: BAR CHART (Rows 1-3, Col 3) ===
    ax3 = nexttile(3, [3 1]);
    b = barh(O_C);
    b.FaceColor = 'flat'; b.CData = PR_C;
    yticks(1:n_clust); yticklabels(u_clusters); 
    xlabel('$O_{\mathcal{C}}$', 'Interpreter','latex', 'FontSize', f_size_lbl);
    grid on; box on;
    colormap(ax3, magma);
    c3 = colorbar; 
    set(c3.Label, 'String', 'Spectral Participation', 'Interpreter', 'latex', 'FontSize', f_size_ax);
    clim(ax3, [0 1]);
    set(gca, 'FontSize', f_size_ax);
    
    if s_idx == 1
        text(-0.25, 1.05, '\textbf{C}', 'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', label_font_size);
    end
    % === PANEL B: HEATMAP (Rows 2-3, Cols 1-2) ===
    ax2 = nexttile(4, [2 2]);
    imagesc(H_mat_Ell); 
    set(gca, 'YDir', 'normal');
    yticks(1:n_clust); yticklabels(u_clusters);
    c2 = colorbar;
    set(c2.Label, 'String', '$\eta \log \kappa$', 'Interpreter', 'latex', 'FontSize', f_size_ax);
    colormap(ax2, magma);
    clim(ax2, [0 .25]);
    xlabel('Frequency (Hz)', 'Interpreter','latex', 'FontSize', f_size_lbl); 
    
    % Ticks
    [~, unique_idx] = unique(round(sorted_freqs, 2), 'stable');
    if length(unique_idx) > 10, unique_idx = unique_idx(1:2:end); end
    xticks(unique_idx); 
    xticklabels(string(round(sorted_freqs(unique_idx), 2)));
    xtickangle(45);
    set(gca, 'FontSize', f_size_ax);
    
    if s_idx == 1
        text(-0.08, 1.05, '\textbf{B}', 'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', label_font_size);
    end
    
    % Global Title applied to the layout
    title(t, ['\textbf{Subject ' num2str(iSub) '}'], 'Interpreter', 'latex', 'FontSize', 20);
end