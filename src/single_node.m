%% SETUP
clearvars; close all; clc

dataDir   = fullfile(pwd,'data');
distFile  = fullfile(dataDir,'mtx_euc_distance.mat');
outDir    = fullfile(dataDir,'regressed_001_01_sim62131');

% choose your subject and node
iSub   = 1;     % which .mat in outDir
iNode  = 70;    % which row of P to inspect

% lag parameters
max_lag = 40;
nLags   = 100;
tau_vals= linspace(0,max_lag,nLags);

%% PRELOAD
% load Euclidean distances
tmp    = load(distFile,'mtx_euc_dis');
D      = tmp.mtx_euc_dis;
Dfull  = D + D.';       % make it symmetric
n      = size(Dfull,1); % number of nodes

% extract only the distances from iNode \to all j\neq iNode
dist_node = Dfull(iNode,:).';
dist_node(iNode) = []; % drop self‑distance, now dist_node is (n‑1)x1

% bin the distance vector
logd    = log10(dist_node);
IQR     = prctile(logd,75) - prctile(logd,25);
bw      = 2 * IQR / numel(logd)^(1/3); 
nbins   = ceil((max(logd)-min(logd))/bw); % Freedman–Diaconis 
edges   = linspace(min(dist_node),max(dist_node), nbins+1);
centers = (edges(1:end-1)+edges(2:end))/2;

% load one subject
flist = dir(fullfile(outDir,'*.mat'));
S     = load(fullfile(outDir,flist(iSub).name));
A     = S.A;
Q     = eye(n)*S.output.eff_conn.NoiseVar;
Sigma = lyap(A,Q);
P_tau = @(t) expm(A*t)*Sigma;

% compute node-to-all covariances over time
node_cov = nan(n-1, nLags);
for k = 1:nLags
    Pt           = P_tau(tau_vals(k));
    vec_all      = Pt(iNode,:).';    % n×1
    vec_all(iNode) = [];             % drop diagonal
    node_cov(:,k)= vec_all;          % store
end
% now I have a matrix where each row contains the time series of lagged
% covariance with respect to another node.

%% BIN & FIT SLOPE FOR EACH \tau
% averaging over binned distance assumes concentric propagation
slopes = nan(1,nLags);
for k = 1:nLags
    cvec = node_cov(:,k);

    % bin-average
    Bbin = nan(1,nbins);
    for b = 1:nbins
        sel     = dist_node>=edges(b) & dist_node<edges(b+1);
        Bbin(b) = mean(cvec(sel));
    end

    % fit within Deco's [8.13 33.82] mm range
    fitR   = [8.13 33.82];
    idx    = centers>=fitR(1) & centers<=fitR(2);
    x_fit  = centers(idx);
    y_fit  = abs(Bbin(idx));
    p      = polyfit(log10(x_fit), log10(y_fit), 1);
    slopes(k)= p(1);

    % plot
    clf; loglog(centers, abs(Bbin), 'b.-'); hold on;
    loglog(x_fit, 10.^(polyval(p, log10(x_fit))), 'k--', 'LineWidth', 1.5);
    ylim([1e-10, 1e-3]), xlim([0, 70]) ; grid on;
    xlabel('Distance (mm)'); ylabel('|Covariance|');
    title(sprintf('\\tau = %.2f,  slope \\alpha = %.3f', tau_vals(k), slopes(k)));
    legend('binned data','power-law fit', 'Location', 'northeast');
    drawnow;
end

% plot
figure, plot(tau_vals, slopes)
xlabel('\tau'), ylabel('Spatial decay exponent \alpha')
title(sprintf('Subject %d, Node %d', iSub, iNode)); grid on;

%% VISUALIZE NODE‐TO‐NODE COVARIANCES OVER SPACE
% recover 3D coordinates by classical MDS on the distance matrix
Y = cmdscale(Dfull, 3); % Y is n×3: [X Y Z] coordinates

% prepare diverging colormap
nC = 256; cmap = parula(nC);

% normalize the node_cov to enhance visibility across \tau
norm_cov = (node_cov - mean(node_cov, 2)) ./ std(node_cov, 0, 2);  % z-score by node
vmax = max(abs(norm_cov(:)));
vmin = -vmax;

% create figure
figure('Color', 'w');
axis equal off
view(3); hold on
xlabel('X'), ylabel('Y'), zlabel('Z');
sgtitle(sprintf('Lagged covariance of node %d', iNode));

% preplot all nodes as gray dots
scatter3(Y(:,1), Y(:,2), Y(:,3), 36, [.7 .7 .7], 'filled');
scatter3(Y(iNode,1), Y(iNode,2), Y(iNode,3), 100, 'k','filled'); % highlight iNode

% add colorbar
cb = colorbar;
colormap(cmap);
clim([vmin vmax]);
cb.Label.String = 'Z-scored covariance from source node';

% loop through \tau and plot connections
for k = 1:nLags
    % grab normalized covariances from node iNode to others at tau_k
    zvals = norm_cov(:,k);  % already z-scored
    zvals = [zvals(1:iNode-1); NaN; zvals(iNode:end)]; % insert NaN at self

    % map z-scores to colormap indices
    ci = round((zvals - vmin) / (vmax - vmin) * (nC - 1)) + 1;
    ci = min(max(ci, 1), nC);  % clip out-of-range indices

    % plot lines from iNode to all other nodes with color
    for j = 1:n
        if j == iNode || isnan(ci(j)), continue; end
        col = cmap(ci(j), :);
        line([Y(iNode,1), Y(j,1)], ...
             [Y(iNode,2), Y(j,2)], ...
             [Y(iNode,3), Y(j,3)], ...
             'Color', col, 'LineWidth', 1);
    end

    % Position dynamic title in a fixed spot
    txt = text(min(Y(:,1)), max(Y(:,2)) + 5, max(Y(:,3)), ...
        sprintf('\\tau = %.2f', tau_vals(k)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    drawnow; pause(0.1); delete(txt);

    if k < nLags
        delete(findobj(gca, 'Type', 'line'));
    end
end
hold off;