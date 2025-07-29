% Toy model (manually defining S and Sigma) where only a subnetwork is
% considered in order to explore the assumption that one would not be able
% to simulate a perfectly flat slope vs. lag plot.

%% SETUP
clearvars, close all; clc

% load data
dataDir    = fullfile(pwd,'data');
distFile   = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir     = fullfile(dataDir,'regressed_001_01_sim62131');

% unique networks & their indices
rname = fullfile(dataDir, 'regions.xlsx');
T     = readtable(rname);
netLabels = T.NETWORK;
[uniqueNets, ~, ic] = unique(netLabels);

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';
n        = size(D,1);

% simulation variables
global lambda sigma;
lambda = 20; sigma = 2;

% single network selection & distance matrix pruning
netId = 13;
netName = uniqueNets{netId};
idx     = find(ic == netId);       % lines for this network
m       = numel(idx);              % number of nodes in this net
[I, J]       = meshgrid(idx, idx);
linIdxFull   = sub2ind([n,n], I, J);
Dsub      = reshape(Dfull(linIdxFull), m, m);
mask_upper   = triu(true(m),  1);
dvals = Dsub(mask_upper);

%% TOY MODEL (WHOLE-BRAIN)
% extracting noise covariance out of one subject
d     = dir(fullfile(outDir,'*.mat')); % list of subject‐model .mat files
files = {d.name};
iSub  = 13; % select subject number
subj  = load(fullfile(outDir,files{iSub}));
% A = subj.A;
Sigma_w = eye(n) * subj.output.eff_conn.NoiseVar; % noise covariance
% Sigma = lyap(A,Sigma_w); % instantaneous stationary covariance
% S = 0.5 * (A * Sigma - Sigma * (A.')); % differential cross-covariance

% starting with regions totally uncorrelated (want to see how correlations
% pop up due to curl flow)
Sigma = sigma * eye(n);
stds       = sqrt(diag(Sigma));
normMat    = stds * stds.'; % normalization matrix (Pearson correlation)

% toy differential cross-covariance: sequence of nodes
% omega = 1; % strength of rotation
% S = zeros(n);
% for i = 1:n-1
%     S(i, i+1) = omega;
%     S(i+1, i) = -omega;
% end

% omega scaling with distance
S = zeros(n);
for i = 1:n-1
    for j = i+1:n
        omega_ij = exp(-Dfull(i,j) / lambda);
        S(i, j) = omega_ij;
        S(j, i) = -omega_ij;
    end
end

% effective connectivity parametrized by (Sigma, S)
A = (-0.5 * Sigma_w + S) / Sigma;

%% TIME-LAGGED COVARIANCE (SINGLE NETWORK)
% time-lagged covariance
Sigma_tau = @(t) expm(A * t) * Sigma;

% simulation hyperparameters (evolution of time-lagged covariance)
max_lag  = 50;
nLags    = 500;
taus = linspace(0, max_lag, nLags);

% create array of time-lagged covariance matrices (includes tau=0)
Csub       = nan(m,m,nLags);
% figure; tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
for k = 1:nLags
    t = taus(k);
    Ctens = Sigma_tau(t) ./ normMat;
    Csub(:,:,k) = reshape(Ctens(linIdxFull), m, m);

    % plot
    % nexttile(1), imshow(Ctens, [])
    % title(sprintf('Time-Lagged Covariance \\Sigma(\\tau) at lag \\tau = %.2f', taus(k)));
    % nexttile(2), imshow(Csub(:,:,k), [])
    % title(netName, 'Interpreter','none')
    % drawnow
    % pause(0.01)
end

%% SPATIO-TEMPORAL COVARIANCE SCALING
close all
% fitting range and indices
[d_sorted, d_idx] = sort(dvals);
% fitR     = [2.90, 14.75];   % optimal from regression
fitR     = [max(4.48, min(dvals)), min(12.18, max(dvals))];   % Benozzo
% fitR     = [8.13 33.82];    % Deco
fit_sel      = d_sorted >= fitR(1) & d_sorted <= fitR(2);
x = log(d_sorted);
xf = x(fit_sel);

% preallocate slope arrays
slopes = nan(1, nLags);

% axes limits
xlims = log([min(d_sorted), max(d_sorted)]);
ylims = [-15, 1];

% differential cross-covariance threshold (given that those pairs do not
% contribute to detecting travelling-wave-like behavior)
threshold = 0.5;

% log-log plot over distance
% figure;
for k = 1:nLags
    % extract correlation matrix at lag k
    Ck = Csub(:,:,k);
    
    % UPPER TRIANGLE
    % nexttile(1)
    cv_u = Ck(mask_upper);
    cv_su = cv_u(d_idx);
    x_u = x;
    y_u = log(abs(cv_su));

    % filter out pairs with null differential cross-covariance
    Su      = S(mask_upper);
    Su_s    = Su(d_idx);
    % valid_u = Su_s ~= 0;
    valid = abs(Su_s) >= threshold;
    % x_u = x(valid);
    % y_u = log(abs(cv_su(valid)));

    % fit line    
    sel_u = x_u >= log(fitR(1)) & x_u <= log(fitR(2));
    xf    = x_u(sel_u);
    yf    = y_u(sel_u);
    p  = polyfit(xf, yf, 1);
    y_fit_u = polyval(p, xf);
    slopes(k) = p(1); % save triu slope

    % plot
    plot(x_u, y_u); hold on;
    plot(xf, y_fit_u, 'r-', 'LineWidth', 1.5);
    title(sprintf('\\tau = %.2f', taus(k)));
    xlim(xlims); ylim(ylims); grid on;
    hold off;

    drawnow;
    pause(0.02);
end

%% SLOPE VS LAG PLOT
% given the time invariance of the time-lagged covariance
% taus_neg   = -flip(taus);
% slopes_neg = flip(slopes_lower);
% taus_all   = [taus_neg, taus];
% slopes_all = [slopes_neg, slopes_upper];

% plot
figure('Color','w');
% plot(taus_all, slopes_all, 'LineWidth', 1);
plot(taus, slopes, 'LineWidth', 1);
grid on; xlabel('\tau'); ylabel('slope');
title('Spatio-Temporal Scaling Exponent vs. Lag');
