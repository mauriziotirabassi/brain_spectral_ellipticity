%% SETUP
clearvars, close all; clc

% load data
dataDir    = fullfile(pwd,'data');
distFile   = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir     = fullfile(dataDir,'regressed_001_01_sim62131');

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';
n        = size(D,1);
mask_upper     = triu(true(n),1);
dvals    = D(mask_upper); % array of distances i to j

%% TOY MODEL
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
Sigma = 20 * eye(n);
stds       = sqrt(diag(Sigma));
normMat    = stds * stds.'; % normalization matrix (Pearson correlation)

% TODO: play with different values of omega: see if making it scale with
% distance e.g. exp(-d) exhibits "travelling waves"

% toy differential cross-covariance: sequence of nodes
% omega = 1; % strength of rotation
% S = zeros(n);
% for i = 1:n-1
%     S(i, i+1) = omega;
%     S(i+1, i) = -omega;
% end

% omega scaling with distance
lambda = 20;
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

%% TIME-LAGGED COVARIANCE
% time-lagged covariance
Sigma_tau = @(t) expm(A * t) * Sigma;

% simulation hyperparameters (evolution of time-lagged covariance)
max_lag  = 100;
nLags    = 1001;
taus = linspace(0, max_lag, nLags);

% create array of time-lagged covariance matrices (includes tau=0)
Ctens       = nan(n,n,nLags);
for k = 1:nLags
    t = taus(k);
    Ctens(:,:,k) = Sigma_tau(t) ./ normMat;
    imshow(Ctens(:,:,k), [], 'InitialMagnification', 'fit')
    title(sprintf('Time-Lagged Covariance \\Sigma(\\tau) at lag \\tau = %.2f', taus(k)));
    drawnow
    pause(0.02)
end

%% SPATIO-TEMPORAL COVARIANCE SCALING
% fitting range and indices
[d_sorted, idx] = sort(dvals);
% fitR     = [2.90, 14.75];   % optimal from regression
fitR     = [4.48, 12.18];   % Benozzo
% fitR     = [8.13 33.82];    % Deco
fit_sel      = d_sorted >= fitR(1) & d_sorted <= fitR(2);
x = log(d_sorted);
xf = x(fit_sel);

% preallocate slope arrays
slopes_upper = nan(1, nLags);
slopes_lower = nan(1, nLags);

% mask for lower triangles
mask_lower = tril(true(n),-1);

% axes limits
xlims = log([min(d_sorted), max(d_sorted)]);
ylims = [-15, 1];

% differential cross-covariance threshold (given that those pairs do not
% contribute to detecting travelling-wave-like behavior)
threshold = 0.1;

% log-log plot over distance
figure;
% tiledlayout(2, 1, 'TileSpacing','compact','Padding','compact');
for k = 1:nLags
    % extract correlation matrix at lag k
    Ck = Ctens(:,:,k);
    
    % UPPER TRIANGLE
    % nexttile(1)
    cv_u = Ck(mask_upper);
    cv_su = cv_u(idx);

    % filter out pairs with null differential cross-covariance
    Su      = S(mask_upper);
    Su_s    = Su(idx);
    % valid_u = Su_s ~= 0;
    valid_u = abs(Su_s) >= threshold;
    x_u = x(valid_u);
    y_u = log(abs(cv_su(valid_u)));
    % x_u = x;
    % y_u = log(abs(cv_su));

    % plot + fit line
    plot(x_u, y_u); hold on;
    
    sel_u = x_u >= log(fitR(1)) & x_u <= log(fitR(2));
    xf_u    = x_u(sel_u);
    yf_u    = y_u(sel_u);
    p  = polyfit(xf_u, yf_u, 1);
    y_fit_u = polyval(p, xf_u);

    plot(xf_u, y_fit_u, 'r-', 'LineWidth', 1.5);
    title(sprintf('\\tau = %.2f', taus(k)));
    xlim(xlims); ylim(ylims); grid on;
    hold off;

    % save triu slope
    slopes_upper(k) = p(1);

    % % LOWER TRIANGLE
    % nexttile(2)
    % cv_l = Ck(mask_lower);
    % cv_sl = cv_l(idx);
    % Sl      = S(mask_lower);
    % Sl_s    = Sl(idx);
    % valid_l = Sl_s ~= 0;
    % x_l     = x(valid_l);
    % y_l = log(abs(cv_sl(valid_l)));
    % 
    % plot(x_l, y_l); hold on;
    % 
    % sel_l    = x_l >= log(fitR(1)) & x_l <= log(fitR(2));
    % xf_l     = x_l(sel_l);
    % yf_l     = y_l(sel_l);
    % p_l      = polyfit(xf_l, yf_l, 1);
    % y_fit_l  = polyval(p_l, xf_l);
    % 
    % plot(xf_l, y_fit_l, 'b-', 'LineWidth', 1.5);
    % xlim(xlims); ylim(ylims); grid on;
    % hold off;
    % 
    % % save tril slope
    % slopes_lower(k) = p_l(1);
    
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
plot(taus, slopes_upper, 'LineWidth', 1);
grid on; xlabel('\tau'); ylabel('slope');
title('Spatio-Temporal Scaling Exponent vs. Lag');
