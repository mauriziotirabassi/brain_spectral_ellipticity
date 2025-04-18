clear all, close all; clc
global D A Q S P P_tau BINS

% Data directory
data = dir("data\");
data = data(~ismember({data.name}, {'.', '..'})); % Filter /. and /.. subfolders

% Structural (transposed for coherence with model matrices)
structural = load(fullfile(pwd, "data", data(1).name)).full_connectome_no_symm;
structural = structural';

% Directory batch 1
regressed1 = dir(fullfile(pwd, "data", data(4).name));
regressed1 = regressed1(~ismember({regressed1.name}, {'.', '..'}));

% First model batch 1
data1 = load(fullfile(pwd, "data", data(4).name, regressed1(1).name));
samples = 1:data1.N; em_iter = 1:size(data1.A_hat, 3);

%% MATRICES DEFINITION
D = load(fullfile(pwd, "data", data(3).name)).mtx_euc_dis; % Distance
A = data1.A; % Effective connectivity
Q = eye(size(A)) * data1.output.eff_conn.NoiseVar; % Noise covariance
P = lyap(A, Q); % Zero-lag covariance
S = 0.5 * (A*P - P*(A')); % dC-Cov
P_tau = @(tau) (expm(A*tau) * P); % Lagged covariance

TR = data1.TR;

%% MAX TIME LAG ESTIMATION
max_lag = 50; num_lags = 50;
tau_values = linspace(0, max_lag, num_lags);
% TODO: 0 is zero-lag covariance, don't know what the smallest one could be though

P_matrices = lagcov(tau_values);
figure, tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
time_axis = linspace(0, max_lag, num_lags);
for i = 1:size(P_matrices, 1)
    for j = 1:size(P_matrices, 2)
        conn = squeeze(P_matrices(i, j, :));
        plot(time_axis, conn), hold on
    end
end
xlabel('\tau'), ylabel('P_{ij}(\tau)'), title('Covariance evolution \forall(i,j)')
% TODO: Decided upon 45

%% DISTANCE BINNING (CLUSTERING)
n = size(A, 1); mask = ~eye(n); % True for all i\neq j
D_full = D + D'; % Making distances symmetric
BINS.dvals  = D_full(mask); % d_ij for every ordered pair

% Calculate number of bins i.e. distance clusters
logd = log10(BINS.dvals);
IQR = prctile(logd, 75) - prctile(logd, 25); % Interquartile range
bin_width = 2 * (IQR / numel(logd)^(1/3)); % Freedman-Diaconis bin width
range_data = max(logd) - min(logd); % Range of distance values
BINS.num = ceil(range_data / bin_width); % Number of bins

% The distance value corresponding to a distance bin is interpolated as
% the geometric mean of the boundaries of each bin

% Linearly binned
BINS.edges = linspace(min(BINS.dvals), max(BINS.dvals), BINS.num+1);
BINS.centers = 0.5 * (BINS.edges(1:end-1) + BINS.edges(2:end));

% % Logarithmically binned
% edges = logspace(log10(min(dvals)), log10(max(dvals)), num_bins+1);
% centers = sqrt(edges(1:end-1) .* edges(2:end));

%% COVARIANCE OVER DISTANCE
% Array of discrete evolution of covariance matrix
max_lag = 45; num_lags = 100;
tau_values = linspace(0, max_lag, num_lags);
P_matrices = lagcov(tau_values);

slopes = nan(1, num_lags); figure
for i = 1:num_lags
    B = distcov(P_matrices(:, :, i));

    % TODO: negative values ignored, gotta rescale?
    loglog(BINS.centers, abs(B)); % Log-log space plot
    xlabel('log(r)'), ylabel('log|B(r;\tau)|')
    % xlim([0 100]), ylim([0 1e-3])
    title(['Lagged covariance at \tau=' int2str(tau_values(i))]);

    % Range linear fitting
    fit_range = [8.13, 33.82];  % Same as Deco paper
    fit_idx = BINS.centers >= fit_range(1) & BINS.centers <= fit_range(2);
    p = polyfit(log10(BINS.centers(fit_idx)), log10(abs(B(fit_idx))), 1);
    slopes(i) = p(1); intercept = p(2);
    fit_line = 10.^(intercept + slopes(i)*log10(BINS.centers(fit_idx)));
    
    % Plot fit line on top
    hold on, loglog(BINS.centers(fit_idx), fit_line, 'k--');
    legend('binned data', 'power-law fit');
    drawnow, %hold off
end

% Slope evolution over time
figure, plot(tau_values, slopes)
xlabel('\tau'), ylabel('\alpha(\tau)')
title('Evolution of spatial decay exponent \alpha with \tau')

%% FUNCTIONS

function P_tau_values = lagcov(tau_values)
% LAGCOV returns an array of time-delayed covariance matrices calculated 
% for each value in tau_values
    global P_tau
    [rown, coln] = size(P_tau(0));
    P_tau_values = zeros(rown, coln, length(tau_values));
    for i = 1:length(tau_values)
        P_tau_values(:, :, i) = P_tau(tau_values(i));
    end
end

function B = distcov(covariance)
% DISTCOV returns the empirical covariance over distance function output
% for a time-lagged covariance matrix covariance
    global A BINS
    n = size(A, 1); mask = ~eye(n); % True for all i\neq j
    cvals  = covariance(mask); % P_ij(tau) for the same ordered pairs
    % Array with place correspondance of distance and covariance
    
    % Covariance average in each bin
    B = nan(1, BINS.num);
    for k = 1:BINS.num
        sel  = BINS.dvals >= BINS.edges(k) & BINS.dvals < BINS.edges(k+1); % Indices of values within bin
        B(k) = mean(cvals(sel)); % Mean covariance for that distance‐bin
    end
end