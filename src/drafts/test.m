clear all, close all; clc

% Data directory
data = dir("data\");
data = data(~ismember({data.name}, {'.', '..'})); % Filter /. and /.. subfolders

% Structural
structural = load(fullfile(pwd, "data", data(1).name)).full_connectome_no_symm;
% la matrice strutturale va trasposta per essere coerente con le matrici 
% del modello
structural = structural';

% Distance
distance = load(fullfile(pwd, "data", data(3).name)).mtx_euc_dis;

% Directory batch 1
regressed1 = dir(fullfile(pwd, "data", data(4).name));
regressed1 = regressed1(~ismember({regressed1.name}, {'.', '..'}));

% First model batch 1
data1 = load(fullfile(pwd, "data", data(4).name, regressed1(1).name));
samples = 1:data1.N; em_iter = 1:size(data1.A_hat, 3);

%% VISUALIZE FITTING PROCESS
% fMRI volumetrico. parcellizzazione anatomica 74 basata su L'Allen
% (atlante famoso) per tradeoff computazionale modello e coprire modello.
% Voxel mediati stessa regione. TR sampling period 1s. w 1870 therefore
% 30m. Laboratorio Alessandro Gozzi.
figure, tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile(1), imshow(data1.A), title('A'), colormap("parula")
c = colorbar; c.Layout.Tile = 'south'; c.Label.String = 'Connection Strength';
for i = 1:data1.output.iter_em
    nexttile(1), imshow(data1.A_hat(:, :, i)), colormap("parula")
    title(['Effective Connectivity Estimate ' num2str(i)])
    nexttile(2), plot(data1.h_hat(:, :, i))
    title(['Haemodynamic Response Estimate ' num2str(i)])
    drawnow
end

% E-M estimated connection strength 1 \to 2
em_conn12 = squeeze(data1.A_hat(1, 2, :));
figure, plot(em_iter, em_conn12)
title('Estimated Connection Strength 1 to 2 Evolution')

%% MATRICES DEFINITION
global A Q S P P_tau
A = data1.A; % Effective connectivity
Q = eye(size(A)) * data1.output.eff_conn.NoiseVar; % Noise correlation
P = lyap(A, Q); % Zero-lag covariance
S = 0.5 * (A*P - P*A'); % dC-Cov
P_tau = @(tau) expm(A*tau) * P; % Lagged covariance

TR = data1.TR;

%% DETERMINING TIME INTERVAL FOR LAGGED COVARIANCE
time_int = 40; time_axis = linspace(0, time_int);
lag_cov = calc_lagged_covariance(time_int);
figure, tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
conn34 = squeeze(lag_cov(3, 4, :));
nexttile, plot(time_axis, conn34)
conn12 = squeeze(lag_cov(1, 2, :));
nexttile, plot(time_axis, conn12)

%% FUNCTIONS
function P_tau_values = calc_lagged_covariance(time)
    global P_tau
    tau_values = linspace(0, time);
    [rown, coln] = size(P_tau(0));
    P_tau_values = zeros(rown, coln, length(tau_values));
    for i = 1:length(tau_values)
        P_tau_values(:, :, i) = P_tau(tau_values(i));
    end
end