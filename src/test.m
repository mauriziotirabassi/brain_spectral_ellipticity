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

%% A PARAMETRIZATION
global A  S
A = data1.A;
Q = eye(size(A)) * data1.output.eff_conn.NoiseVar;
TR = data1.TR;
Sigma = lyap(A, Q);
S = 0.5 * (A*Sigma - Sigma*A');

%% PLAYING WITH TIME LAG
conn_num = 3;%size(S, 1);
time_int = 50;
figure, tiledlayout(conn_num, conn_num, 'TileSpacing', 'compact', 'Padding', 'compact')
for i = 1:conn_num
    for j = 1:conn_num
        nexttile, plot(lagged_covariance_values(i, j, time_int)); grid on;
    end
end

%% PLAYING WITH DISTANCE


%% FUNCTIONS
function S_tau_values = lagged_covariance_values(row, col, time)
    global A S
    lagged_covariance = @(tau) expm(A*tau) * S;
    tau_values = linspace(0, time);
    S_tau_values = zeros(length(tau_values), 1);
    for i = 1:length(tau_values)
        S_tau_matrix = lagged_covariance(tau_values(i));
        S_tau_values(i) = S_tau_matrix(row, col);
    end
end