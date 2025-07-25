%% SETUP
clearvars, close all; clc

% load data
dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
n        = size(D,1);
mask      = triu(true(n),1);
dvals    = D(mask); % array of distances i to j

%% SUBJECT DATA
% list of subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};
iSub = 13; % select subject number

% load this subject’s fitted A, TR, noise variance
subj = load(fullfile(outDir,files{iSub}));
% A = subj.A;
Sigma_w = eye(n)*subj.output.eff_conn.NoiseVar; % noise covariance
% Sigma = lyap(A,Sigma_w); % instantaneous stationary covariance
% S = 0.5 * (A * Sigma - Sigma * (A.')); % differential cross-covariance

% stds and normalization matrix (want Pearson correlation)
stds       = sqrt(diag(Sigma));
normMat    = stds * stds.';

%% TIME-LAGGED COVARIANCE
% time-lagged covariance
% Sigma_tau = @(t) expm(A*t)*Sigma;

% simulation hyperparameters
max_lag  = 100;
nLags    = 1001;
tau_vals = linspace(0,max_lag,nLags);

% create array of time-lagged covariance matrices (includes tau=0)
Ctens       = nan(n,n,nLags);  
for k = 1:nLags
    t = tau_vals(k);
    Ctens(:,:,k) = Sigma_tau(t) ./ normMat;
end

%% SPATIO-TEMPORAL COVARIANCE SCALING
% fitting range and indices
% fitR     = [2.90, 14.75];   % optimal from regression
fitR     = [4.48, 12.18];   % Benozzo
% fitR     = [8.13 33.82];    % Deco
fit_sel      = dvals >= fitR(1) & dvals <= fitR(2);

