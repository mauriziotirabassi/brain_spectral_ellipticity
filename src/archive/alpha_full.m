% clear workspace and set paths
clearvars; close all; clc

dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% list of subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';         % make symmetric
n        = size(Dfull,1);

% masks for upper and lower triangles (ensuring alignment)
mask_up     = triu(true(n), 1);
mask_low    = tril(true(n), -1);
dvals_up    = Dfull(mask_up);
dvals_low   = Dfull(mask_low);

% select distance range for fitting
fitR    = [4.48, 12.18];
sel_up  = dvals_up >= fitR(1) & dvals_up <= fitR(2);
sel_low = dvals_low >= fitR(1) & dvals_low <= fitR(2);

% define signed lags
max_lag = 100;
nLags   = 1001;
tau_vals= linspace(0, max_lag, nLags);

% load model for selected subject
iSub    = 13;
subj    = load(fullfile(outDir,files{iSub}));
A       = subj.A;
Q       = eye(n) * subj.output.eff_conn.NoiseVar;
P       = lyap(A, Q);
P_tau   = @(t) expm(A*t)*P;
stds    = sqrt(diag(P));     % n×1 vector of standard deviations
normMat = stds * stds.';     % n×n matrix where (i,j)=σ_i*σ_j

% calculate lagged-covariance
Ctens   = nan(n,n,nLags); % preallocate
for k = 1:nLags
    t = tau_vals(k);
    C    = P_tau(t) ./ normMat;
    Ctens(:,:,k) = abs(C);
end

%% TIME-LAGGED COVARIANCE DISTANCE SCALING
% preallocate exponent vectors
a_vec_up       = nan(1,nLags);
a_vec_low      = nan(1,nLags);

% regression
for k = 1:nLags
    c_off     = Ctens(:,:,k);

    % fit upper triangular
    c_off_up     = c_off(mask_up);
    p_up         = polyfit(log(dvals_up(sel_up)), log(abs(c_off_up(sel_up))), 1);
    a_vec_up(k)  = p_up(1);

    % fit lower triangular
    c_off_low     = c_off(mask_low);
    p_low         = polyfit(log(dvals_low(sel_low)), log(abs(c_off_low(sel_low))), 1);
    a_vec_low(k)  = p_low(1);
end

tau_neg      = -tau_vals(end:-1:1);        % = -[3,2,1,0] = [-3,-2,-1,0]
tau   = [tau_neg, tau_vals];        % = [-3,-2,-1,0, 0,1,2,3]

a_neg        = a_vec_low(end:-1:1);        % = [4,3,2,1]
alpha = [a_neg, a_vec_up];          % = [4,3,2,1,7,8,9,10]

% plot
figure;
plot(tau, alpha);
xlabel('\tau');
ylabel('\alpha(\tau)');
title('Scaling exponent \alpha vs signed time-lag \tau');
grid on;

%% FREQUENCY CONTENT OF ALPHA
% % TODO: better spectrum estimation
% % sampling parameters
% M       = length(tau);
% dt      = tau(2) - tau(1);        % assumed uniform spacing
% Fs      = 1/dt;                   % sampling frequency in 1/units_of_tau
% 
% % FFT
% Y       = fft(alpha);
% P2      = abs(Y/M);               % two-sided spectrum (normalized)
% P1      = P2(1:floor(M/2)+1);     % one-sided spectrum
% P1(2:end-1) = 2*P1(2:end-1);      % account for energy in negative freqs
% 
% % frequency axis
% f       = Fs*(0:floor(M/2))/M;    % in cycles per unit τ
% 
% % plot
% figure, tiledlayout(2, 1)
% nexttile, plot(tau, alpha);
% xlabel('\tau'); ylabel('\alpha(\tau)');
% title('Original signal');
% nexttile, plot(f, P1);
% xlabel('Frequency (cycles per unit \tau)'); ylabel('|A(f)|');
% title('Single-sided amplitude spectrum of \alpha(\tau)');
% grid on;