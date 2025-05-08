%% SETUP
clearvars, close all; clc

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
% mask     = ~eye(n);
mask      = triu(true(n),1); % working with upper triangular
dvals    = Dfull(mask);

% fitting range and indices
fitR     = [2.90, 14.75];   % optimal from regression
% fitR     = [4.48, 12.18];   % Benozzo
% fitR     = [8.13 33.82];    % Deco
sel      = dvals >= fitR(1) & dvals <= fitR(2);

% simulation hyperparameters
max_lag  = 75;
nLags    = 200;
tau_vals = linspace(0,max_lag,nLags);

%% SINGLE SUBJECT DATA
% select subject number
iSub = 1;

% load this subject’s fitted A, TR, noise variance
subj = load(fullfile(outDir,files{iSub}));
A = subj.A;
Q = eye(n)*subj.output.eff_conn.NoiseVar;
P = lyap(A,Q);
S = 0.5 * (A * P - P * (A.'));
P_tau = @(t) expm(A*t)*P;

% symmetric and skew-symmetric generator decomposition
A_sym   = -0.5 * (Q / P);
A_skew  = S / P + inv(P) * S; % TODO: correct inv

% lagged dynamics
P_tau_sym  = @(t) expm(A_sym * t)  * P;
P_tau_skew = @(t) expm(A_skew * t) * P;

% stds and normalization matrix
stds       = sqrt(diag(P));     % n×1 vector of standard deviations
normMat    = stds * stds.';     % n×n matrix where (i,j)=σ_i*σ_j

% preallocate
Ctens       = nan(n,n,nLags);  
Ctens_sym   = nan(n,n,nLags);
Ctens_skew  = nan(n,n,nLags);

% % abs two points per distance value
% for k = 1:nLags
%     t = tau_vals(k);
%     Ctens(:,:,k)      = P_tau(t)      ./ normMat; % full covariance
%     Ctens_sym(:,:,k)  = P_tau_sym(t)  ./ normMat; % symmetric part
%     Ctens_skew(:,:,k) = P_tau_skew(t) ./ normMat; % skew-symmetric part
% end

% abs max between two symmetric parts
for k = 1:nLags
    t = tau_vals(k);

    % full
    C    = P_tau(t)      ./ normMat;
    C    = max(C, C.');         % elementwise max of C_ij and C_ji
    Ctens(:,:,k) = abs(C);

    % symmetric part
    Csym = P_tau_sym(t)  ./ normMat;
    Csym = max(Csym, Csym.');
    Ctens_sym(:,:,k) = abs(Csym);

    % skew‐symmetric part
    Csk  = P_tau_skew(t) ./ normMat;
    Csk  = max(Csk, Csk.');
    Ctens_skew(:,:,k) = abs(Csk);
end

% TODO: seaprate upper and lower triu and plot the lower triu values as
% values of negative tau given that property.

%% ZERO-LAG COVARIANCE DISTANCE SCALING
y_all           = Ctens(:,:,1);
y_all           = y_all(mask);
[d_sorted, idx] = sort(dvals);
y_sorted        = y_all(idx);
xmin            = log(fitR(1));
xmax            = log(fitR(2));

figure, plot(log(d_sorted), log(abs(y_sorted))), hold on;
xline(xmin, 'r--'); xline(xmax, 'r--');
xlabel('log(distance)'); ylabel('log(|covariance|)');
title('Zero-lag Covariance over Distance'); grid on;

%% TIME-LAGGED COVARIANCE DISTANCE SCALING
% preallocate exponent vectors
a_vec       = nan(1,nLags);
a_sym_vec   = nan(1,nLags);
a_skew_vec  = nan(1,nLags);

% regression
for k = 1:nLags
    % FULL dynamics
    c_off     = Ctens(:,:,k);
    c_off     = c_off(mask);
    p         = polyfit(log(dvals(sel)), log(abs(c_off(sel))), 1);
    a_vec(k)  = p(1);

    % SYMMETRIC dynamics
    c_off_sym      = Ctens_sym(:,:,k);
    c_off_sym      = c_off_sym(mask);
    p_sym          = polyfit(log(dvals(sel)), log(abs(c_off_sym(sel))), 1);
    a_sym_vec(k)   = p_sym(1);

    % SKEW‐SYMMETRIC dynamics
    c_off_skew     = Ctens_skew(:,:,k);
    c_off_skew     = c_off_skew(mask);
    p_skew         = polyfit(log(dvals(sel)), log(abs(c_off_skew(sel))), 1);
    a_skew_vec(k)  = p_skew(1);
end

% plot the three together
figure;
plot(tau_vals, a_vec, 'k'); hold on
plot(tau_vals, a_sym_vec, 'b--');
plot(tau_vals, a_skew_vec, 'r--');
legend('full A', 'symmetric part', 'skew-symmetric part');
xlabel('\tau'); ylabel('Slope a(\tau)');
title(sprintf('Distance-decay exponent vs lag in subject %d', iSub));
grid on;