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

%% SOURCE/SINK MEASURES
net_flow = sum(S,1);     % sum columns: net outflow per region (sources positive, sinks negative)

% source/sink strengths
source_strength = max(net_flow, 0); 
sink_strength   = max(-net_flow, 0);

% plot
figure;
bar(net_flow);
xlabel('Brain region'); ylabel('Net flow (positive = source)');
title('Net source/sink flow per region (columns of S)');
grid on;

% map net_flow to colors; normalize net_flow to [1,256]
cm = parula(256);
flow_norm = (net_flow - min(net_flow)) / (max(net_flow) - min(net_flow));
cidx      = round(1 + flow_norm * 255);
node_colors = cm(cidx,:);

% plot
figure('Color','w');
axis equal off; view(3); hold on
Y  = cmdscale(Dfull,3);
scatter3(Y(:,1), Y(:,2), Y(:,3), 36, [0.8 0.8 0.8], 'filled'); % gray nodes
scatter3(Y(:,1), Y(:,2), Y(:,3), 60, node_colors, 'filled');   % overlay colored nodes
colormap(cm);
hcb = colorbar; %  sources (red) vs sinks (blue)
hcb.Label.String = 'Net flow (pos = source, neg = sink)';

%% NESS MEASURES (ENTROPY PRODUCTION)
phi = -2 * trace(A*S/Q);     % entropy production rate (scale of irreversibility)
S_norm = norm(S,'fro');      % Frobenius norm = sqrt(sum_ij S_ij^2)
node_irrev = sum(abs(S),1);  % node-level irreversibility (sum of |flows|)