% clear workspace and set paths
clearvars; close all; clc

dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% list of subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

% unique networks & their indices
rname = fullfile(dataDir, 'regions.xlsx');
T     = readtable(rname);
netLabels = T.NETWORK;
[uniqueNets, ~, ic] = unique(netLabels);

% Euclidean distance matrix
tmp      = load(distFile,'mtx_euc_dis');
D        = tmp.mtx_euc_dis;
Dfull    = D + D.';         % make symmetric
n        = size(Dfull,1);
fitR    = [4.48, 12.18]; % select distance range for fitting

% define signed lags
max_lag = 50;
nLags   = 200;
tau_vals= linspace(0, max_lag, nLags);

% load model for selected subject
iSub    = 1;
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

%% SINGLE NETWORK TRIAL
netId = 5;
netName = uniqueNets{netId};
idx     = find(ic == netId);           % lines for this network
m       = numel(idx);              % number of nodes in this net

% build linear indices for the sub-block in Dfull and Ctens
[I, J]       = meshgrid(idx, idx);
linIdxFull   = sub2ind([n,n], I, J);

% extract distance vectors for upper and lower triangles of the sub-block
subMaskUp    = triu(true(m),  1);
subMaskLow   = tril(true(m), -1);
Dsub         = reshape(Dfull(linIdxFull), m, m);
dvals_up     = Dsub(subMaskUp);
dvals_low    = Dsub(subMaskLow);

% now threshold *these* for fitting
sel_up       = dvals_up  >= fitR(1) & dvals_up  <= fitR(2);
sel_low      = dvals_low >= fitR(1) & dvals_low <= fitR(2);

% preallocate slope arrays
a_up  = nan(1,nLags);
a_low = nan(1,nLags);

% loop over lags, slice Ctens into the same sub-block, and regress
for t = 1:nLags
    c_off = Ctens(:,:,t);
    Csub = reshape(c_off(linIdxFull), m, m);
    cu   = abs( Csub(subMaskUp) );
    cl   = abs( Csub(subMaskLow) );

    p1      = polyfit(log(dvals_up(sel_up)),   log(cu(sel_up)),   1);
    p2      = polyfit(log(dvals_low(sel_low)), log(cl(sel_low)), 1);

    a_up(t)  = p1(1);
    a_low(t) = p2(1);
end

% assemble signed-lag result
tau_neg  = -tau_vals(end:-1:1);
alpha    = [a_low(end:-1:1), a_up];
tau = [tau_neg, tau_vals];

% plot in next tile
figure
plot(tau, alpha);
title(netName, 'Interpreter','none');
xlabel('\tau');
ylabel('\alpha(\tau)');
grid on;

% plot zero-lag covariance over distance
C0_full   = Ctens(:,:,1);                  % first lag = zero
C0_sub    = reshape(C0_full(linIdxFull), m, m);
mask_up   = triu(true(m),1);               % same as subMaskUp
dvals_up  = Dsub(mask_up);                 % you already have this
y_all     = abs( C0_sub(mask_up) );        % 1×(m*(m-1)/2) vector

% plot
[d_sorted, sortIdx] = sort(dvals_up);
y_sorted            = y_all(sortIdx);

figure;
plot(log(d_sorted), log(y_sorted));
hold on;
xline(log(fitR(1)),'r--');
xline(log(fitR(2)),'r--');
xlabel('log(distance)');
ylabel('log(|covariance|)');
title(['Zero-lag Covariance vs Distance — ' netName], 'Interpreter','none');
grid on;

%% FITTING FOR EACH NETWORK
% set up tiled layout for all networks
nNets = numel(uniqueNets);
nCols = ceil(sqrt(nNets));
nRows = ceil(nNets/nCols);
figure('Units','normalized','Position',[.1 .1 .8 .8]);
tiledlayout(nRows, nCols, 'Padding','compact','TileSpacing','compact');

results = struct();
for k = 1:numel(uniqueNets)
    netName = uniqueNets{k};
    idx     = find(ic == k);           % lines for this network
    m       = numel(idx);              % number of nodes in this net
    
    % build linear indices for the sub-block in Dfull and Ctens
    [I, J]       = meshgrid(idx, idx);
    linIdxFull   = sub2ind([n,n], I, J);

    % extract distance vectors for upper and lower triangles *of the sub-block*
    subMaskUp    = triu(true(m),  1);
    subMaskLow   = tril(true(m), -1);
    Dsub         = reshape(Dfull(linIdxFull), m, m);
    dvals_up     = Dsub(subMaskUp);
    dvals_low    = Dsub(subMaskLow);

    % now threshold *these* for fitting
    sel_up       = dvals_up  >= fitR(1) & dvals_up  <= fitR(2);
    sel_low      = dvals_low >= fitR(1) & dvals_low <= fitR(2);

    % preallocate slope arrays
    a_up  = nan(1,nLags);
    a_low = nan(1,nLags);

    % loop over lags, slice Ctens into the same sub-block, and regress
    for t = 1:nLags
        c_off = Ctens(:,:,t);
        Csub = reshape(c_off(linIdxFull), m, m);
        cu   = abs( Csub(subMaskUp) );
        cl   = abs( Csub(subMaskLow) );

        p1      = polyfit(log(dvals_up(sel_up)),   log(cu(sel_up)),   1);
        p2      = polyfit(log(dvals_low(sel_low)), log(cl(sel_low)), 1);

        a_up(t)  = p1(1);
        a_low(t) = p2(1);
    end
    
    % assemble signed-lag result
    tau_neg  = -tau_vals(end:-1:1);
    alpha    = [a_low(end:-1:1), a_up];
    tau = [tau_neg, tau_vals];

    % plot in next tile
    nexttile;
    plot(tau, alpha);
    title(netName, 'Interpreter','none');
    xlabel('\tau');
    ylabel('\alpha(\tau)');
    grid on;
end
