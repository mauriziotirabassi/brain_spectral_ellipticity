clearvars; close all; clc;

%% ————— PATHS & FILES —————
dataDir = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% get list of your subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

%% ————— LOAD DISTANCES —————
fprintf('Loading distance matrix from:\n  %s\n', distFile);
tmp    = load(distFile,'mtx_euc_dis');
Dfull  = tmp.mtx_euc_dis + tmp.mtx_euc_dis.';    % symmetric
n      = size(Dfull,1);
mask   = ~eye(n);
dvals  = Dfull(mask);
fprintf('n = %d, extracted %d distance‐pairs\n', n, numel(dvals));

%% ————— BUILD Ctens(:,:,:) —————
iSub = 1;
fprintf('Loading subject %d from file:\n  %s\n', iSub, files{iSub});
subj = load(fullfile(outDir,files{iSub}));
A    = subj.A;
Q    = eye(n)*subj.output.eff_conn.NoiseVar;
P    = lyap(A,Q);       % Zero-lag covariance
stds = sqrt(diag(P));
normMat = stds * stds.';
P_norm = P ./ normMat;  % Normalized zero-lag covariance i.e. correlation

%% ————— GRID SEARCH FOR BEST RANGE —————
c_off_all = abs(P_norm);
c_off     = c_off_all(mask);
x_log     = log(dvals);
y_log     = log(c_off);

d_unique  = sort(unique(dvals));
N         = numel(d_unique);
minPts    = 100;
maxPts    = 500;
bestR2    = -Inf;
bestLR    = [NaN NaN];

fprintf('Starting grid‐search over %d unique distances (minPts=%d, maxPts=%d)\n', ...
        N, minPts, maxPts);
tic;
for iL = 1:(N-minPts)
    L = d_unique(iL);
    if mod(iL,50)==0
        fprintf('  L‐index %d/%d (L = %.4f)\n', iL, N-minPts, L);
    end
    for iR = (iL+minPts):N
        R   = d_unique(iR);
        sel = (dvals >= L) & (dvals <= R);
        nPts = nnz(sel);
        % only consider intervals with between minPts and maxPts
        if nPts < minPts || nPts > maxPts
            continue;
        end
        p    = polyfit(x_log(sel), y_log(sel), 1);
        yhat = polyval(p, x_log(sel));
        SSres= sum((y_log(sel) - yhat).^2);
        SStot= sum((y_log(sel) - mean(y_log(sel))).^2);
        R2   = 1 - SSres/SStot;
        if R2 > bestR2
            bestR2 = R2;
            bestLR = [L R];
            fprintf('   → New best R²=%.4f on [%0.4f, %0.4f] (nPts=%d)\n',...
                    bestR2, bestLR(1), bestLR(2), nPts);
        end
    end
end
toc;
fprintf('Optimal fit range: [%.4f, %.4f], R² = %.4f\n', bestLR(1), bestLR(2), bestR2);

%% ————— PLOT LOG–LOG WITH OPTIMAL RANGE —————
bestLR = [2.9073, 14.7513];
fitR = bestLR;
y_all           = P(mask);
[d_sorted, idx] = sort(dvals);
y_sorted        = y_all(idx);
xmin            = log(fitR(1));
xmax            = log(fitR(2));

figure, plot(log(d_sorted), log(abs(y_sorted))), hold on;
xline(xmin, 'r--'); xline(xmax, 'r--');
xlabel('log(distance)'); ylabel('log(|covariance|)');
title('Lag 1 (highlighted fitR interval)'); grid on;

%% ————— SLOPE vs τ USING OPTIMAL RANGE —————
% Simulation hyperparameters
max_lag  = 75;
nLags    = 200;
tau_vals = linspace(0,max_lag,nLags);
sel_all = (dvals >= bestLR(1) & dvals <= bestLR(2));

% Calculate lagged covariances (correlations)
P_tau = @(t) expm(A*t)*P;
Ctens       = nan(n,n,nLags);  
for k = 1:nLags
    t = tau_vals(k);
    Ck                = P_tau(t);
    Ctens(:,:,k)      = Ck ./ normMat;
end

a_vec = nan(1,nLags);
fprintf('Computing distance‐decay exponent a(τ) for each lag...\n');
tic;
for k = 1:nLags
    if mod(k,50)==0
        fprintf('  ... lag %d/%d\n', k, nLags);
    end
    c_k = abs(Ctens(:,:,k));
    yk  = log(c_k(mask));
    pk  = polyfit(log(dvals(sel_all)), yk(sel_all), 1);
    a_vec(k) = pk(1);
end
toc;

figure;
plot(tau_vals, a_vec);
xlabel('\tau'); ylabel('a(\tau)');
title(sprintf('Subject %d: a(\\tau) using optimal [%.2f, %.2f]', ...
      iSub, bestLR(1), bestLR(2)));
grid on;
