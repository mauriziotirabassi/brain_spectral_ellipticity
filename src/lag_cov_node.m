%% SETUP
clearvars; close all; clc

dataDir   = fullfile(pwd,'data');
distFile  = fullfile(dataDir,'mtx_euc_distance.mat');
outDir    = fullfile(dataDir,'regressed_001_01_sim62131');

% lag simulation parameters
max_lag = 50;
nLags   = 200;
tau_vals= linspace(0, max_lag, nLags);

% subject and node
iSub   = 1;             % which .mat in outDir
iNode  = 1;            % which row of P to inspect

% load Euclidean distances
tmp    = load(distFile,'mtx_euc_dis');
D      = tmp.mtx_euc_dis;
Dfull  = D + D.';       % make it symmetric
n      = size(Dfull,1); % number of nodes 
dist_node        = Dfull(iNode,:).'; % distances from iNode \to all j\neq iNode
dist_node(iNode) = []; % drop self‑distance, now dist_node is (n‑1)x1

% fitting range and indices
fitR     = [4.48, 17.18];    % TODO: range choice matters
% fitR     = [4.48, 12.18];   % Benozzo
% fitR     = [8.13 33.82];    % Deco
leftR    = max(fitR(1), min(dist_node));
sel      = dist_node >= leftR & dist_node <= fitR(2);

% load one subject
flist       = dir(fullfile(outDir,'*.mat'));
subj        = load(fullfile(outDir,flist(iSub).name));
A           = subj.A;
Q           = eye(n)*subj.output.eff_conn.NoiseVar;
Sigma       = lyap(A,Q);
S           = 0.5 * (A*Sigma - Sigma*A.');
stds        = sqrt(diag(Sigma)); % n×1 vector of standard deviations
normMat     = stds * stds.';     % n×n matrix where (i,j)=σ_i*σ_j

% lagged and decoupled dynamics
P_tau       = @(t) expm(A*t)      * Sigma;
A_sym       = -0.5 * (Q / Sigma);
A_skew      = S / Sigma;
P_tau_sym   = @(t) expm(A_sym*t)  * Sigma;
P_tau_skew  = @(t) expm(A_skew*t) * Sigma;

% compute node-to-all covariances over time
nCov      = nan(n-1, nLags);
nCov_sym  = nan(n-1, nLags);
nCov_skew = nan(n-1, nLags);

for k = 1:nLags
    t = tau_vals(k);

    % full part
    C    = P_tau(t)      ./ normMat;
    C    = max(C, C.');              % elementwise max of C_ij and C_ji
    C    = abs(C);                   % absolute value
    vec       = C(iNode,:).';        % extract row iNode as column vector
    vec(iNode)= [];                  % drop the diagonal entry
    nCov(:,k) = vec;

    % symmetric part
    Csym  = P_tau_sym(t)  ./ normMat;
    Csym  = max(Csym, Csym.');
    Csym  = abs(Csym);
    vec_sym       = Csym(iNode,:).';
    vec_sym(iNode)= [];
    nCov_sym(:,k) = vec_sym;

    % skew‐symmetric part
    Csk   = P_tau_skew(t) ./ normMat;
    Csk   = max(Csk, Csk.');
    Csk   = abs(Csk);
    vec_skew       = Csk(iNode,:).';
    vec_skew(iNode)= [];
    nCov_skew(:,k) = vec_skew;
end

%% PLOT FULL FOR SINGLE TAU
lag    = 1;    % choose
y_all  = nCov(:, lag);
[d_sorted, idx] = sort(dist_node);
y_sorted        = y_all(idx);
xmin = log(leftR);
xmax = log(fitR(2));

figure;
plot(log(d_sorted), log(y_sorted)), hold on;
xline(xmin, 'r--');
xline(xmax, 'r--');
xlabel('log(distance)');
ylabel('log(|covariance|)');
title(sprintf('Lag %d (node %d)', lag, iNode));
grid on;

%% PLOT FIT FOR FULL
% preallocate
slopes = nan(1, nLags);

figure;
for k = 1:nLags
    cvec = nCov_skew(:, k);
    x = log(dist_node);
    y = log(abs(cvec));

    % Fit in log–log space
    p = polyfit(x(sel), y(sel), 1);
    slopes(k) = p(1);

    % Prepare data for plotting
    x_data = dist_node;
    y_data = abs(cvec);

    % Sort raw data by distance
    [x_data_sorted, sortIdx] = sort(x_data);
    y_data_sorted = y_data(sortIdx);

    % Fit curve in log–log, back-transform for plotting
    x_fit_log = x(sel);
    y_fit_log = polyval(p, x_fit_log);
    x_fit = exp(x_fit_log);
    y_fit = exp(y_fit_log);

    % Sort fit points too (so the red line is smooth)
    [x_fit_sorted, fitIdx] = sort(x_fit);
    y_fit_sorted = y_fit(fitIdx);

    % Plot
    clf;
    loglog(x_data_sorted, y_data_sorted, 'b-'); hold on;
    loglog(x_fit_sorted, y_fit_sorted, 'r--');
    xlim([2, 1e2]);
    ylim([1e-6, 1]);
    grid on;
    xlabel('Distance (mm)');
    ylabel('|Covariance|');
    title(sprintf('\\tau = %.2f,  slope \\alpha = %.3f', tau_vals(k), slopes(k)));
    legend('raw data','power-law fit', 'Location', 'northeast');
    drawnow;
    hold off;
end

% finally, plot slope vs lag
figure;
plot(tau_vals, slopes);
xlabel('\tau');
ylabel('Spatial decay exponent \alpha(\tau)');
title(sprintf('Subject %d, Node %d', iSub, iNode));
grid on;

%% REGRESSION & SLOPE VS LAG
% preallocate exponent vectors
a_vec       = nan(1, nLags);
a_sym_vec   = nan(1, nLags);
a_skew_vec  = nan(1, nLags);

% regression loop
for k = 1:nLags
    % FULL dynamics
    c_off       = nCov(:, k);
    p           = polyfit(log(dist_node(sel)), log(c_off(sel)), 1);
    a_vec(k)    = p(1);

    % SYMMETRIC dynamics
    c_off_sym   = nCov_sym(:, k);
    p_sym       = polyfit(log(dist_node(sel)), log(c_off_sym(sel)), 1);
    a_sym_vec(k)= p_sym(1);

    % SKEW‐SYMMETRIC dynamics
    c_off_skew  = nCov_skew(:, k);
    p_skew      = polyfit(log(dist_node(sel)), log(c_off_skew(sel)), 1);
    a_skew_vec(k)= p_skew(1);
end

% plot all three exponents against lag
figure;
plot(tau_vals, a_vec,       'k'  ); hold on;
plot(tau_vals, a_sym_vec,   'b--');
plot(tau_vals, a_skew_vec,  'r--');
legend('full A', 'symmetric part', 'skew-symmetric part', 'Location', 'best');
xlabel('\tau');
ylabel('Distance‐decay exponent \alpha(\tau)');
title(sprintf('Node-%d distance‐decay exponent vs lag', iNode));
grid on;
