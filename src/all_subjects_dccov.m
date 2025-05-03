%% ------------------------------------------------------------------------
%  SETTINGS & DATA‐FOLDER LAYOUT
%  data/
%    ├─ mtx_euc_distance.mat       % contains mtx_euc_dis
%    └─ regressed_001_01_sim62131/ % one .mat per subject: A, output.eff_conn.NoiseVar, etc.
%        ├─ sub001.mat
%        ├─ sub002.mat
%        └─ ... (20 total)
%  ------------------------------------------------------------------------
clearvars; close all; clc

dataDir  = fullfile(pwd,'data');
distFile = fullfile(dataDir,'mtx_euc_distance.mat');
outDir   = fullfile(dataDir,'regressed_001_01_sim62131');

% Subject files
d      = dir(fullfile(outDir,'*.mat'));
files  = {d.name};
nSubs  = numel(files);

%% PRELOAD CONSTANTS
tmp     = load(distFile,'mtx_euc_dis');
Dfull   = tmp.mtx_euc_dis + tmp.mtx_euc_dis.';  % symmetric
n       = size(Dfull,1);
mask    = ~eye(n);
dvals   = Dfull(mask);

% Fit‐range and fixed selection mask
fitR    = [4.48, 12.18];
sel     = dvals >= fitR(1) & dvals <= fitR(2);

% Lags
max_lag = 75;
nLags   = 200;
tau_vals= linspace(0, max_lag, nLags);

% Preallocate
a_full = nan(nSubs, nLags);
a_sym  = nan(nSubs, nLags);
a_skew = nan(nSubs, nLags);

%% LOOP OVER SUBJECTS
for iSub = 1:nSubs
    S   = load(fullfile(outDir, files{iSub}));
    A   = S.A;
    Q   = eye(n) * S.output.eff_conn.NoiseVar;
    P   = lyap(A, Q);

    % Generators
    A_sym  = -0.5 * (Q / P);
    Ssk    = 0.5*(A*P - P*A.');
    A_skew = Ssk / P;

    % Lagged‐covariance functions
    P_full = @(t)    expm(A*t)       * P;
    P_sym  = @(t)    expm(A_sym*t)   * P;
    P_skew = @(t)    expm(A_skew*t)  * P;

    % Regression at each \tau
    for k = 1:nLags
        t = tau_vals(k);

        % full
        cf      = P_full(t);  vf = cf(mask);
        p       = polyfit(log(dvals(sel)), log(abs(vf(sel))), 1);
        a_full(iSub,k) = p(1);

        % sym
        cs      = P_sym(t);   vs = cs(mask);
        p       = polyfit(log(dvals(sel)), log(abs(vs(sel))), 1);
        a_sym(iSub,k)  = p(1);

        % skew
        ck      = P_skew(t);  vk = ck(mask);
        p       = polyfit(log(dvals(sel)), log(abs(vk(sel))), 1);
        a_skew(iSub,k) = p(1);
    end
end

%% PLOT
nCols = ceil(sqrt(nSubs));
nRows = ceil(nSubs / nCols);
figure('Position',[100 100 1200 800]);
t = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

for iSub = 1:nSubs
    ax = nexttile;
    plot(ax, tau_vals, a_full(iSub,:),  'k-'); hold(ax,'on');
    plot(ax, tau_vals, a_sym(iSub,:),   'b--');
    plot(ax, tau_vals, a_skew(iSub,:),  'r-.');
    xlim(ax, [tau_vals(1), tau_vals(end)]);
    ylim(ax, [min([a_full(:);a_sym(:);a_skew(:)]), max([a_full(:);a_sym(:);a_skew(:)])]);
    xlabel(ax,'\tau'); ylabel(ax,'\alpha');
    title(ax,sprintf('Subj %d',iSub));
    grid(ax,'on');
end

sgtitle(t,'Distance‐Decay Exponent \alpha(\tau) Across All Subjects');