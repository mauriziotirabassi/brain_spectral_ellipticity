clear; clc; close all

rng(42);

% Data
dataDir    = fullfile(pwd,'data');
% distFile   = fullfile(dataDir,'mtx_euc_distance.mat'); % distance matrix
% structFile = fullfile(dataDir,'asymm_ncd_no_thr_N_74.mat'); % structural conn matrix
outDir     = fullfile(dataDir,'regressed_001_01_sim62131');
d = dir(fullfile(outDir,'*.mat')); % list of subject‐model .mat files
files = {d.name};

% Subject data
iSub = 1; % subject number
subj = load(fullfile(outDir,files{iSub}));
A = subj.A; % effective connectivity
n = size(A, 1);
I = eye(n);
Sigma_w = I * subj.output.eff_conn.NoiseVar; % noise covariance
Sigma = lyap(A, Sigma_w); % zero-lag covariance
S = 0.5 * (A * Sigma - Sigma * (A.')); % dC-Cov
hr = subj.h; % haemodynamic response

% Simulation parameters
n_time = 2e3; %1e4; % Simulation time steps
transient_length = 1e3;
tr = 0.1; % Sampling period
t = 0 : tr : (n_time + transient_length - 1) * tr;
t_sim = t(transient_length + 1:end);

%% SIMULATE
sys = ss(A, eye(n), eye(n), zeros(n));
w = mvnrnd(zeros(1, n), Sigma_w, length(t)); % Generate noise input
x = lsim(sys, w, t); % Simulate LTI response to noise
x = x(transient_length + 1:end,:); % Remove transient

% BOLD response
bold = zeros(n_time, n);
for k = 1:n
    h = hr(:,k);
    h = h / sum(h); % normalize
    bold(:,k) = conv(x(:,k), h, "same");
end

% Plot the norm of the state vector
figure, plot(t_sim, vecnorm(x, 2, 2))
xlabel('time'); ylabel('||x||_2');
title('Stochastic Input Response')

% State evolution in 2D/3D after PCA
[coeff, score, latent, tsquared, explained, mu] = pca(x);
% x_reduced3 = score(:, 1:3); animate3(x_reduced3)
% x_reduced2 = score(:, 1:2); animate2(x_reduced2)

%% SYNCHRONIZATION
% Filter to narrow band
fpass = [0.1, 1.0];
fs = 1 / tr;
[b, a] = butter(4, fpass / (fs / 2), 'bandpass'); % 4th-order Butterworth filter
bold_filt = filtfilt(b, a, bold);
z = hilbert(bold_filt);
phi = angle(z);

% figure, axis equal
% theta = linspace(0, 2*pi, 100);
% for k = 1:size(phi, 1)
%     cla
%     Rvec = mean(exp(1i*phi(k,:)));
%     plot(cos(theta), sin(theta), 'k--'), hold on
%     quiver(zeros(1,n), zeros(1,n), cos(phi(k,:)), sin(phi(k,:)), 0, 'b-')
%     hold on
%     quiver(0, 0, real(Rvec), imag(Rvec), 0, 'r', 'LineWidth', 2)
%     title(sprintf('t = %.2f, |R| = %.2f', t(k), abs(Rvec)))
%     drawnow, pause(.1)
% end

R = abs(mean(exp(1i*phi(:,:)), 2));
x = t_sim(:);
r = R(:);
win_sec = 30;
dt = x(2)-x(1);
w = max(1, round(win_sec / dt));
mu = movmean(r, w); % centered moving mean
s  = movstd (r, w, 0); % centered moving std (sample std)

figure, hold on
fill([x; flipud(x)], [mu+s; flipud(mu-s)], [0.9 0.9 1])
plot(x, r, 'b'), plot(x, mu, 'r-', 'LineWidth', 2)
yline(mean(r), 'k--', 'LineWidth', 2)                      
xlabel('Time (s)'); ylabel('R(t)'); title('Moving mean ± moving std')

%% COVARIANCE MATRIX
% Max lag has to be inferior to total simulation time (in transient)
% because otherwise there wouldn't be enough data points to form the lagged
% pair. If I choose maxLag > length(t) then I have empty plot.
maxLag = 400;  % number of lags to evaluate
lags = (0:maxLag) * tr;

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Sigma_th = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Sigma_th(:,:,k) = Sigma_tau(tau) ./ normMat_th; % Pearson correlation
end

% Empirical covariance
Sigma_emp0 = (x' * x) / size(x,1);
stds_emp = sqrt(diag(Sigma_emp0));
normMat_emp = stds_emp * stds_emp';
Sigma_emp = nan(n,n,numel(lags));
for k = 0:maxLag
    X = x(1:end-k, :);
    X_lag = x(1+k:end, :);
    Sigma_emp(:,:,k+1) = (X_lag' * X) / (size(x, 1) - k) ./ normMat_emp;
end

% Extend to negative lags based on property C_ij(-tau) = C_ji(tau)
lags_full = [-fliplr(lags(2:end)), lags]; % Symmetric lag vector
Sigma_th_full = nan(n,n,numel(lags_full));
Sigma_emp_full = nan(n,n,numel(lags_full));
for i = 1:n
    for j = 1:n
        Cth_pos = squeeze(Sigma_th(i,j,:));
        Cth_neg = squeeze(Sigma_th(j,i,2:end));
        Sigma_th_full(i,j,:) = [flipud(Cth_neg); Cth_pos];

        Cemp_pos = squeeze(Sigma_emp(i,j,:));
        Cemp_neg = squeeze(Sigma_emp(j,i,2:end));
        Sigma_emp_full(i,j,:) = [flipud(Cemp_neg); Cemp_pos];
    end
end

%% PLOT COVARIANCE
figure;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n+j);
        plot(lags_full, squeeze(Sigma_th_full(i,j,:)),'r-'); hold on;
        plot(lags_full, squeeze(Sigma_emp_full(i,j,:)),'b--');
        title(sprintf('\\Sigma_{%d%d}(\\tau)',i,j));
        xlabel('\tau'); grid on;
        if i==1 && j==1
            legend('Theory','Empirical');
        end
    end
end
sgtitle('Auto/Cross-Covariance');

figure;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n+j);
        Cth = squeeze(Sigma_th_full(i,j,:));
        Cth_sym = 0.5 * (Cth + flipud(Cth));
        Cth_asym = 0.5 * (Cth - flipud(Cth));
        Cemp = squeeze(Sigma_emp_full(i,j,:));
        Cemp_sym = 0.5 * (Cemp + flipud(Cemp));
        Cemp_asym = 0.5 * (Cemp - flipud(Cemp));

        plot(lags_full, Cth_sym, 'r-', 'LineWidth', 1.2); hold on;
        plot(lags_full, Cth_asym, 'b-', 'LineWidth', 1.2);
        plot(lags_full, Cemp_sym, 'r--');
        plot(lags_full, Cemp_asym, 'b--');
        title(sprintf('Sym/Asym C_{%d%d}(\\tau)', i, j));
        xlabel('\tau'); grid on;
        if i==1 && j==1
            legend('Sym (th)','Asym (th)','Sym (emp)','Asym (emp)');
        end
    end
end
sgtitle('Symmetric (Even) vs Asymmetric (Odd) Covariance Components');

%% FUNCTIONS
function animate3(data)
    figure; 
    h = animatedline('LineWidth', 1);
    axis([min(data(:,1)) max(data(:,1)) ...
          min(data(:,2)) max(data(:,2)) ...
          min(data(:,3)) max(data(:,3))]);
    grid on
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    view(3)
    
    for i = 1:size(data, 1)
        addpoints(h, data(i, 1), data(i, 2), data(i, 3));
        drawnow
    end
end

function animate2(data)
    figure; 
    h = animatedline('LineWidth', 1);
    axis([min(data(:,1)) max(data(:,1)) ...
          min(data(:,2)) max(data(:,2))]);
    grid on
    xlabel('PC_1'); ylabel('PC_2');
    title('2D Trajectory Animation');
    
    for i = 1:size(data, 1)
        addpoints(h, data(i, 1), data(i, 2));
        drawnow
    end
end