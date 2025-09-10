clear; clc; close all

rng(42);
n = 3;
I = eye(n);
Sigma_w = 0.1 * I; % Uncorrelated noise

% Derive A
S = randn(n);
S = 0.5 * (S - S'); % Ensure skew-symmetry
Sigma = [7 0 3; 0 1 0; 3 0 3];
% Sigma = I;
A = (-1/2 * Sigma_w + S) / Sigma;

% Derive Sigma
% S = randn(n); S = 0.5*(S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q'*Q; % Random negative definite part
% A = -Q + S;            
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * A.');

ev = eig(A);
fprintf('max real part = %.4g\n', max(real(ev)));

% Simulation parameters
n_time = 2e3; %1e4; % Simulation time steps
transient_length = 1e3;
tr = 0.1; % Sampling period
t = 0 : tr : (n_time + transient_length - 1) * tr;

%% SIMULATE W/STOCHASTIC INPUT
sys = ss(A, eye(n), eye(n), zeros(n));
w = mvnrnd(zeros(1, n), Sigma_w, length(t)); % Generate noise input
y_stoch = lsim(sys, w, t); % Simulate LTI response to noise
y_stoch = y_stoch(transient_length + 1:end,:); % Remove transient

% Plot the norm of the state vector
figure, plot(t(transient_length + 1:end), vecnorm(y_stoch, 2, 2))
xlabel('time'); ylabel('||x||_2');
title('Stochastic Input Response')

% Animate state vector evolution
% animate3(y_stoch)

% State evolution in 2D/3D after PCA
% [coeff, score, latent, tsquared, explained, mu] = pca(y_sim);
% y_reduced3 = score(:, 1:3); animate3(y_reduced3)
% y_reduced2 = score(:, 1:2); animate2(y_reduced2)

%% STATE EVOLUTION & STEADY-STATE COVARIANCE
[U, D] = eig(Sigma); % eigenvectors & eigenvalues
radii = sqrt(diag(D)); % principal axes
mean_x = mean(y_stoch,1)'; % mean of the samples

% Plot 1std ellipsoid
[xs, ys, zs] = sphere(20); % unit sphere
ellipsoid_pts = U * diag(radii) * [xs(:)'; ys(:)'; zs(:)'];
X = reshape(ellipsoid_pts(1,:) , size(xs));
Y = reshape(ellipsoid_pts(2,:) , size(ys));
Z = reshape(ellipsoid_pts(3,:) , size(zs));

figure; hold on; grid on; view(3);
surf(X, Y, Z, 'FaceAlpha',0.2, 'EdgeColor','none');
plot3(mean_x(1), mean_x(2), mean_x(3), 'k+', 'MarkerSize',10);
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');

% Plot eigenvectors scaled by eigenvalues
for i = 1:3
    vec = radii(i) * U(:,i);   % scale eigenvector by sqrt of eigenvalue
    quiver3(mean_x(1), mean_x(2), mean_x(3), vec(1), vec(2), vec(3), ...
        'LineWidth', 1);
end

% Animate trajectory
h = animatedline('LineWidth',1);
for i = 1:n_time
    addpoints(h, y_stoch(i,1), y_stoch(i,2), y_stoch(i,3));
    drawnow;
end

%% SYNCHRONIZATION
% Filter to narrow band
fpass = [0.1, 1.0];
fs = 1 / tr;
[b, a] = butter(4, fpass / (fs / 2), 'bandpass'); % 4th-order Butterworth filter
y_filt = filtfilt(b, a, y_stoch);
z = hilbert(y_filt);
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

t_steady = t(transient_length + 1:end);
R = abs(mean(exp(1i*phi(:,:)), 2));
x = t_steady(:);
r = R(:);
win_sec = 30;
dt = x(2)-x(1);
w = max(1, round(win_sec / dt));
mu = movmean(r, w);       % centered moving mean
s  = movstd (r, w, 0);    % centered moving std (sample std)

figure, hold on
fill([x; flipud(x)], [mu+s; flipud(mu-s)], [0.9 0.9 1])
plot(x, r, 'b'), plot(x, mu, 'r-', 'LineWidth', 2)
yline(mean(r), 'k--', 'LineWidth', 2)                      
xlabel('Time (s)'); ylabel('R(t)'); title('Moving mean ± moving std')

%% SIMULATE AUTONOMOUS
x0 = 2 * ones(n,1);  % initial condition
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t_aut, y_aut] = ode45(@(t,x) A*x, t, x0, options);

% Plot the norm of the state vector
figure, plot(t_aut, vecnorm(y_aut, 2, 2))
xlabel('time'); ylabel('||x||_2');
title('Autonomous response (ode45)')

% Animate state vector evolution
animate3(y_aut)

%% COVARIANCE MATRIX
% Max lag has to be inferior to total simulation time (in transient)
% because otherwise there wouldn't be enough data points to form the lagged
% pair. If I choose maxLag > length(t) then I have empty plot. It is not
% useful to consider a number of lags close to the simulated time points
% because estimate variance will be too high, given the fewer number of
% points
maxLag = min(size(y_stoch, 1)/2, 2000);  % number of lags to evaluate
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
Sigma_emp0 = (y_stoch' * y_stoch) / size(y_stoch,1);
stds_emp = sqrt(diag(Sigma_emp0));
normMat_emp = stds_emp * stds_emp';
Sigma_emp = nan(n,n,numel(lags));
for k = 0:maxLag
    X = y_stoch(1:end-k, :);
    X_lag = y_stoch(1+k:end, :);
    Sigma_emp(:,:,k+1) = (X_lag' * X) / (size(y_stoch, 1) - k) ./ normMat_emp;
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

%% PLOT COVARIANCES
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

% for k = transient_length:(n_time + transient_length - 1)
%     imshow(Ctens(:,:,k), [], 'InitialMagnification', 'fit')
%     title(sprintf('Time-Lagged Covariance \\Sigma(\\tau) at lag \\tau = %.2f', t(k)));
%     drawnow
%     pause(0.01)
% end

%% FCD-STYLE SIMILATIRY ACROSS LAGS
% Seeing at what lag the dynamics are the same.
triuIdx = find(triu(ones(n),1));

% Theoretical
vecs_th = [];
for k = 1:length(lags_full)
    Ck = Sigma_th_full(:,:,k);
    vecs_th(:,k) = Ck(triuIdx);
end
FCD_th = corr(vecs_th);

% Empirical
vecs_em = [];
for k = 1:length(lags_full)
    Ck = Sigma_emp_full(:,:,k);
    vecs_em(:,k) = Ck(triuIdx);
end
FCD_emp = corr(vecs_em);

% Plot
figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, imagesc(lags_full, lags_full, FCD_th);
axis square; colorbar; colormap jet;
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title('Theoretical Time-Lagged Covariance Dynamics');

nexttile, imagesc(lags_full, lags_full, FCD_emp);
axis square; colorbar; colormap jet;
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
title('Empirical Time-Lagged Covariance Dynamics');

%% PSD MATRIX
Fs = 1 / tr;
% FFT algorithm faster with power of 2
Nfft = 2^nextpow2(maxLag + 1);
freqs = Fs * (0:(Nfft/2)) / Nfft;
omega = 2 * pi * freqs;

% Theoretical PSD
Phi_th = nan(n,n,length(freqs));
for k = 1:length(freqs)
    H = 1i * omega(k) * eye(n) - A;
    H_inv = inv(H);
    Phi_th(:,:,k) = H_inv * Sigma_w * (H_inv)';
end

% Empirical PSD
window = hamming(256);
noverlap = 128;
Phi_emp = nan(n, n, Nfft/2+1);
for i = 1:n
    for j = 1:n
        Pij = cpsd(y_stoch(:,i), y_stoch(:,j), window, noverlap, Nfft, Fs);
        Phi_emp(i,j,:) = Pij;  % Complex cross-spectrum
    end
end

figure;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n+j);
        plot(omega, abs(squeeze(Phi_th(i,j,:))),'r-'); hold on;
        plot(omega, abs(squeeze(Phi_emp(i,j,:))),'b--');
        title(sprintf('\\Phi_{%d%d}(\\omega)',i,j));
        xlabel('\omega'); grid on;
        if i==1 && j==1
            legend('Theory','Empirical');
        end
    end
end
sgtitle('Cross-Spectral Density Magnitude');

figure;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n + j);
        plot(omega, unwrap(angle(squeeze(Phi_th(i,j,:)))),'r-'); hold on;
        plot(omega, unwrap(angle(squeeze(Phi_emp(i,j,:)))),'b--');
        title(sprintf('Phase of \\Phi_{%d%d}(\\omega)', i, j));
        xlabel('\omega'); grid on;
        if i==1 && j==1
            legend('Theory','Empirical');
        end
    end
end
sgtitle('Cross-Spectral Density Phase');

figure;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n+j);
        co_th = real(squeeze(Phi_th(i,j,:)));
        quad_th = imag(squeeze(Phi_th(i,j,:)));
        co_emp = real(squeeze(Phi_emp(i,j,:)));
        quad_emp = imag(squeeze(Phi_emp(i,j,:)));

        plot(omega, co_th, 'r-', 'LineWidth', 1.2); hold on;
        plot(omega, quad_th, 'b-', 'LineWidth', 1.2);
        plot(omega, co_emp, 'r--', 'LineWidth', 1);
        plot(omega, quad_emp, 'b--', 'LineWidth', 1);
        title(sprintf('Co/Quad Spectrum \\Phi_{%d%d}(\\omega)', i, j));
        xlabel('\omega'); grid on;
        if i==1 && j==1
            legend('Co-spectrum (theory)','Quad (theory)','Co-spectrum (emp)','Quad (emp)');
        end
    end
end
sgtitle('Co- and Quadrature Spectra');

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