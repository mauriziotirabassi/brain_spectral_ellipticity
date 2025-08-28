clear; clc; close all

rng(42);
n = 3;
I = eye(n);

S = [ 0  -.1  0 ;
      .1   0 -1;
      0  1  0 ];
% S = randn(n);
S = 0.5 * (S - S'); % Ensure skew-symmetry
Sigma_w = 0.2 * I;
Sigma = I;
A = (-1/2 * Sigma_w + S) / Sigma;

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
title('Autonomous response (ode45)')

% Animate state vector evolution
animate3(y_stoch)

% State evolution in 2D/3D after PCA
% [coeff, score, latent, tsquared, explained, mu] = pca(y_sim);
% y_reduced3 = score(:, 1:3); animate3(y_reduced3)
% y_reduced2 = score(:, 1:2); animate2(y_reduced2)

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
Sigma_emp0 = (y_stoch' * y_stoch) / size(y_stoch,1);
stds_emp = sqrt(diag(Sigma_emp0));
normMat_emp = stds_emp * stds_emp';
Sigma_emp = nan(n,n,numel(lags));
for k = 0:maxLag
    X = y_stoch(1:end-k, :);
    X_lag = y_stoch(1+k:end, :);
    Sigma_emp(:,:,k+1) = (X_lag' * X) / (size(y_stoch, 1) - k) ./ normMat_emp;
end

figure;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n+j);
        plot(lags, squeeze(Sigma_th(i,j,:)),'r-'); hold on;
        plot(lags, squeeze(Sigma_emp(i,j,:)),'b--');
        title(sprintf('\\Sigma_{%d%d}(\\tau)',i,j));
        xlabel('\tau'); grid on;
        if i==1 && j==1
            legend('Theory','Empirical');
        end
    end
end
sgtitle('Auto/Cross-Covariance');

% for k = transient_length:(n_time + transient_length - 1)
%     imshow(Ctens(:,:,k), [], 'InitialMagnification', 'fit')
%     title(sprintf('Time-Lagged Covariance \\Sigma(\\tau) at lag \\tau = %.2f', t(k)));
%     drawnow
%     pause(0.01)
% end

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
y_freq = fft(y_stoch, Nfft);
% Extract only positive frequencies
halfIdx = 1:(Nfft/2 + 1);
y_freq = y_freq(halfIdx, :);

Phi_emp = nan(n,n,length(halfIdx));
for k = 1:length(halfIdx)
    Yf = y_freq(k, :).';
    Phi_emp(:,:,k) = (Yf * Yf') / size(y_stoch, 1);
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