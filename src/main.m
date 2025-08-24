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

%% COVARIANCE DYNAMICS
Sigma0 = cov(y_stoch);  % steady-state covariance (empirical estimate)

Sigma_tau = @(tau) expm(A * tau) * Sigma0;
stds       = sqrt(diag(Sigma0));
normMat    = stds * stds.';

maxLag = 500;  % number of lags to evaluate
lags = (0:maxLag) * tr;

% Theoretical covariance
Sigma_th = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Sigma_th(:,:,k) = Sigma_tau(tau) ./ normMat; % Pearson correlation
end

% Empirical covariance
Sigma_emp = nan(n,n,numel(lags));
for k = 1:numel(lags)
    tauSteps = round(lags(k)/tr);
    X1 = y_stoch(1:end-tauSteps,:);
    X2 = y_stoch(1+tauSteps:end,:);
    Sigma_emp(:,:,k) = (X2')*X1 / (size(X1,1)-1); % E[x(t+tau)x(t)^T]
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
sgtitle('Lagged Covariance Evolution');

% for k = transient_length:(n_time + transient_length - 1)
%     imshow(Ctens(:,:,k), [], 'InitialMagnification', 'fit')
%     title(sprintf('Time-Lagged Covariance \\Sigma(\\tau) at lag \\tau = %.2f', t(k)));
%     drawnow
%     pause(0.01)
% end

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