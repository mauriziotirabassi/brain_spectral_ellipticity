clear; clc; %close all

rng(42); % This affects simulated noise injections
n = 10; I = eye(n);

figure

loop = logspace(log10(0.001), log10(100), 50);

for i = 1:length(loop)

    s = loop(i);

    Sigma_w = I; % Uncorrelated noise
    
    % Topology
    topology = 'Ring';
    S = s * buildS(n, topology);
    % showtop(S)
    
    % S = zeros(n);
    
    % Initial energy distribution
    Sigma = I; % Balanced
    % eps = 1e-1; Sigma = eps * I; Sigma(5,5) = 10; % Unbalanced
    
    % Dynamics
    A = (-0.5 * Sigma_w + S) / Sigma;
    
    % % Define A
    % topology = 'Random';
    % S = randn(n); S = 0.5 * (S - S'); % Random skew-symmetric part
    % Q = randn(n); Q = Q' * Q; % Random negative definite part
    % A = -Q + S;
    % Sigma = lyap(A, Sigma_w);
    % S = 0.5 * (A * Sigma - Sigma * (A.'));
    % showtop(S);
    
    % ev = eig(A); fprintf('max real part = %.4g\n', max(real(ev)));
    
    % COVARIANCE MATRIX
    maxLag = 5000 / 2; % Number of lags to evaluate
    delta_tau = 0.1;
    lags = (0:maxLag) * delta_tau;
    
    % Theoretical covariance
    Sigma_tau = @(tau) expm(A * tau) * Sigma;
    Cov_th = nan(n,n,numel(lags));
    for k = 1:numel(lags)
        tau = lags(k);
        Cov_th(:,:,k) = Sigma_tau(tau);
    end

    % Theoretical correlation
    stds_th = sqrt(diag(Sigma));
    normMat_th = stds_th * stds_th.';
    Corr_th = Cov_th ./ normMat_th; 

    % Extend to negative lags based on property C_ij(-tau) = C_ji(tau)
    lags_full = [-fliplr(lags(2:end)), lags]; % Symmetric lag vector
    Cov_th_neg = flip(Cov_th(:,:,2:end), 3);
    Cov_th_neg = permute(Cov_th_neg, [2 1 3]);
    Cov_th_full = cat(3, Cov_th_neg, Cov_th);

    Corr_th_neg = flip(Corr_th(:,:,2:end), 3);
    Corr_th_neg = permute(Corr_th_neg, [2 1 3]);
    Corr_th_full = cat(3, Corr_th_neg, Corr_th);
    
    % CROSS-LAG COVARIANCE (CLC) <--------------------------------------
    % mask = find(triu(ones(n),0)); % Upper triangular for negative
    mask = find(ones(n));
    clc = crosslagcov(Corr_th, mask);

    clc_ut = triu(clc, 1);
    clc_ut(clc_ut == 0) = NaN;
    h = imagesc(lags, lags, clc_ut);
    axis square; clim([-1 1]); colorbar; colormap(magma);
    xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
    set(h, 'AlphaData', ~isnan(clc_ut));
    set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
    % title(sprintf('Theoretical %s', topology));

    drawnow, pause(0.0001)
end

%% AUTO/CROSS-COVARIANCE FUNCTIONS
figure;
n = 4;
for i = 1:n
    for j = 1:n
        subplot(n,n,(i-1)*n+j);
        plot(lags_full, squeeze(Corr_th_full(i,j,:)),'b-');
        title(sprintf('\\Sigma_{%d%d}(\\tau)',i,j));
        xlabel('\tau'); grid on;
    end
end
sgtitle('Auto/Cross-Covariance');

%% SINGLE ROWS
% figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
% clc = nan(length(lags_full),length(lags_full),numel(n));
% for i = 1:n
%     mask = false(size(A)); mask(:, i) = true;
%     clc(:,:,i) = crosslagcov(Sigma_th_full, mask);
%     nexttile(1), imagesc(mask);
%     nexttile(2), imagesc(lags_full, lags_full, clc(:,:,i))
%     colorbar, colormap(magma)
%     drawnow, pause(0.2)
% end
% mask = find(triu(ones(length(lags_full)), 0));
% figure, imagesc(crosslagcov(clc, mask)), colorbar, colormap(magma)

%% FUNCTIONS
function S = buildS(n, topology)
%BUILD S Construct skew-symmetric adjacency matrix S
%   n        : number of nodes
%   topology : 'chain', 'uni_chain', 'ring', 'star'

S = zeros(n);
switch lower(topology)
    case 'chain'  % chain
        for i = 1:n-1
            S(i,i+1) = 1;
            S(i+1,i) = -1;
        end
    case 'ring'
        for i = 1:n
            j = mod(i,n) + 1;
            S(i,j) = 1;
            S(j,i) = -1;
        end
    case 'hub'
        for j = 2:n
            S(1,j) = 1;
            S(j,1) = -1;
        end
    otherwise
        error('Unknown topology: %s', topology)
end
end

function showtop(S)
    G = digraph(S);
    figure; h = plot(G, 'Layout','circle', 'EdgeLabel',G.Edges.Weight);
    h.LineWidth = abs(G.Edges.Weight)/max(abs(G.Edges.Weight));
end

function clc = crosslagcov(matrixvec, mask)
%CLC Construct cross-lag covariance matrix (CLC)
%   matrixvec : series of matrices whose CLC to calculate
%   mask      : mask for the matrices
    vecs = reshape(matrixvec, [], size(matrixvec,3)); % vectorize each matrix
    vecs = vecs(mask(:), :); % apply mask
    clc = corr(vecs);
end