clear; clc; %close all
rng(42);

n = 2; I = eye(n);
Sigma_w = I;
A = [-.02, 0; 0, -.2];
Sigma = lyap(A, Sigma_w);
S = 0.5 * (A * Sigma - Sigma * (A.'));

% Topology
% topology = 'Chain'; S = buildS(n, topology); % showtop(S)
% topology = 'Hub'; S2 = buildS(n, topology); % showtop(S)
% S = S + 3 * S2; showtop(S)
% S = zeros(n);

% Dynamics
% Sigma = I;
% A = (-0.5 * Sigma_w + S) / Sigma;

% % Define A
% topology = 'Random';
% S = randn(n); S = 0.5 * (S - S'); % Random skew-symmetric part
% Q = randn(n); Q = Q' * Q; % Random negative definite part
% A = -Q + S;
% Sigma = lyap(A, Sigma_w);
% S = 0.5 * (A * Sigma - Sigma * (A.'));
% showtop(S);

% TIME-LAGGED AUTO/CROSS-COVARIANCE & CORRELATION FUNCTIONS
maxLag = 1000; % Number of lags to evaluate (excluding zero lag)
delta_tau = 0.1;
lags = (0:maxLag) * delta_tau;

% Theoretical covariance
Sigma_tau = @(tau) expm(A * tau) * Sigma;
Cov_th = nan(n, n, numel(lags));
for k = 1:numel(lags)
    tau = lags(k);
    Cov_th(:,:,k) = Sigma_tau(tau);
end

% Theoretical correlation
stds_th = sqrt(diag(Sigma));
normMat_th = stds_th * stds_th.';
Corr_th = Cov_th ./ normMat_th;

figure;
if n > 4, npl = 5; else, npl = n; end
for i = 1:npl
    for j = 1:npl
        subplot(npl, npl, (i - 1) * npl + j);
        plot(lags, squeeze(Corr_th(i, j, :)), 'b-');
        title(sprintf('\\Sigma_{%d%d}(\\tau)', i, j));
        xlabel('\tau'); grid on;
    end
end
sgtitle('Auto/Cross-Covariance Functions');

% CROSS-LAG COVARIANCE (CLC)
mask = find(ones(n)); % All
% mask = find(~eye(n)); % Only off-diagonal elements
mask = find(eye(n)); % Only diagonal elements
C_corr = crosslag_corr(Corr_th, mask);
C_corr(tril(true(size(C_corr)), -1)) = NaN;

figure, tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
nexttile, h = imagesc(lags, lags, C_corr);
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
set(h, 'AlphaData', ~isnan(C_corr));
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title('Pearson Correlation')

C_cos = crosslag_cos(Corr_th, mask);
C_cos(tril(true(size(C_cos)), -1)) = NaN;

nexttile, h = imagesc(lags, lags, C_cos);
axis square; clim([-1 1]); colorbar; colormap(magma);
xlabel('Lag \tau_1'); ylabel('Lag \tau_2');
set(h, 'AlphaData', ~isnan(C_cos));
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
title('Cosine Similarity')

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

%%
% Put this right after vecs = reshape(...) inside crosslag_cos or in the workspace
vecs = reshape(Corr_th, [], size(Corr_th,3));   % same reshape you use
fprintf('size(vecs) = [%d %d]\n', size(vecs,1), size(vecs,2));

% show first few columns
disp('first 4 columns of vecs (rows x cols):');
disp(vecs(:, 1:min(4,size(vecs,2))));

% check if any column is exactly proportional to the first (exact equality of ratio)
ratios = nan(size(vecs,2),1);
for k=1:size(vecs,2)
    % avoid division by zero: find index of the first nonzero entry in col1
    idx = find(vecs(:,1) ~= 0, 1);
    if isempty(idx)
        ratios(k) = NaN;
    else
        ratios(k) = vecs(idx,k)/vecs(idx,1);   % if proportional, this ratio should satisfy vecs(:,k)==ratio*vecs(:,1)
    end
end
fprintf('ratio of each column to col1 at first nonzero row: (first 8)\\n');
disp(ratios(1:min(8,end))');

% check exact proportionality boolean
isProp = false(1,size(vecs,2));
for k=1:size(vecs,2)
    r = ratios(k);
    if ~isnan(r)
        isProp(k) = all(vecs(:,k) == r * vecs(:,1));
    end
end
fprintf('any column exactly proportional to column1? %d (1=yes)\n', any(isProp));

% check pairwise cosine numerically and min/max
dotp = vecs' * vecs;
norms = sqrt(sum(vecs.^2,1));
den = norms' * norms;
Ccheck = dotp ./ den;
fprintf('min(Ccheck) = %.17g, max(Ccheck) = %.17g\n', min(Ccheck(:)), max(Ccheck(:)));

% check rank of vector set
fprintf('rank(vecs) = %d\n', rank(vecs));


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

function C = crosslag_corr(matrixvec, mask)
%CROSSLAG_CORR Construct cross-lag covariance matrix (CLC) Pearson correaltion
%   matrixvec : series of matrices whose CLC to calculate
%   mask      : mask for the matrices
    vecs = reshape(matrixvec, [], size(matrixvec,3)); % vectorize each matrix
    if nargin >= 2 && ~isempty(mask), vecs = vecs(mask(:), :); end
    C = corr(vecs);
end

function C = crosslag_cos(matrixvec, mask)
%CROSSLAG_COS cosine similarity
    vecs = reshape(matrixvec, [], size(matrixvec,3));
    if nargin >= 2 && ~isempty(mask), vecs = vecs(mask(:), :); end
    C = (vecs' * vecs) ./ (sqrt(sum(vecs.^2,1))' * sqrt(sum(vecs.^2,1)));
end
