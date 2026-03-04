function results = get_kappa_spectrum(A, Sw)
% GET_KAPPA_SPECTRUM  Computes spectral properties via Reordered Schur.
%
%   METHOD: Iterative Reordered Schur decomposition.
%           1. Whitens the system using the Lyapunov solution.
%           2. Isolates every mode (Real and Oscillatory) to the top-left block.
%           3. Calculates geometric invariants and energy weights.
%
%   INPUTS:
%       A  : (NxN) System interaction matrix (must be Hurwitz).
%       Sw : (NxN) Noise covariance matrix.
%
%   OUTPUTS:
%       results : Struct containing two fields:
%           .Oscillatory : Struct for complex conjugate pairs
%               .kappas       : (M_osc x 1) Intrinsic non-normality factor.
%               .omegas       : (M_osc x 1) Intrinsic solenoidal drive.
%               .decay_rates  : (M_osc x 1) Mean isotropic decay.
%               .frequencies  : (M_osc x 1) Rotation frequency.
%               .eigenvalues  : (M_osc x 1) Complex eigenvalue (one per pair).
%               .node_weights : (N x M_osc) Node participation matrix.
%           .Real : Struct for purely decaying modes
%               .eigenvalues  : (M_real x 1) Real eigenvalues.
%               .node_weights : (N x M_real) Node participation matrix.

    %% 1. Geometry Setup (Whitening)
    Sigma = lyap(A, Sw);
    Sigma_half = sqrtm(Sigma);
    Sigma_half_inv = inv(Sigma_half);
    
    % Whitened drift matrix
    A_tilde = Sigma_half_inv * A * Sigma_half;

    %% 2. Modal Analysis (Initial Decomposition)
    [U_init, T_init] = schur(A_tilde, 'real');
    N = size(A, 1);
    
    % Identify blocks on the diagonal (Oscillatory pairs vs Real singletons)
    block_indices = {}; 
    skip_next = false;
    
    for i = 1:N
        if skip_next
            skip_next = false;
            continue;
        end
        
        % Check subdiagonal for 2x2 block (Strict Non-Zero)
        if i < N && T_init(i+1, i) ~= 0
            block_indices{end+1} = [i, i+1]; %#ok<AGROW>
            skip_next = true;
        else
            block_indices{end+1} = [i]; %#ok<AGROW>
        end
    end
    
    %% 3. Initialize Output Arrays
    osc.kappas       = [];
    osc.omegas       = [];
    osc.decay_rates  = [];
    osc.frequencies  = [];
    osc.deltas       = [];
    osc.eigenvalues  = [];
    osc.node_weights = []; % Dimensions: [N x M_osc]
    
    real_modes.eigenvalues  = [];
    real_modes.node_weights = []; % Dimensions: [N x M_real]
    
    %% 4. Iterative Isolation Loop
    for k = 1:length(block_indices)
        idx_range = block_indices{k};
        
        % A. Isolation
        % Reorder block to top-left (1,1) for independent weight calculation
        select = false(N, 1);
        select(idx_range) = true;     
        
        % Rotate system to place mode at (1,1) or (1:2, 1:2)
        [U_iso, T_iso] = ordschur(U_init, T_init, select);
        
        % B. Extraction
        if length(idx_range) == 2
            % --- OSCILLATORY MODE (2x2) ---
            B_iso = T_iso(1:2, 1:2);
            props = calculate_kappa_properties(B_iso);
            eigs  = eig(B_iso);
            
            osc.kappas(end+1, 1)      = props.kappa;
            osc.omegas(end+1, 1)      = props.omega;
            osc.decay_rates(end+1, 1) = props.mu;
            osc.deltas(end+1, 1)      = props.delta;
            osc.frequencies(end+1, 1) = props.frequency;
            osc.eigenvalues(end+1, 1) = eigs(1); % Store one of the pair
            
            % Weights (Sum of 2 vectors)
            z_basis = U_iso(:, 1:2);
            x_basis = Sigma_half * z_basis;
            w = sum(x_basis.^2, 2); 
            osc.node_weights = [osc.node_weights, w];
            
        else
            % --- REAL MODE (1x1) ---
            lambda = T_iso(1,1);
            
            real_modes.eigenvalues(end+1, 1) = lambda;
            
            % Weights (Single vector)
            z_basis = U_iso(:, 1);
            x_basis = Sigma_half * z_basis;
            w = x_basis.^2; 
            real_modes.node_weights = [real_modes.node_weights, w];
        end
    end
    
    %% 5. Pack Results
    results.Oscillatory = osc;
    results.Real = real_modes;
end

function props = calculate_kappa_properties(B)
% CALCULATE_KAPPA_PROPERTIES Extracts invariants from a 2x2 matrix block.
%
%   Block B form implicitly: [ mu + h            gamma - omega ]
%                            [ gamma + omega     mu - h        ]
    
    % 1. Decomposition
    mu    = (B(1,1) + B(2,2)) / 2;
    h     = (B(1,1) - B(2,2)) / 2;
    
    off_sum  = B(2,1) + B(1,2);
    off_diff = B(2,1) - B(1,2);
    
    gamma = off_sum / 2;
    omega = abs(off_diff) / 2;
    
    % 2. Invariants
    delta_sq = h^2 + gamma^2;
    delta = sqrt(delta_sq);
    omega_sq = omega^2;
    
    % 3. Kappa Calculation
    num = omega_sq + delta_sq;
    den = abs(omega_sq - delta_sq);
    
    if den < 1e-9 * num
        kappa = Inf; 
    else
        kappa = num / den;
    end
    
    % 4. Output
    props.kappa     = kappa;
    props.omega     = omega;
    props.mu        = mu;
    props.delta     = delta;
    props.frequency = sqrt(den);
end