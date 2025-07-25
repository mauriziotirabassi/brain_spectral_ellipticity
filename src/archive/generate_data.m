%% silulate new data
function y_sim = generate_data(A, n_time, pwd_data, filename, only_state)
    n = size(A, 1);
    sys = load([pwd_data, filename.name], 'sys_est').sys_est;
    em_params = load([pwd_data, filename.name], 'em_params').em_params;
    %em_params.eff_conn.Ts = load([pwd_data, filename.name], 'TR').TR; em_params.hemo.m = size(load([pwd_data, filename.name], 'h').h, 1); %LEMON
    sys.A = [expm(A*em_params.eff_conn.Ts) zeros(n,(em_params.hemo.m-1)*n); eye(n*(em_params.hemo.m-1)) zeros((em_params.hemo.m-1)*n,n)];
    if only_state
        sys.C = [eye(n) zeros(n,(em_params.hemo.m-1)*n)]; %to remove hemo contribution
    end
    TR_sys = em_params.eff_conn.Ts;
    transient_length = 1e3;
    t = 0 : TR_sys : (n_time + transient_length - 1) * TR_sys;
    w = mvnrnd(zeros(n_time + transient_length, size(sys.Q, 1)), sys.Q);
    y_sim = lsim(ss(sys.A, sys.B, sys.C, sys.D, TR_sys), w, t);
    y_sim = y_sim(transient_length + 1:end,:);
end
