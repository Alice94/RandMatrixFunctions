%% Running standard Arnoldi, Balabanov-Grigori, and Nakatsukasa-Tropp to 
% find a decomposition of the form A V_m = V_{m+1} H_{m}. Then if V_m does
% not have orthonormal columns we find the last column of M := V_m^+ A V_m
% using Blenenpik or LSQR for solving the least squares problem. Finally we get
% f(A)b \approx V_m f(M) e_1, with some scaling factors if necessary.

errArnoldi = [];
errNTB = [];
errNTW = [];
errNTtrid = [];
errBG = [];
errsFOM = [];
errArnoldiRestart = [];

timeArnoldi = [];
timeNTB = [];
timeNTW = [];
timeNTtrid = [];
timeBG = [];
timesFOM = [];
timeArnoldiRestart = [];

for m = ms
    fprintf("-----------------------------m = %d-----------------------\n", m);
    
    if (algorithms(1) == 1)
        % Standard Arnoldi (modified Gram-Schmidt without reorthogonalization)
        tic;
        [V, H] = Basis_Arnoldi(n, Afun, b, m);
        mm = size(H,2);
        y = V(:,1:mm) * f(H(1:mm, :), [1; zeros(mm-1, 1)]);
        timeArnoldi = [timeArnoldi, toc];
        errArnoldi = [errArnoldi, norm(y - y_ex)/norm(y_ex)];
        fprintf("Error of Arnoldi for dimension %d is %1.2e, time is %1.2f\n", mm, errArnoldi(end), timeArnoldi(end));
    end
    
    if (algorithms(2) == 1)
        % Balabanov/Grigori with LSQR for the last column
        tic;
        if NAME == "invsqrt3DLap"
            [V, H, R] = Basis_BG_fast(n, Afun, b, round(m*1.05));
        else
            [V, H, R] = Basis_BG_fast(n, Afun, b, m*2);
        end
        x = lsqr(V(:,1:m), H(m+1,m)*V(:,m+1), 1e-6, 1000);
        H = H(1:m, 1:m); H(:,end) = H(:,end) + x;
        y = V(:,1:m) * f(H, [R(1,1); zeros(m-1, 1)]);
        timeBG = [timeBG, toc];
        errBG = [errBG, norm(y - y_ex)/norm(y_ex)];
        fprintf("Error of Balabanov-Grigori LSQR (FWHT) for dimension %d is %1.2e, time is %1.2f\n", m, errBG(end),...
            timeBG(end));
    end
   
    if (algorithms(3) == 1)
        % Nakatsukasa/Tropp + Blendenpik
        tic;
        [V, H] = Basis_NakatsukasaTropp(n, Afun, b, m);
        x = blendenpik(V(:,1:end-1), V(:,end)) * H(end,end);
        M = H(1:m, :); M(:,end) = M(:,end) + x;
        y = V(:,1:m) * f(M, [1; zeros(m-1, 1)]);
        timeNTB = [timeNTB, toc];
        errNTB = [errNTB, norm(y - y_ex)/norm(y_ex)];
        fprintf("Error of Nakatsukasa-Tropp for dimension %d is %1.2e, time is %1.2f\n", m, errNTB(end), timeNTB(end));
        fprintf("norm(A*V(:,1:m) - V*H) = %1.2e\n", norm(A*V(:,1:m) - V*H)) 
        fprintf("******* Condition number of V is %1.2e\n", cond(V))
    end
    
    if (algorithms(4) == 1)
        % Nakatsukasa/Tropp without least squares
        tic;
        [V, H] = Basis_NakatsukasaTropp(n, Afun, b, m);
        M = H(1:m, :);
        y = V(:,1:m) * f(M, [1; zeros(m-1, 1)]);
        timeNTtrid = [timeNTtrid, toc];
        errNTtrid = [errNTtrid, norm(y - y_ex)/norm(y_ex)];
        fprintf("Error of NT without least-squares for dimension %d is %1.2e, time is %1.2f\n", m, errNTtrid(end), timeNTtrid(end));
        fprintf("norm(A*V(:,1:m) - V*H) = %1.2e\n", norm(A*V(:,1:m) - V*H)) 
        fprintf("******* Condition number of V is %1.2e\n", cond(V))
    end

    if (algorithms(5)==1)
        % sFOM from Guettel/Schweitzer paper
        tic;
        k = 2; % how many blocks to orthogonalise against
        hS = setup_sketching_handle(n,min(n, 2*m));
        [SV,SAV,Vtrunc] = bta_mod(A,b,m,k,hS); % number of mat-vec products = m
        % whitening the basis
        [SV, SAV, Rw] = whiten_basis(SV, SAV);
        % Compute errors
        errsFOM = [errsFOM, sfom_closed_eval_error_mod(Vtrunc,SV,...
            SAV,hS(b),Rw,y_ex,f,m)];
        timesFOM = [timesFOM, toc];
        fprintf("Error of sFOM for dimension %d is %1.2e, time is %1.2f\n", m, errsFOM(end), timesFOM(end));
    end

    if (algorithms(6) == 1)
        % Nakatsukasa/Tropp + whitening + Blendenpik
        tic;
        [V, H, whitened] = Basis_NakatsukasaTroppWhite(n, Afun, b, m);
        if (whitened)
            x = lsqr(V(:,1:end-1), H(end,end)*V(:,end), 1e-6, 1000);
        else
            x = blendenpik(V(:,1:end-1), V(:,end)) * H(end,end);
        end
        M = H(1:m, :); M(:,end) = M(:,end) + x;
        y = V(:,1:m) * f(M, [1/norm(V(:,1)); zeros(m-1, 1)]);
        timeNTW = [timeNTW, toc];
        errNTW = [errNTW, norm(y - y_ex)/norm(y_ex)];
        fprintf("Error of Nakatsukasa-Tropp-whitening for dimension %d is %1.2e, time is %1.2f\n", m, errNTW(end), timeNTW(end));
        fprintf("norm(A*V(:,1:m) - V*H) = %1.2e\n", norm(A*V(:,1:m) - V*H)) 
        fprintf("******* Condition number of V is %1.2e\n", cond(V))
    end

    if (algorithms(7) == 1)
        % Restarted Arnoldi (every 20 iterations)
        tic;
        param.function = ff; 
        param.max_restarts = round(m/param.restart_length);                % perform at most 20 restart cycles
        if (mod(m, param.restart_length) ~= 0)
            timeArnoldiRestart = [timeArnoldiRestart, NaN];
            errArnoldiRestart = [errArnoldiRestart, NaN];
            continue;
        end
        param.transformation_parameter = 1;     % parameter for the integral transformation
        param.hermitian = 0;                    % set 0 if A is not Hermitian
        param.V_full = 0;                       % set 1 if you need Krylov basis
        param.H_full = 0;                       % do not store all Hessenberg matrices
        param.exact = [];                       % Exact solution. If not known set to []
        param.stopping_accuracy = 0    ;        % stopping accuracy
        param.inner_product = @(a,b) b'*a;      % use standard euclidean inner product
        param.thick = [];                       % no implicit deflation is performed
        param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95 
        param.waitbar = 0;                      % show waitbar 
        param.reorth_number = 0;                % #reorthogonalizations
        param.truncation_length = inf;          % truncation length for Arnoldi
        param.verbose = 1;                      % print information about progress of algorithm
        [y,out1] = funm_quad(A, b, param);
%         [V, H] = Basis_ArnoldiRestart(n, Afun, b, m, restart_dim);
%         mm = size(H,2);
%         y = V(:,1:mm) * f(H(1:mm, :), [1; zeros(mm-1, 1)]);
        timeArnoldiRestart = [timeArnoldiRestart, toc];
        errArnoldiRestart = [errArnoldiRestart, norm(y - y_ex)/norm(y_ex)];
        fprintf("Error of restarted Arnoldi for dimension %d is %1.2e, time is %1.2f\n", ...
            mm, errArnoldiRestart(end), timeArnoldiRestart(end));
    end

end

%%




save(NAME, 'algorithms', 'ms', 'timeArnoldi', 'errArnoldi', ...
    'timeBG', 'errBG', 'timeNTB', 'errNTB', 'timeNTtrid', ...
    'errNTtrid', 'timeNTW', 'errNTW', 'timesFOM', 'errsFOM', ...
    'errArnoldiRestart', 'timeArnoldiRestart');
