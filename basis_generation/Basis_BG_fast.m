function [Q, H, R] = Basis_BG_fast(n,A,b,m,thetagenfun,lssolver,lowprecision)
%  Randomized Gram-Schmidt process
    if (nargin < 2)
        error(message('randgmres: not enough input variables'));
    end
    
    if isa(A,'sparse')||isa(A,'double')
        A = @(x) A*x;
    elseif ~isa(A,'function_handle')
        error(message('randgmres: l.h.s. must be a square matrix or function handle'));
    end
    
    
    % Assign default values to unspecified parameters.    
    if (nargin < 5) || isempty(thetagenfun)
        k = min(n, ceil(2*m*log(n)/log(m)));
        k = min(n, 2*m);
        thetagenfun = @() SRHT(n,k);
        % Generate random sketching matrix.
        Theta = thetagenfun();
    else
        Theta = thetagenfun;
        k = size(Theta(zeros(n,1)),1);
    end
    
    if (nargin < 6) || isempty(lssolver)
        lssolver = '5reorth';
        lssolver = '3reorth';
    end
    
    if (nargin < 7) || isempty(lowprecision)
        lowprecision = 0;
    end
    
    if (nargin > 7)
        error(message('randgmres: too many input variables'));
    end
    
    % Normalize r.h.s.
    n2b = norm(b);
    rhs = b/n2b;
    
    % Allocate space for Krylov basis matrix.
    Q = zeros(n,m,'double');
    
    % randomized GMRES
        
    % Allocate space for small matrices and vectors.
    R = zeros(m,m);
    S = zeros(k,m);
    P = zeros(k,m);
    
    % Perform first iteration.
    s = Theta(rhs);
    p = s;
    P(:,1) = p;
    r = [norm(s);zeros(m-1,1)];
    s = s/r(1);
    q = rhs/r(1);

    for initer=1:m-1

        Q(:,initer) = q;
        S(:,initer) = s;
        R(:,initer) = r;

        % Perform RGS iteration.
        w = A(q);
    
%         tic;
        p = Theta(w);
%         tottime = tottime + toc;
        P(:,initer+1) = p;

        r = leastsquares(S,p,lssolver,initer);
        r(initer+1:m) = zeros(m-initer,1);

        if lowprecision
            q = single(w) - Q*single(r);
            s = Theta(double(q));
        else
            q = w -Q*r;
%             tic;
            % s = Theta(q);
            s = p - S*r; % CHANGED HERE
%             tottime = tottime + toc;
        end

        r(initer+1) = norm(s);
        q = q/r(initer+1);
        s = s/r(initer+1);
    end
    H = R(:,2:end);
end

% Solving nearly orthogonal least-squares problem
function r = leastsquares(S,p,lssolver,initer)
    if strcmp(lssolver,'3reorth')
        j = 3;
    elseif strcmp(lssolver,'5reorth')
        j = 5;
    elseif strcmp(lssolver,'20reorth')
        j = 20;
    else
        j = 0;
    end
    if j ~= 0
        % Richardson iterations
        r = (p'*S)';
        p = p - S*r;
        for i=1:j-1
            dr = (p'*S)';
            p = p - S*dr;
            r = r + dr;
        end
    else
        % Conjugate Gradient
        Stemp = S(:,1:initer);
        [rtemp,~,~,~,~] = pcg(@(x) ((Stemp*x)'*Stemp)',Stemp'*p,1.0e-14,20);
        r = [rtemp; zeros(size(S,2) - initer,1)];
    end
end


% Subsampled randomized Hadamard Transform (SRHT).
function Theta = SRHT(n,k)
    D = randi([0 1], n,1)*2 - 1;
    N = 2^ceil(log(n)/log(2));
    perm = randperm(N,k);
    select = @(t,ind) t(ind);
    Theta = @(t) (1/sqrt(k)) * select(myfwht(D.*t),perm);
end


% Fast Walsh Hadamard Transform
function z = myfwht(a)
    h = 1;
    n = length(a);
    N = 2^ceil(log(n)/log(2));
    z = zeros(N,1);
    z(1:n) = a;
    while h < N
        for i = 1:2*h:N
            for j = i:(i + h - 1)
                x = z(j);
                y = z(j + h);
                z(j) = x + y;
                z(j + h) = x - y;
            end
        end
        h = 2*h;
    end
end