function [V, H, whitened] = Basis_NakatsukasaTroppWhite(n, Afun, b, m, k)
if (~exist('k', 'var'))
    k = 2;
end
whitened = false;
tol = 1e-3; % tolerance after which we pass to Balabanov/Grigori algorithm
nit = 2;
sketchsize = min(n, 2*m);
scaling = norm(b);
V = zeros(n, m+1);
V(:,1) = b / scaling;
H = zeros(m+1, m);

% Initialize sketching stuff and sketched basis
wellCond = -1;
R = zeros(m+1,m+1);
invR = zeros(m+1,m+1);
SV = zeros(sketchsize, m+1);
rademacher = sign(randn(n,1));
sample = randsample(n, sketchsize);
p = fwht(rademacher .* V(:,1));
p = p(sample);
R(1,1) = norm(p);
invR(1,1) = 1/R(1,1); % initialize inverse of R
SV(:,1) = p/R(1,1);

for i = 1:m
    q = V(:,max(i-k+1,1):i); % local orthogonalization
    w = Afun(q(:,end));
    alpha = w' * q(:,end);
    if (i > 1)
        dim = size(q, 2);
        beta = w' * q(:,1:max(1, dim-k+1));
        w = w - alpha*q(:,end) - beta*q(:,1:max(1, dim-k+1));
        H(max(1, i-k+1):i-1,i) = beta';
    else
        w = w - alpha*q(:,end);
    end
    gamma = norm(w);
    H(i,i) = alpha;
    H(i+1, i) = gamma;
    
    % update basis 
    w = w/gamma;
    V(:,i+1) = w;   

    % Update sketched basis
    p = fwht(rademacher .* w);
    p = p(sample);
    x = zeros(m+1,1);
    for j = 1:nit
        % Richardson iteration
        x = x + SV'*(p - SV*x);
    end
    R(:,i+1) = x;
    p = p - SV*x;
    alpha = norm(p);
    R(i+1,i+1) = alpha;
    SV(:,i+1) = p/alpha;
    invR(:,i+1) = 1/alpha * invR * x;
    invR(i+1,i+1) = 1/alpha;

    % check whether we are exceeding the threshold
    if (max(abs(invR(i+1,:))) > 1/tol)
        warning("Ill-conditioning detected at iteration %d\n", i);
        wellCond = i;
        break;
    end

end

% If ill-conditioning was detected, change the basis and from now on
% proceed as in Balabanov/Grigori
if (wellCond > 0)
    V(:,wellCond+1) = zeros(n, 1);
    V(:,1:wellCond) = V(:,1:wellCond)/R(1:wellCond, 1:wellCond);
    H(1:wellCond,1:wellCond-1) = R(1:wellCond,1:wellCond) * ...
        H(1:wellCond,1:wellCond-1) / R(1:wellCond-1, 1:wellCond-1);
    R(1:wellCond,1:wellCond) = eye(wellCond);
    R(:,wellCond+1) = zeros(m+1,1);
    whitened = true;

    for i = wellCond+1:m+1
        w = Afun(V(:,i-1));
        p = fwht(rademacher .* w);
        p = p(sample);
        
        x = zeros(m+1,1);
        for j = 1:nit
            % Richardson iteration
            x = x + SV'*(p - SV*x);
        end
        R(:,i) = x;
        q = w - V*x;
        s = fwht(rademacher .* q);
        s = s(sample);
        R(i,i) = norm(s);
        SV(:,i) = s/R(i,i);
        V(:,i) = q/R(i,i);
    end

    H = [H(:,1:wellCond-1),R(:,wellCond+1:m+1)];
    % Remember that the first vectors of the basis are not normalized!!
end
