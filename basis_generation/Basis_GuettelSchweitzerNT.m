function [V, H, SV, SAV] = Basis_GuettelSchweitzerNT(n, Afun, b, m, k, s)
sample = randsample(n, s);
rademacher = sign(randn(n,1));

scaling = norm(b);
V = zeros(n, m+1);
V(:,1) = b / scaling;
H = zeros(m+1, m);
SV = zeros(s, m);
SAV = zeros(s, m);

for i = 1:m
    q = V(:,max(i-k+1,1):i); % local orthogonalization
    w = Afun(q(:,end));
    
    p = dct(rademacher .* w);
    p = p(sample);
    SAV(:,i) = p;
    p = dct(rademacher .* V(:,i));
    p = p(sample);
    SV(:,i) = p;

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
end