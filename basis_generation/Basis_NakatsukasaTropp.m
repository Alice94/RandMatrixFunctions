function [V, H] = Basis_NakatsukasaTropp(n, Afun, b, m, k)
if (~exist('k', 'var'))
    k = 2;
end
scaling = norm(b);
V = zeros(n, m+1);
V(:,1) = b / scaling;
H = zeros(m+1, m);

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
end