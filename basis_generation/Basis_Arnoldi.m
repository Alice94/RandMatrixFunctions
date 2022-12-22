function [V, H] = Basis_Arnoldi(n, Afun, b, m)
beta = norm(b);
V = zeros(n, m+1);
H = zeros(m+1, m);
V(:,1) = b / beta;

for i = 1:m
    w = Afun(V(:,i));
    h = V'*w;
    H(:, i) = h;
    w = w - V*h;
    
    beta2 = norm(w);
    w = w/norm(w);
    V(:,i+1) = w;
    H(i+1, i) = beta2;
    
    if (abs(beta2) < 1e-14)
        V = V(:,1:i+1);
        H = H(1:i+1, 1:i);
        break;
    end
end
