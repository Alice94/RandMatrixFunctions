function [SV,SAV,Vfull] = bta_mod(A,v,m,t,hS)
% block truncated rational Arnoldi with randomized reduction
% A     - matrix
% v     - vector
% m     - dimension of the Krylov subspace
% t     - truncation parameter (t=2 will give Lanczos for Hermitian AB)
% hS    - handle to left-side basis reduction
% 
% SV    - left-reduced Krylov basis matrix
% SAV   - left-reduced A times basis matrix
% Vfull - also return full unreduced Krylov basis
% Slightly modified from bta.m from https://github.com/MarcelSchweitzer/sketched_fAb

N = length(v);

V = zeros(N, t); % truncated orthonormal basis (last t blocks)
v = v/norm(v);
V = [ V(:,2:end) , v ];
SV = hS(v); s = size(SV,1);
SV = [ SV, zeros(s,m) ];
SAV = zeros(s,m+1); 
Vfull = zeros(N, m+1);
Vfull(:,1) = v;
for j = 1:m
    w = V(:,end);
    Aw = A*w; 

    % compute these retrospectively
    SAV(:,j) = hS(Aw); 
    w = Aw;
    w = w - V*(V'*w);
    v = w/norm(w);
    V = [ V(:,2:end) , v ];
    SV(:,1+j) = hS(v); 
    Vfull(:,1+j) = v;
end

% final block column of SAV
SAV(:,1+j) = hS(A*v); 
end

