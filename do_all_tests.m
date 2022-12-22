close all

% 1: Arnoldi (modified Gram-Schmidt)
% 2: Algorithm 3 -- Balabanov-Grigori + LSQR for the least-squares solve
% 3: Algorithm 4 -- Nakatsukasa-Tropp + Blendenpik
% 4: Algorithm 5 -- Nakatsukasa-Tropp without least-squares solve
% 5: sFOM from Guettel/Schweitzer preprint
% 6: Nakatsukasa-Tropp + whitening (i.e. switch to Algorithm 3 if things are bad)


%% Figure 1
% Convection-Diffusion equation (example from section 3.1 in
% [Botchev/Knizhnerman 2020], but a bit smaller).
N = 350;
Pe = 200; 
A = -ConvectionDiffusionExample(N+2,Pe);
n = size(A, 1);
vv = linspace(0,1,N);
vv = sin(pi*vv);
b = kron(vv, vv);
b = b'/norm(b);
[y_ex, ~] = phiRT_exp2(-A,b,1,1e-15,30);
Afun = @(x) A*x;
f = @(X, y) expm(full(X))*y;
ff = "exp";
ms = 50:50:500;
param.tol = 1e-9;                           % tolerance for quadrature rule
algorithms = [1, 1, 1, 1, 1, 0];

NAME = "ConvDiff1";
run_example;
PlotWithTime;

%% Figure 2
% Exponential of toeplitz matrix from option princing problem
n = 7000; 
[C, R] = exampleOptionPricing(n);
A = toeplitz(C, R);
Afun = @(x) ttimes(C, R, x);
f = @(X, y) expm(X)*y;
ff = "exp";
b = randn(n, 1);
b = b/norm(b);
y_ex = f(A, b);
ms = 400:400:3600;
param.tol = 1e-9;
algorithms = [1, 1, 1, 1, 1, 0]; 

NAME = "toeplitz";
run_example;
PlotWithTime;

%% Figure 3: 
% sqrt of Graph Laplacian of p2p-Gnutella08
load("p2p-Gnutella08.mat");
A = Problem.A;
A = diag(sum(A)) - A;
n = size(A, 1);
Afun = @(x) A*x;
b = randn(n, 1);
b = b/norm(b);
f = @(X, y) sqrtm(full(X))*y;
y_ex = f(A, b);
ms = 20:20:360;
algorithms = [1, 1, 1, 1, 1, 1];

NAME = "gnutella";
run_example;
PlotNoTime;

%% Figure 4
% Sign function
A = mmread("bfw782a.mtx");
n = size(A, 1);
Afun = @(x) A*x;
b = randn(n, 1);
b = b/norm(b);
f = @(X, y) sqrtm(full(X*X))\(X*y);
y_ex = f(A, b);
ms = 20:20:360;
algorithms = [1, 1, 1, 1, 1, 1];

NAME = "spectralProj";
run_example;
PlotNoTime;

%% Figure 5:
% inverse sqrt of 3D laplacian
N = 80;
T = spdiags(ones(N,1) * [-1, 2, -1], -1:1, N, N);
EYEN = speye(N);
A = kron(kron(T, EYEN), EYEN) + kron(kron(EYEN, T), EYEN) + kron(kron(EYEN, EYEN), T);
Afun = @(x) A*x;
n = size(A, 1);
b = randn(n,1); b = b / norm(b);
f = @(X, y) sqrtm(full(X))\y;
[V, H] = Basis_Arnoldi(n, Afun, b, 650);
mm = size(H,2);
y_ex = V(:,1:mm) * f(H(1:mm, :), [1; zeros(mm-1, 1)]);
ms = [50:50:300, 400:100:1000];
algorithms = [1, 1, 0, 0, 0, 0];

NAME = "invsqrt3DLap";
run_example;
PlotWithTime;

%% Figure 6: 
% Convection-diffusion from Guettel/Schweitzer paper
N = 100;
n = N^2;
L = spdiags(ones(n,1) * [-1, 2, -1], -1:1, N, N);
C = spdiags(ones(N,1) * [-1, 1, 0], -1:1, N, N);
h = 1/(N+1);
A2 = 1e-3 / h^2 * (kron(speye(N,N),L) + kron(L,speye(N,N))) ...
    + 1/h*(kron(speye(N,N), C) + kron(C, speye(N,N)));

load('data/convdiff_matrix.mat'); load('data/convdiff_sol.mat');
n = size(A,1); 
y_ex = ex_convdiff;

Afun = @(X) A*X;
f = @(X, y) sqrtm(full(X))\y;
b = ones(n,1);
b = b/norm(b);
[V, H] = Basis_Arnoldi(n, Afun, b, 800);
ms = 5:5:250;
algorithms = [1, 1, 1, 1, 1, 0];

NAME = "ConvDiff2";
run_example;
PlotWithTime;

%% Figure 7
% Exponential of wiki-vote
load('wiki-Vote.mat');
A = -Problem.A;
n = size(A, 1);
b = ones(n,1);
b = b / norm(b);
y_ex = expm(full(A))*b;
Afun = @(x) A*x;
f = @(X, y) expm(full(X))*y;
ms = 2:50;
algorithms = [1, 1, 1, 1, 1, 1];

NAME = "wikivote";
run_example;
PlotNoTime;




