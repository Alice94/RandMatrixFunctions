function [col, row] = exampleOptionPricing(n)
% [c,r] = exampleOptionPricing(n)
%
% Example 3 from Lee, Hong-Kui Pang, and Hai-Wei Sun. Shift-invert Arnoldi
% approximation to the Toeplitz matrix exponential. SIAM J. Sci. Comput.,
% 32(2):774-792, 2010
%
% Time steps for expm(t*T) used in their work are t=0.5 and t=1.0

% Parameters used in their experiment
xi_min = -2;
xi_max = 2;
Delta_xi = (xi_max - xi_min) / (n+1);
K = 100;
nu = 0.25;
r = 0.05;
lambda = 0.1;
mu = -0.9;
sigma = 0.45;

kappa = exp( mu + sigma^2/2) - 1;

% TODO use Matlab's built in Gaussian pdf functionality
theta = @(eta) exp( - (eta-mu).^2 / (2*sigma^2) ) / (sqrt(2 * pi) * sigma);

% Tridiagonal matrix corresponding to central differences
mid = nu^2/Delta_xi^2;
off = (2*r - 2*lambda*kappa - nu^2)/(4*Delta_xi);

% Integral part of the equation
Ic = theta( ( 0:-1:(1-n) )' * Delta_xi );
Ir = theta( ( 0:(n-1) )  * Delta_xi );

% Assemble total matrix and rhs
col = [-mid - r - lambda; mid/2 - off; zeros(n-2, 1)] + lambda * Delta_xi * Ic;
row = [-mid - r - lambda; mid/2 + off; zeros(n-2, 1)]' + lambda * Delta_xi * Ir;
end
