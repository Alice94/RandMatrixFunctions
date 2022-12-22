function y = ttimes(c,r,x)
%TTIMES	Toeplitz matrix product.
%   y = TTIMES(c,r,x) computes y=Tx, where c and r are the first row and
%   column of the rectangular Toeplitz matrix T.
%
%   See also cltimes, toeplitz.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 3, 2010

if nargin~=3, error('ttimes:nargin','Wrong number of arguments.'); end
if ~isvector(c) || ~isvector(r), error('ttimes:vector','c and r must be vectors'); end
if c(1)~=r(1), warning('ttimes:diagconflict','col wins'); end

c = c(:);
r = r(:);
m = size(c,1);
n = size(r,1);
if n~=size(x,1), error('ttimes:size','Inner dimensions must agree.'); end

% two choices for the size of a circulant matrix C such that C(1:n,1:n)==T.
%N = pow2(nextpow2(m+n-1)); % N is the minimal power of 2, so k=N-(2n-1).
N = m+n-1; % N is the minimal integer value, so k=0.

k = N-(m+n-1); % k = # of extra blank diagonals in C

% y = T*x  via  [y;?] = C * [x;0], where C = [T ?;? ?]
t = [ c; zeros(k,1); r(end:-1:2) ]; % t = C(:,1)

% Don't invoke sparse machinery, assume sufficient memory. RL 2016.
%y = ifft(full( spdiags(fft(t),0,N,N) * fft([x;zeros(N-n,size(x,2))]) ));
y = ifft( repmat(fft(t),1,size(x,2)) .* fft([x;zeros(N-n,size(x,2))] ));


y = y(1:m,:); % take Tx by [y;?]=C*[x;0]

% output is real/integer if input is real/integer
if isreal(t) && isreal(x), y=real(y); end
if all(t==round(t)) && all(all(x==round(x))), y=round(y); end
