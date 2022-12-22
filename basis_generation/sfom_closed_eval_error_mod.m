function err = sfom_closed_eval_error_mod(V,SV,SAV,Sv,R,ex,f,m)
% Assumes that SV has orthonormal columns
SVm = SV(:,1:m);
M = SVm'*SVm;

coeffs = M\( f( (SVm'*SAV(:,1:m))/M , (SVm'*Sv)) );
appr = V(:,1:m)*(R(1:m,1:m)\coeffs);
err = norm(appr - ex)/norm(ex);