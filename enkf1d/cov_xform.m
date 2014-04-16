
%%
%   Function for estimating sample covariance in spectral space 
%   and reducing it by given mask. 
%
%   function Q=cov_xform(X,F,M)
%
%   X - enseble (ensemble members in columns)
%   F - transformation matrix
%   M - mask matrix
%%

function Q=cov_xform(XF,F,M)
    C=cov(XF');
    CF=F*C*F';
    Q=CF.*M;
    Q=real(F'*Q*F);
end