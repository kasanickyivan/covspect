function XA=enkf_update_diag(XF,XO,A)
% Function compute analysis from the forecast using diagonal approximation
% of cross-covarinace between observed ann updated variable.
% 
%       XA = XF - E_{FO}.A
%       A = (Q + R)^(-1)(XO-D)
%   
%   XF - forecast state of the variable to be updated
%   XO - forecast state of the observed variable
%   A -  weighted inovation (see enkf_winov_diag.m)
 
    N = size(XF,2);
    e = cov_diag(XF,XO);
    XA = XF - repmat(e,1,N) .* A;
    
end