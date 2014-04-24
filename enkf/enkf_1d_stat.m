function [XM,bias,mae,rmse] = enkf_1d_stat(X,Y)
%   Computrte mean enseble member and some basic statistics for 1d state.
%
%   in:
%   X   :   ensemble
%   Y   :   reality
%
%   out:
%   XM  :   mean ensemble member
%   bias :  1/(nx*ny)*(X-Y)
%   mae  :  1/(nx*ny)*|X-Y|
%   rmse :  (1/(nx*ny)*(X-Y)^2)^.5
    
    [n,nvar,~,nac] = size(X);
    XM = squeeze(mean(X,3));
    
    R = reshape(XM,n,nvar,nac) -  reshape(Y,n,nvar,nac);
    
    bias = squeeze(mean(R,1));
    mae = squeeze(mean(abs(R),1));
    rmse = squeeze(mean(R.^2,1).^.5);
    
end