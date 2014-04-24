function [XM,bias,mae,rmse] = enkf_2d_stat(X,Y)
%   Computrte mean enseble member and some basic statistics.
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
    
    [nx,ny,nvar,~,nac] = size(X);
    XM = squeeze(mean(X,4));
    
    R = reshape(XM,nx*ny,nvar,nac) -  reshape(Y,nx*ny,nvar,nac);
    
    bias = squeeze(mean(R,1));
    mae = squeeze(mean(abs(R),1));
    rmse = squeeze(mean(R.^2,1).^.5);
    
end