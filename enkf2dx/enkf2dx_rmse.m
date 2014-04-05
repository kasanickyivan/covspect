%%
%Function count RMSE of analysis. 
%
%   function RMSE = enkf2dx_rmse(XA,X)
%   XA - analysis
%   X - true state
function RMSE = enkf2dx_rmse(XA,X)
    [nx,ny,nvar,r] = size(XA);
    
    XA = reshape(XA,nx*ny,nvar,r);
    X = reshape(X,nx*ny,nvar,r);
   
    RMSE = sqrt(mean((XA-X).^2,1));
end