
function XA = enkf(XF,var_dims,winov_f,update_f)
% Function for EnKF update with observed full one (1st) variable. 
% 
% Kalman formula for i-th variable and each ensemble:
%
%   X^a_i = X_i - P_{i1}(P{11}+R)^{-1}(X_1-D)
%         = X_i - P_{i1}A, where A = (P{11}+R)^{-1}(X_1-D)
%       P_{ij} covariance between i-th and j-th variable
%       R data error matrix%     
%
% Inputs:
%   XF - forecast ensemble - each column is the full state
%   var_dims - dimension, e.g. number of rows, of each variable, so  sum(var_dims) == size (XF,1)
%   winov_f - function handle - fuction that creates weighted innovations,
%             returnig matrix A
%   update_f - cell function handle - functions that make update to one variable,
%              when given XF and A, ith cell function to update ith variable  
%
 
    % variable index to
    vit = cumsum(var_dims);
    % variable index from
    vif = [1 vit(1:length(vit)-1)+1]; 
    
    % analysis enseble 
    XA = zeros(size(XF));
    % Kalman filter update
    i = vif(1):vit(1); % indexes of 1st (observed) variable
    % weighted innovations 
    A = winov_f(XF(i,:));
    % updating all variables
    for var_ind = 1:length(var_dims)
        j = vif(var_ind):vit(var_ind); % indexed of variables to be updated
        XA(j,:) = update_f{var_ind}(XF(j,:),XF(i,:),A);
    end
end