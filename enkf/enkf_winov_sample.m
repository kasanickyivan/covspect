function A=enkf_winov_sample(X,D,r)
% Function compute weighted invation using  
% sample cavariance as a approximation to forecast covariation.
%   
% A = (C + rI)^(-1)(X-D)
%   C = sample covariance 
%   r - variance of obseravtions
%   D - perturbed data
%

    n = size(X,1);
    IN = X - D;
    % sample covariance
    C = cov_sample(X,X);
    % innovation weigthed by (P + R)^1 
    A = (C+r*eye(n))\IN;

end