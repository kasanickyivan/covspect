function A=enkf_winov_sample(X,H,D,r)
% Function compute weighted invation using  
% sample cavariance as a approximation to forecast covariation.
%   
%   A = H'(HCH' + rI)^(-1)(HX-HD)
%   H = observation operator
%   C = sample covariance 
%   r - variance of obseravtions
%   D - perturbed data
%

    %n = size(X,1);
    IN = H*X - H*D;
    % sample covariance
    C = cov_sample(X,X);
    HCH = H*C*H';
    % innovation weigthed by (P + R)^1 
    A = (HCH+r*eye(size(H,1)))\IN;
    A = H'*A;
end