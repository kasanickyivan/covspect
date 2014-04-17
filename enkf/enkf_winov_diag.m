function A=enkf_winov_diag(X,D,r)
% Function compute weighted inovation using diagonal elemens of sample 
% sample cavariance as a approximation to forecast covariation.
%   
% A = (E + rI)^(-1)(X-D)
%   E = sample covariance with zeros outside main diagonal
%   r - variance of obseravtions
%   D - perturbed data
%

    N = size(X,2);
    % innovation
    IN = X - D;
    % diagonal elements of the covariance of the observed variable
    e = cov_diag(X,X);
    % innovation weigthed by (P + R)^1 
    A = repmat((e + r).^(-1),1,N) .* IN;

end