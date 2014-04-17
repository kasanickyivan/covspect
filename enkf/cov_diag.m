function d = cov_diag(X,Y)
% Returns diagonal elements of sample covariance between 2 variables (must have same dimension). 
% Replication in columns.
% 
    N = size(X,2);
    EX = mean(X,2);    
    % centering
    X = X - repmat(EX,1,N);
    EY = mean(Y,2);
    Y = Y - repmat(EY,1,N);
    d = 1/(N-1) * sum(X.*conj(Y),2);
end