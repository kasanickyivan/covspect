function d = cov_sample(X,Y)
% Returns sample covariance between 2 variables. 
% Replication in columns.
% 
    N = size(X,2);
    EX = mean(X,2);    
    % centering
    X = X - repmat(EX,1,N);
    EY = mean(Y,2);
    Y = Y - repmat(EY,1,N);
    d = 1/(N-1) * X*Y';
end