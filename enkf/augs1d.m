function Y=augs1d(X,m)
%   Function augment the state (1st variable) ensemble, and puts zeroes where no
%   observations.
%   
%   in:
%   X   : ensemble 
%   m   : mask vector (1 where observation is, 0 where is not)
%
%   out:
%   Y   : ensemble with augmented first variable.
%
    [n,nvar,N] = size(X);
    Y = zeros(n,nvar+1,N);
    Y(:,2:nvar+1,:) = X;
    Y(:,1,:) = squeeze(X(:,1,:)) .* repmat(m,1,N);
end