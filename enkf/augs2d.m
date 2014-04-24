function Y=augs2d(X,M)
%   Function augment the 2 dimensional state (1st variable) ensemble,
%   and puts zeroes where no observations.
%   
%   in:
%   X   : ensemble 
%   M   : mask vector (1 where observation is, 0 where is not)
%
%   out:
%   Y   : ensemble with augmented first variable.
%
    [nx,ny,nvar,N] = size(X);
    Y = zeros(nx,ny,nvar+1,N);
    Y(:,:,2:nvar+1,:) = X;
    Y(:,:,1,:) = squeeze(X(:,:,1,:)) .* repmat(M,1,1,N);
end