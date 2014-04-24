function Y=reds2d(X)
%   Reduce augmented state of 2d state by dropping first variable.
%   
%   in:
%   X   : ensemble 
%   
%   out:
%   Y   : ensemble with dropped first variable
    nvar = size(X,3);
    Y = X(:,:,2:nvar,:);
end