function Y=reds1d(X)
%   Reduce augmented state by dropping first variable.
%   
%   in:
%   X   : ensemble 
%   
%   out:
%   Y   : ensemble with dropped first variable

    nvar = size(X,2);
    Y = X(:,2:nvar,:);
end