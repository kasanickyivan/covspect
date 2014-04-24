function X=transf1d(X,trans_f)
%   %Two dimensional transformation. 
%
%   in:
%   X - 3 dimensional array (2nd dim variables, 3rd dim replications)
%   trans_f - function handle - 1D transformation  
    
    [~,nvar,reps]=size(X);
    for rep_ind = 1:reps
        for var_ind = 1:nvar
            X(:,var_ind,rep_ind) = trans_f(X(:,var_ind,rep_ind));
        end
     end
end