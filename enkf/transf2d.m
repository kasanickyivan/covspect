function X=transf2d(X,trans_f)
%Two dimensional transformation. 
%
%   function X=transf2d(X,trans_f)
%   X - 4 dimensional array (3rd dim variables, 4th dim replications)
%   trans_f - function handle - 1D transformation  

    [nx,ny,nvar,reps]=size(X);
    for rep_ind = 1:reps
        for var_ind = 1:nvar
            %apply transformation on columns
            for ny_ind = 1:ny
                X(:,ny_ind,var_ind,rep_ind) = ...
                    trans_f(squeeze(X(:,ny_ind,var_ind,rep_ind)));
            end
            %apply transformation on rows
            for nx_ind = 1:nx
                X(nx_ind,:,var_ind,rep_ind) = ...
                    trans_f(transpose(squeeze(X(nx_ind,:,var_ind,rep_ind))));
            end
        end
    end
    
end