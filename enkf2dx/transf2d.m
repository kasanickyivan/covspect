%%
%Two dimensional transformation. 
%
%   function X=transf2d(X,trans_f)
%   X - 2 dimensional array
%   trans_f - function handle - 1D transformation  
function X=transf2d(X,trans_f)
    [nx,ny]=size(X);
    %apply transformation on columns
    for ny_ind = 1:ny
        X(:,ny_ind) = trans_f(squeeze(X(:,ny_ind)));
    end
    %apply transformation on rows
    for nx_ind = 1:nx
        X(nx_ind,:) = trans_f(squeeze(X(nx_ind,:)));
    end
end