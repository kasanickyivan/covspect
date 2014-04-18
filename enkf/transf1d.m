function X=transf1d(X,trans_f)
%   transform columns using trans_f function
%
%   in:
%   X   :   2d array
%   trans_f     :   (function handle) transfomration function
    [n,N]=size(X);
    for N_ind = 1:N
        X(:,N_ind) = trans_f(X(:,N_ind));
     end
end