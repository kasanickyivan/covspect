%   Counts covariances, corelations, local covariancies and local
%   correlations for 4 dimensional array, where 
%  
%       - 1st and 2nd dimension is gridpoint position
%               - 3rd is variable index
%               - 4th is replications index
%
function [cov,loc_cov,cor,loc_cor] = cov_cor_2d(Y)
    [n,~,nvar,reps] = size(Y);
    Y = pack_state(Y);
    % 
    Y = Y - repmat(mean(Y,2),1,reps);
    var_index = cell(4,1);
    for var_ind = 1:nvar
        var_index{var_ind} = (1:n*n)+((var_ind-1)*n*n);
    end
    
    
    cov = zeros(nvar,nvar,n*n,n*n);
    for var_ind_1 = 1:nvar
        for var_ind_2 = 1:nvar
            cov(var_ind_1,var_ind_2,:,:)= ...
                1/(reps-1)*Y(var_index{var_ind_1},:)*Y(var_index{var_ind_2},:)'; 
        end
    end

    cor = zeros(nvar,nvar,n*n,n*n);
    for var_ind_1 = 1:nvar
        for var_ind_2 = 1:nvar
            var1 = diag(squeeze(cov(var_ind_1,var_ind_1,:,:)));
            var2 = diag(squeeze(cov(var_ind_2,var_ind_2,:,:)));
            cor(var_ind_1,var_ind_2,:,:) = ...
               squeeze(cov(var_ind_1,var_ind_2,:,:)) ./ (var1*var2').^(.5);
        end
    end

    loc_cov = zeros(nvar,nvar,n,n,n,n);
    loc_cor = zeros(nvar,nvar,n,n,n,n);

    for x_ind = 1:n
        for y_ind =1:n
            xy_gp_ind = coor_grid2pack(x_ind,y_ind,n);
            for var_ind_1 = 1:nvar
                for var_ind_2 = 1:nvar
                    c12 = squeeze(cov(var_ind_1,var_ind_2,xy_gp_ind,:));
                    loc_cov(var_ind_1,var_ind_2,x_ind,y_ind,:,:) = reshape(c12,n,n);
                    c12 = squeeze(cor(var_ind_1,var_ind_2,xy_gp_ind,:));
                    loc_cor(var_ind_1,var_ind_2,x_ind,y_ind,:,:) = reshape(c12,n,n);    
                end
            end
        end
    end
end