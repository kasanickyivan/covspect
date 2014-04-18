function [Y,XA]=enkf_lorenz_augmented(Yi,Xi,m,r,ts,nac,l_f,t_f,it_f,winov_f,update_f)
%   EnKF simulation using lorenz model.
%   
%   in:
%   Yi  :   initial observation
%   Xi  :   initial ensemble 
%   m   :   mask vector
%   r   :   variance of observed data
%   ts  :   time steps between each assimilation cycle, vector of length nac, 
%           between i an i+1 analysis the model is evolved usinngs ts(i) time
%           steps, if size differs, only first values is used
%   nac :   number of assimilation cycles
%   l_f :   (function handle) function to evolve model
%   t_f :   (function handle) transform function
%   it_f:   (function handle) inverse transform function
%   winov_f :   (function handle) function creating weighted innovations
%               see enkf function fore details
%   update_f:   (function handle) function creating update to variables,
%               see enkf function to details
%
%   out:
%   Y   :   (size(Y)=[n,nac+1])the "true" state, dim n x nac + 1
%   XA  :   (size(XA)=[n,N,nac+1]) analysis ensemble, 
%           (:,:,1) is initial perturbed ensemble 
    if length(ts) ~= nac
        ts = repmat(ts(1),1,nac);
    end
    [n,N] = size(Xi);
    
    XA = zeros(n,N,nac+1);
    Y = zeros(n,nac+1);
    %mask matrix
    M = repmat(m,1,N);
    Xaug = zeros(n,2,N);
    XA(:,:,1) = Xi;
    Y(:,1) = Yi;
    for ac_ind = 2:nac+1
       Y(:,ac_ind) = l_f(Y(:,ac_ind-1),ts);
       for N_ind = 1:N
            XA(:,N_ind,ac_ind) = l_f(XA(:,N_ind,ac_ind-1),ts(ac_ind-1));
       end
       D = repmat(Y(:,ac_ind),1,N)+randn(n,N)*sqrt(r);
       %augment state  
       Xaug(:,1,:) = XA(:,:,ac_ind);
       Xaug(:,2,:) = XA(:,:,ac_ind); 
       D = (D.*M) + (squeeze(Xaug(:,1,:)).*(1-M));
       D = t_f(D);
       for var_ind = 1:2
            Xaug(:,var_ind,:) = t_f(squeeze(Xaug(:,var_ind,:)));
       end
       % pack state 
       X = reshape(Xaug,2*n,N);
       X = enkf(X,[n n] ,@(x) winov_f(x,D), update_f );
       j = n+1:2*n;
       % transform just the second variable 
       XA(:,:,ac_ind) = real(it_f(squeeze(X(j,:))));
    end
end