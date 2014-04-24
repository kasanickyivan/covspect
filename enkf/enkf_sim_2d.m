function [X,Y] = enkf_sim_2d(Xi,Yi,nac,scn)
%   Function to simulate assimilation of some 1d model.
%
%   in:
%   Xi  :   initial perturbed ensemble
%   Yi  :   initial 'true' state
%   nac :   number of assimilation cycles
%   scn :   'scenario' for simulation, cell array, content
%           {1} :   mask matrix of the assimilated obseravations
%           {2} :   r - variance of the observations
%           {3} :   (function handle) function to evelve the state vector
%           {4} :   (function handle) function to augment the state 
%           {5} :   (function handle) function to reduce the state (remove 
%           augmented variable)
%           {6} :   (function handle) transform function
%           {7} :   (function handle) inverse transform function
%           {8} :   (function handle) weighted innovation function (see
%           enkf.m for deteaels)
%           {9} :   (cell of function handle) update function (see enkf.m
%           for details)
%
    M = scn{1};     %mask matirx
    r = scn{2};     %variance of observation    
    m_f = scn{3};   %model evolution function
    a_f = scn{4};   %augment function
    r_f = scn{5};   %reduce function
    t_f = scn{6};   %transform function
    it_f = scn{7};  %inverse transform function
    w_f = scn{8};   %weighted innovation function 
    u_fc = scn{9};   %cell array of update function
    

    [nx,ny,nvar,N] = size(Xi);
    X = zeros(nx,ny,nvar,N,nac+1);
    X(:,:,:,:,1)=Xi;
    Y = zeros(nx,ny,nvar,nac+1);
    Y(:,:,:,1)=Yi;
    
    
    for ac_ind = 2:nac+1
        % model advance 
        Y(:,:,:,ac_ind) = m_f(Y(:,:,:,ac_ind-1));
        for N_ind =1:N
            X(:,:,:,N_ind,ac_ind) = m_f(X(:,:,:,N_ind,ac_ind-1));
        end
        % data perturbation
        D = repmat(Y(:,:,1,ac_ind),1,1,N) + randn(nx,ny,N)*sqrt(r);
        % null, where is no observation
        D = D .* repmat(M,1,1,N);
        % augmentation of state vector
        XAG= a_f(X(:,:,:,:,ac_ind),M);
        % transformation
        XAG = t_f(XAG);                 % X augmented
        D = t_f(D);
        nvar = size(XAG,3);
        vd = repmat(nx*ny,nvar,1);          % variable dimensions  
        XP = reshape(XAG,nx*ny*nvar,N);     % X packed
        % pack observations
        DP = reshape(D,nx*ny,N);
        XP = enkf(XP,vd,@(x) w_f(x,DP,r),u_fc);
        XAG = reshape(XP,nx,ny,nvar,N);
        XAG = it_f(XAG);
        XA = real(r_f(XAG));
        X(:,:,:,:,ac_ind) = XA;
    end
end