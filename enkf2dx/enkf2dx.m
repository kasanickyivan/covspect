%%
%   EnKF with observation of entire one variable. 
%
%   Implemenation based on J. Mandel: A  Brief Tutorial on 
%   the Ensemble Kalman Filter 
%   http://arxiv.org/abs/0901.3725
%   
%   Full algorithm (protected access): 
%   http://mathweb.ucdenver.edu/~wikiuser/grantwiki/doku.php?id=spectral_enkf_for_multidimensional_data
%
%   function XA=enkf2dx(XF,d,r,i,trans_f,itrans_f)
%
%   XF - forecast ensemble
%   d - observed data (matrix)
%   r - variance of observed data - covariance matrix of observed data is
%   then r*I
%   i - index of observed variable
%   trans_f - function handle -  2D transformation  
%   itrans_f - function handle - 2D inverse transformation
%  
%%
function XA=enkf2dx(XF,d,r,i,trans_f,itrans_f)
    % nx,ny - grid dimension
    % N - number of ensemble members
    [nx,ny,m,N] = size(XF);
    

    X = zeros(size(XF));
    % transform ensemble to spectral space
    for ens_ind = 1:N
        for m_ind = 1:m
            X(:,:,m_ind,ens_ind)=trans_f(squeeze(XF(:,:,m_ind,ens_ind)));
        end
    end
    X = pack_state(X);
    
    % transform, pack and perturbate data to spectral space
    D = repmat(pack_state(trans_f(d)),1,N)+randn(nx*ny,N)*sqrt(r);
    
    % i-th variable row from
    ri_f = (i-1)*nx*ny+1; 
    % i-th variable row to
    ri_t = i*nx*ny;
    
    % diagonal approximation of covariance matrix
    %E = block_cov_vec(X,m,i);
    E = mvbdca(X,nx*ny,i);
    E_i = E(ri_f:ri_t);
    
    
    % inovation
    INOV = X(ri_f:ri_t,:) - D;
    A = INOV .* repmat((E_i + r).^(-1),1,N);
    B = repmat(A,m,1) .* repmat(E,1,N);
    X = X - B;
    
    % unpack and transform back to spatial space
    XA = unpack_state(X,nx,ny);
    for ens_ind = 1:N
        for m_ind = 1:m
            XA(:,:,m_ind,ens_ind)=itrans_f(squeeze(XA(:,:,m_ind,ens_ind)));
        end
    end
    
    XA = real(XA);
end