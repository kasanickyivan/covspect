%   Assimilation run for 2d model on square grid (n x n)
%
%   
%   Y - true state
%   X - initial ensemble
%   m_f - model advance (function handle)
%   ass_f - method of assimilation (function handle)

function XA = assim2d(ens,obs,m_f,ass_f)
    %   nx,ny - grid size
    %   nvar - number of variables
    %   rl - run length (number of assimilation cycles) 
    [nx,ny,nvar,N] = size(ens);
    rl = size(obs,3);
    XA = zeros(nx,ny,nvar,rl);
    fprintf('assimilation run, %g cycles\n',rl);
    for rl_ind = 1:rl
        fprintf('.')
        ens = ass_f(ens,squeeze(obs(:,:,rl_ind)));
        XA(:,:,:,rl_ind)=unpack_state(mean(pack_state(ens),2),nx,ny);
        for ens_ind = 1:N
            ens(:,:,:,ens_ind) = m_f(squeeze(ens(:,:,:,ens_ind)));
        end 
    end 
    fprintf('\n');
end
