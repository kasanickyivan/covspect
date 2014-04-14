% generate random random ensemble of state of shallow wate equations
%
% input
%   n - dimension of grid
%   reps - number of replications
%   h - initial water level 
%   dw_min,dw_max,dh_min,dh_max - argumenst to add_rand_drop
%   ts - time step, argument for waterwave2 function
%   init_ts - initial time step
%   init_d - initial nuber of drops
%   dper - period for adding drops
%
%
% output
%   Y - 4d array dimension [x,y,variable,replication]

function Y=generate_waterwave(n,reps,h,dw_min,dw_max,dh_min,dh_max,ts,init_ts,init_d,dper)

    nvar=3; %number of variables foe EnKF
    dx=1;   %argument to waterwave2
    dy=1;   %argument to waterwave2
    dt=0.01;%argument to waterwave2

    Y=zeros(n,n,nvar,reps);Y(:,:,1,1)=ones(n,n)*h;
    for init_ind = 1:init_d
        Y(:,:,1,1)=add_rand_drop(squeeze(Y(:,:,1,1)),dw_min,dw_max,dh_min,dh_max);
        Y(:,:,:,1)=waterwave2(squeeze(Y(:,:,:,1)),dt,dx,dy,init_ts);
    end
    
    fprintf('.');
    for rep_ind = 2:reps
        Y(:,:,:,rep_ind) = waterwave2(squeeze(Y(:,:,:,rep_ind-1)),dt,dx,dy,ts);
        fprintf('.');
        % adding drop
        if(mod(rep_ind,dper)==0)
            Y(:,:,1,rep_ind)=add_rand_drop(squeeze(Y(:,:,1,rep_ind)),dw_min,dw_max,dh_min,dh_max);
            Y(:,:,:,rep_ind) = waterwave2(squeeze(Y(:,:,:,rep_ind)),dt,dx,dy,ts);
            fprintf('\n')
        end
    end
    fprintf('\n');
end