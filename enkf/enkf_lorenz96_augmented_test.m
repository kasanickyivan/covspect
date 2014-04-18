function[BIAS,MAE,RMSE]=enkf_lorenz96_augmented_test(m,n,N,ipv,r,nac,reps,ts)
%   Function for testing diferent estimates of forecast matrix in EnKF
%   usinf lorenz96 model.
%   
%   in:
%   m   :   mask vector
%   n   :   length of state
%   N   :   number of ensemble members
%   nac :   number of assimilation cycles
%   reps:   number of replications
%   ipv :   initial perturbation variance
%   r   :   variance of observations
%   ts  :   number of time steps between assimilations
    
    %arguments lorenz96 function
    F=8;kappa=1;dt=0.01;
    l_f = @(u,ts) lorenz96(u,dt,F,kappa,ts);

    qmf=MakeONFilter('Coiflet',2);
    L=4;

    % diferent scenarios
    % assimilation function
    % transformation scenarios
    
    sc_t = {@(x) x,...
            @(x) n^(-.5)*fft(x),...
            @(x) transf1d(x,@(y)FWT_PO(y,L,qmf))};%,...
%             @(x) n^(-.5)*fft(x), ...
%             @(x) transf1d(x,@(y)FWT_PO(y,L,qmf))};
    sc_it = {@(x) x, ...
             @(x) n^(.5)*ifft(x), ...
             @(x) transf1d(x,@(y) IWT_PO(y,L,qmf))};%;,...
%              @(x) n^(.5)*ifft(x), ...
%              @(x) transf1d(x,@(y) IWT_PO(y,L,qmf))};
    % wegihted innovation scenarios
    sc_w = {@(x,y) enkf_winov_sample(x,y,r),...
            @(x,y) enkf_winov_diag(x,y,r),...
            @(x,y) enkf_winov_diag(x,y,r)};%,...
%             @(x,y) enkf_winov_diag(x,y,r),...
%             @(x,y) enkf_winov_diag(x,y,r)};
    % update scenarios
    sc_u= { {@(x,y,z) enkf_update_sample(x,y,z),@(x,y,z) enkf_update_sample(x,y,z)},...
            {@(x,y,z) enkf_update_diag(x,y,z),@(x,y,z) enkf_update_diag(x,y,z)},...
            {@(x,y,z) enkf_update_diag(x,y,z),@(x,y,z) enkf_update_diag(x,y,z)}};%,...
%             {@(x,y,z) enkf_update_diag(x,y,z),@(x,y,z) enkf_update_sample(x,y,z)},...
%             {@(x,y,z) enkf_update_diag(x,y,z),@(x,y,z) enkf_update_sample(x,y,z)} };
   
    sc_names = {'Sample','FFT d CC','DWT d CC'};

    rmse = zeros(length(sc_t),nac+1,reps);
    mae = zeros(length(sc_t),nac+1,reps);
    bias = zeros(length(sc_t),nac+1,reps);


    for rep_ind = 1:reps
        %initialization:
        Xi = rand(n,1)*0.001-0.0005;
        Yi = rand(n,1)*0.001-0.0005;
        % initial number of time steps 
        ints = 1000;
        Xi = l_f(Xi,ints);
        Yi = l_f(Yi,ints);
        Xens=repmat(Xi,1,N) + randn(n,N)*sqrt(ipv);
        for sc_ind = 1:length(sc_t)
            [Y,X]=enkf_lorenz_augmented(Yi,Xens,m,r,ts,nac,l_f,sc_t{sc_ind},sc_it{sc_ind},sc_w{sc_ind},sc_u{sc_ind});
            XA = squeeze(mean(X,2));
            R = Y - XA;
            bias(sc_ind,:,rep_ind) = mean(R,1);
            mae(sc_ind,:,rep_ind) = mean(abs(R),1);
            rmse(sc_ind,:,rep_ind) = mean(R.^2,1);
        end
        fprintf('.');
    end
    fprintf('\n');
    PX = repmat((0:nac)',1,length(sc_t));
    BIAS = mean(bias,3);
    figure()
    plot(PX,BIAS');title(sprintf('%s %.f','BIAS augmented',N));
    legend(sc_names);
    xlabel('assimilation steps');
    RMSE = mean(rmse,3);
    figure();
    plot(PX,RMSE');title(sprintf('%s %.f','RMSE augmented',N));
    legend(sc_names);
    xlabel('assimilation steps');
    MAE = mean(bias,3);
    figure();
    plot(PX,MAE');title(sprintf('%s %.f','MAE augmented',N));
    xlabel('assimilation steps');
    legend(sc_names);
end