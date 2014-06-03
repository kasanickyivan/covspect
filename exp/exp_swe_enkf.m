function exp_swe_enkf(N,ts,no)
%   Experiments with diferent covariance approximation using diffent
%   approximation of covariance matrices.
%   in :
%   N   :   number of ensembles
%   ts  :   time step between assimilations
%   no  :   number of observations, observations are no matrix no x no
%
%
%   SWE model is assumed to be a symplification of real atmospheric flow.
%   The dimensions are set so that the size of the area is simimar to the
%   area of a wholw globe. 

    reps = 50;
    n = 64;        %length of state vector          
    r = 100;        %variance of the observations
    r_pert = 1000;   % variance of initial perturbation
    M = zeros(n);
    M(1:no,1:no)=1;     % mask matrix  
    nac = 25;       % number of assimilation cyccles


    %arguments swe function
    dt=1;
    dx=150000;dy=150000;
    %initial condition to swe
    ih = 10000; %initial water height (water level)
    dw = 16; %width of drop at begining
    dh = 1000; % height of intitial drop
    mbd = 0; % minimal boundary distance 

    %argument to Coiflets 
    qmf=MakeONFilter('Coiflet',2);
    L=4;


    scn = cell(3);
    scn_names = {'FFT','DWT','EnKF'};%,'EnKF, H=I'};

    scn{1} =cell(9);
    % Transoformation to fourier space and diagonal approximation.
    scn{1}{1} = M;                                         %mask vector
    scn{1}{2} = r;                                         %variance of observation
    scn{1}{3} = @(x) waterwave2(x,dt,dx,dy,ts);            %model evolution function
    scn{1}{4} = @(x,m) augs2d(x,m);                        %augment function
    scn{1}{5} = @(x) reds2d(x);                            %reduce function
    scn{1}{6} = @(x) transf2d(x,@(y) n^(-.5)*fft(y));      %transform function
    scn{1}{7} = @(x) transf2d(x,@(y) n^(.5)*ifft(y));      %inverse transform function
    scn{1}{8} = @(x,d,r) enkf_winov_diag(x,d,r);             %weighted innovation function 
    scn{1}{9} = {@(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a)};        %cell array of update function

    % Transoformation to wavelet space and diagonal approximation.
    scn{2} = cell(9);
    scn{2}{1} = M;                                         
    scn{2}{2} = r;                                         
    scn{2}{3} = @(x) waterwave2(x,dt,dx,dy,ts);            
    scn{2}{4} = @(x,m) augs2d(x,m);                        
    scn{2}{5} = @(x) reds2d(x);                            
    scn{2}{6} = @(x) transf2d(x,@(y) FWT_PO(y,L,qmf));     
    scn{2}{7} = @(x) transf2d(x,@(y) IWT_PO(y,L,qmf));     
    scn{2}{8} = @(x,d,r) enkf_winov_diag(x,d,r);           
    scn{2}{9} = {@(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a)};
    % Standart EnKF with full observations.   
    h = reshape(M,[n*n 1]);
    H = sparse(diag(h));
    Hnull = all(H==0,2);
    H = H(~Hnull,1:n*n);
    
    scn{3}{1} = M;                                         
    scn{3}{2} = r;                                         
    scn{3}{3} = @(x) waterwave2(x,dt,dx,dy,ts);            
    scn{3}{4} = @(x,m) x;                        
    scn{3}{5} = @(x) x;                            
    scn{3}{6} = @(x) x;     
    scn{3}{7} = @(x) x;     
    scn{3}{8} = @(x,d,r) enkf_winov_sample(x,H,d,r);           
    scn{3}{9} = {@(x,o,a) enkf_update_sample(x,o,a),...
                 @(x,o,a) enkf_update_sample(x,o,a),...
                 @(x,o,a) enkf_update_sample(x,o,a),...
                 @(x,o,a) enkf_update_sample(x,o,a)};     


         

    BIAS = zeros(length(scn),reps,3,nac+1);
    MAE = zeros(length(scn),reps,3,nac+1);
    RMSE = zeros(length(scn),reps,3,nac+1);

    for rep_ind = 1:reps
        % initializationa
        U = init_swe(n,ih,dw,dh*(rand()+1),mbd,dt,dx,dy);
        Xi = zeros(n,n,3,N);
        Xi(:,:,:,:) = repmat(U,[1 1 1 N]) + randn(n,n,3,N)*sqrt(r_pert);
        Yi = init_swe(n,ih,dw,dh,mbd,dt,dx,dy); 
        fprintf('%g/%g',rep_ind,reps);
        for scn_ind = 1:length(scn)
            [X,Y] = enkf_sim_2d(Xi,Yi,nac,scn{scn_ind}); 
            [~,bias,mae,rmse]=enkf_2d_stat(X,Y);
            BIAS(scn_ind,rep_ind,:,:) = bias;
            MAE(scn_ind,rep_ind,:,:) = mae;
            RMSE(scn_ind,rep_ind,:,:) = rmse;
            fprintf('.');
        end
        fprintf('\n');
    end
    
    clear bias rmse mae

    bias = squeeze(mean(BIAS,2));
    mae = squeeze(mean(MAE,2));
    rmse = squeeze(mean(RMSE,2));
    
    x = repmat((0:nac)',1,length(scn));
    
    for var_ind = 1:3
        subplot(131);
        plot(x,squeeze(bias(:,var_ind,:))')
        legend(scn_names);
        xlabel('Assimilation step');
        title('Bias');

        subplot(132);
        plot(x,squeeze(mae(:,var_ind,:))')
        legend(scn_names);
        xlabel('Assimilation step');
        title('MAE');


        subplot(133);
        plot(x,squeeze(rmse(:,var_ind,:))')
        legend(scn_names);
        xlabel('Assimilation step');
        title('RMSE');

         axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
        'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        main_title = sprintf('SWE - variable no.%g - n %g - N %g - no. of obs. p %g - ts %g',...
                        var_ind,n,N,no,ts);
        text(0.5, 1,main_title,...
        'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
    
        img_file = sprintf('swe_v%g_n%g_N%g_no%g_ts%g',var_ind,n,N,no,ts);

        %savefig(['img/' img_file '.fig']);
        print('-dpng', ['img/' img_file '.png']);
    end
    
    file_out_mat = sprintf('mat/swe_n%g_N%g_no%g_ts%g.mat',n,N,no,ts); 
    save(file_out_mat,'n','N','no','ts','BIAS','RMSE','MAE'); 

end
% 


