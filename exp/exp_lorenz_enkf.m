function  exp_lorenz_enkf(n,N,no)
%   Experiments file - assimilation to lorenz 96 model, 3 scenarios
%   (assimilation only observations on subgrid with augmenting the state
%   and estimating the covariance by diagonal of spectral sample covaiance
%   - using FFT, DWT, and standart EnKF with fully observed state).
%
%
%       n   :  length of state vector
%       N   :  number of ensebles 
%       no  :  number of observations
%

    reps=50;
    ts = 100;       %time step between assismilations
    r = 0.1;        %variance of the observations
    m = zeros(n,1);
    m(1:no)=1;     % mask vector  
    nac = 10;       % number of assimilation cycles

    %arguments lorenz96 function
    F=8;kappa=1;dt=0.01;
    %argument to Coiflets 
    qmf=MakeONFilter('Coiflet',2);
    L=4;

    %   Scenarios
    scn = cell(3);
    scn_names = {'FFT' 'DWT' 'EnKF, H=I'};
    scn{1} =cell(9);
    % Diagonal approximation in spectral space with diagonal cross-covariances.
    scn{1}{1} = m;                                         %mask vector
    scn{1}{2} = r;                                         %variance of observation
    scn{1}{3} = @(u) lorenz96(u,dt,F,kappa,ts);            %model evolution function
    scn{1}{4} = @(x,m) augs1d(x,m);                        %augment function
    scn{1}{5} = @(x) reds1d(x);                            %reduce function
    scn{1}{6} = @(x) transf1d(x,@(y) n^(-.5)*fft(y));      %transform function
    scn{1}{7} = @(x) transf1d(x,@(y) n^(.5)*ifft(y));      %inverse transform function
    scn{1}{8} = @(x,d,r) enkf_winov_diag(x,d,r);             %weighted innovation function 
    scn{1}{9} = {@(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a)};        %cell array of update function
    % Diagonal approximation in Coiflet space with diagonal cross-covariances.
    scn{2} = cell(9);
    scn{2}{1} = m;                                         
    scn{2}{2} = r;                                         
    scn{2}{3} = @(u) lorenz96(u,dt,F,kappa,ts);            
    scn{2}{4} = @(x,m) augs1d(x,m);                        
    scn{2}{5} = @(x) reds1d(x);                            
    scn{2}{6} = @(x) transf1d(x,@(y)FWT_PO(y,L,qmf));      
    scn{2}{7} = @(x) transf1d(x,@(y)IWT_PO(y,L,qmf));      
    scn{2}{8} = @(x,d,r) enkf_winov_diag(x,d,r);           
    scn{2}{9} = {@(x,o,a) enkf_update_diag(x,o,a),...
                 @(x,o,a) enkf_update_diag(x,o,a)};      
    % Sample covariance - full state state observed
    scn{3} = cell(9);
    scn{3}{1} = ones(n,1);                                         
    scn{3}{2} = r;                                         
    scn{3}{3} = @(u) lorenz96(u,dt,F,kappa,ts);            
    scn{3}{4} = @(x,m) x;                        
    scn{3}{5} = @(x) x;                            
    scn{3}{6} = @(x) x;      
    scn{3}{7} = @(x) x;      
    scn{3}{8} = @(x,d,r) enkf_winov_sample(x,eye(n),d,r);           
    scn{3}{9} = {@(x,o,a) enkf_update_sample(x,o,a),...
                 @(x,o,a) enkf_update_sample(x,o,a)};        

    BIAS = zeros(length(scn),reps,nac+1);
    MAE = zeros(length(scn),reps,nac+1);
    RMSE = zeros(length(scn),reps,nac+1);

    for rep_ind = 1:reps
        % initializationa
        U = init_lorenz96(n,dt,F,kappa);
        Xi = zeros(n,1,N);
        Xi(:,1,:) = repmat(U,1,N) + randn(n,N)*sqrt(r);
        Yi = zeros(n,1);
        Yi(:,1) = init_lorenz96(n,dt,F,kappa);
        fprintf('%g/%g',rep_ind,reps);
        for scn_ind = 1:length(scn)
            [X,Y] = enkf_sim_1d(Xi,Yi,nac,scn{scn_ind}); 
            [~,bias,mae,rmse]=enkf_1d_stat(X,Y);
            BIAS(scn_ind,rep_ind,:) = bias;
            MAE(scn_ind,rep_ind,:) = mae;
            RMSE(scn_ind,rep_ind,:) = rmse;
            fprintf('.');
        end
        fprintf('\n');
    end
    
    %make pictures
    
    clear bias mae rmse
    
    bias = squeeze(mean(BIAS,2));
    mae = squeeze(mean(MAE,2));
    rmse = squeeze(mean(RMSE,2));
    
    x = repmat((0:nac)',1,length(scn));
    
    subplot(131);
    plot(x,bias')
    legend(scn_names);
    xlabel('Assimilation step');
    title('Bias');
    
    subplot(132);
    plot(x,mae')
    legend(scn_names);
    xlabel('Assimilation step');
    title('MAE');
    
    
    subplot(133);
    plot(x,rmse')
    legend(scn_names);
    xlabel('Assimilation step');
    title('RMSE');
    
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
        'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    main_title = sprintf('Lorenz96 - n %g - N %g - number of observed points %g \nEnKF is done using observation of full state!!!',...
                        n,N,no);
    text(0.5, 1,main_title,...
        'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
    
    file_out = sprintf('lorenz96_n%g_N%g_no%g',n,N,no);

    [~,~] = mkdir('img');
    %savefig(['img/' file_out '.fig']);
    
    print('-dpng', ['img/' file_out '.png']);
    
    % save results for future analysis
    [~,~]=mkdir('mat');
    file_out_mat = sprintf('mat/lorenz96_n%g_N%g_no%g.mat',n,N,no); 
    save(file_out_mat,'n','N','no','BIAS','RMSE','MAE');
end
      




