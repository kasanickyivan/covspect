n = 32;        %length of state vector
N = 4;          %number of ensebles 
ts = 1;       %time step between assismilations
r = 0.1;        %variance of the observations
M = zeros(n);
M(1:16,1:16)=1;     % mask matrix  
nac = 5;       % number of assimilation cyccles


%arguments swe function
dt=1;dx=150000;dy=150000;
%initial condition to swe
ih = 10000; %initial water height (water level)
dw = 15; %width of drop at begining
dh = 1000; % height of intitial drop
mbd = 5; % minimal boundary distance 

%argument to Coiflets 
qmf=MakeONFilter('Coiflet',2);
L=4;


scn = cell(1);
scn{1} =cell(9);
 %fourier diagonal
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

scn{3}{1} = M;                                         
scn{3}{2} = r;                                         
scn{3}{3} = @(x) waterwave2(x,dt,dx,dy,ts);            
scn{3}{4} = @(x,m) augs2d(x,m);                        
scn{3}{5} = @(x) reds2d(x);                            
scn{3}{6} = @(x) transf2d(x,@(y) FWT_PO(y,L,qmf));     
scn{3}{7} = @(x) transf2d(x,@(y) IWT_PO(y,L,qmf));     
scn{3}{8} = @(x,d,r) enkf_winov_diag(x,d,r);           
scn{3}{9} = {@(x,o,a) enkf_update_diag(x,o,a),...
             @(x,o,a) enkf_update_diag(x,o,a),...
             @(x,o,a) enkf_update_sample(x,o,a),...
             @(x,o,a) enkf_update_sample(x,o,a)};     
         
scn{4}{1} = ones(n,n);                                         
scn{4}{2} = r;                                         
scn{4}{3} = @(x) waterwave2(x,dt,dx,dy,ts);            
scn{4}{4} = @(x,m) x;                        
scn{4}{5} = @(x) x;                            
scn{4}{6} = @(x) x;     
scn{4}{7} = @(x) x;     
scn{4}{8} = @(x,d,r) enkf_winov_sample(x,sparse(eye(n*n)),d,r);           
scn{4}{9} = {@(x,o,a) enkf_update_sample(x,o,a),...
             @(x,o,a) enkf_update_sample(x,o,a),...
             @(x,o,a) enkf_update_sample(x,o,a),...
             @(x,o,a) enkf_update_sample(x,o,a)};     
         
         
         
         
U = init_swe(n,ih,dw,dh,mbd,dt,dx,dy);
Xi = zeros(n,n,3,N);
Xi(:,:,:,:) = repmat(U,[1 1 1 N]) + randn(n,n,3,N)*sqrt(r);
Yi = init_swe(n,ih,dw,dh,mbd,dt,dx,dy); 


      
[X,Y] = enkf_sim_2d(Xi,Yi,nac,scn{4}); 
XM = squeeze(mean(X,3));
Y = squeeze(Y);
rmse = mean((Y-XM).^(2),1).^(.5);
plot(rmse);
surface(XM);




