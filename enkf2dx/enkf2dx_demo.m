%  initialization
n=32;   %grid dimension
nvar=3; %number of variables foe EnKF
N=4;    %ensemble size
%   paramaters for shallow water equations
height=4;   
dh=5;   %drop height
dw=4;   %drop width
ds=5;   %initial droplet start point
dx=1;
dy=1;
dt=0.01;

rl=15;  %assimilation run length
%   number of initial steps
init_steps=1000;
%assimiltion period
ap = 100; 

%   Y - 4 dimensional array
%       - 1st,2nd dimension - grid
%       - 3rd dim - variables
%       - 4th dim - analysis steps
%

Y=zeros(n,n,nvar,rl);Y(:,:,1,1)=ones(n,n)*height;
Z=zeros(n,n,nvar,rl);Z(:,:,1,1)=ones(n,n)*height;

%initial drop 
Y(1+ds:dw+ds,1+ds:dw+ds,1,1)=squeeze(Y(1+ds:dw+ds,1+ds:dw+ds,1,1))+droplet(dh,dw);
%   initial drop to reference run
%   drop is moved 
Z(dw/2+1+ds:dw*1.5+ds,dw/2+1+ds:dw*1.5+ds,1,1)= ...
    squeeze(Z(dw/2+1+ds:dw*1.5+ds,dw/2+1+ds:dw*1.5+ds,1,1))+droplet(dh,dw);
%
%   first steps to initialize the state
Y(:,:,:,1)=waterwave2(squeeze(Y(:,:,:,1)),dt,dx,dy,init_steps);
Z(:,:,:,1)=waterwave2(squeeze(Z(:,:,:,1)),dt,dx,dy,init_steps);


%   Y - "true" state, which will be assimilated, series of full states
%   Z - reference run without any assimilation, which the ensemble will
%   start from
for rl_ind = 2:rl
    Y(:,:,:,rl_ind) = waterwave2(squeeze(Y(:,:,:,rl_ind-1)),dt,dx,dy,ap);
    Z(:,:,:,rl_ind) = waterwave2(squeeze(Z(:,:,:,rl_ind-1)),dt,dx,dy,ap);
end

% ensemble intitialization
ens_init = repmat(squeeze(Z(:,:,:,1)),[1 1 1 N]) + randn(n,n,nvar,N)*0.1;
% series of observation
obs = squeeze(Y(:,:,1,:));
fprintf('\nDST ');
Adst = assim2d(ens_init,obs,@(x) waterwave2(x,dt,dx,dy,ap),...
            @(x,y) enkf2dx(x,y,0.1,1,...
                @(x) transf2d(x, @(x) dst_i(x)), ...
                @(x) transf2d(x, @(x) dst_i(x))));
fprintf('\nDCT ');
Adct = assim2d(ens_init,obs,@(x) waterwave2(x,dt,dx,dy,ap),...
            @(x,y) enkf2dx(x,y,0.1,1,...
                @(x) transf2d(x, @(x) dct_iv(x)), ...
                @(x) transf2d(x, @(x) dct_iv(x))));

qmf=MakeONFilter('Coiflet',2);
L=4;
fprintf('\nCoiflets ');
Acoi = assim2d(ens_init,obs,@(x) waterwave2(x,dt,dx,dy,ap),...
            @(x,y) enkf2dx(x,y,0.1,1,...
                @(x) transf2d(x, @(x) FWT_PO(x,L,qmf)), ...
                @(x) transf2d(x, @(x) IWT_PO(x,L,qmf))));
            
%enkf2dx_anim(AF,Y,AW,1);
% 
figure('name','RMSE');
rmse=zeros(3,nvar,rl);
rmse(1,:,:) = enkf2dx_rmse(Adst,Y);
rmse(2,:,:) = enkf2dx_rmse(Adct,Y);
rmse(3,:,:) = enkf2dx_rmse(Acoi,Y);
plot(squeeze(rmse(:,2,:))');
legend('DST','DCT','Coi');

 

% Mask matrix
%M = zeros(n,n);
%M(2:5,2:5) = 1;
%M(20:25,20:25) = 1;
 
%AFsg = assim2d(ens_init,obs,@(x) waterwave2(x,dt,dx,dy,ap),@(x,y) enkf2dx_sgo(x,y,M,0.1,1,F,F));
%AWsg = assim2d(ens_init,obs,@(x) waterwave2(x,dt,dx,dy,ap),@(x,y) enkf2dx_sgo(x,y,M,0.1,1,W,W));
 
% figure('name','RMSE - subgrid observed');
% rmse_sg = enkf2dx_rmse(AFsg,Y,AWsg);
% plot(squeeze(rmse_sg(:,1,:))');
% legend('FFT','Wav');
 
 


