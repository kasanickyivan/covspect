% Experiments with 'true' covariances in SWE model
%
%   definitions, see waterwave2 and generate_waterwave for details:
n=64;
reps=500;
init_h =10000;
dw_min=40;dw_max=40;
dh_min=1;dh_max=1;
bdd=10;
ts=1;
init_ts=10000;
init_d=1;
per_d = 1e6;

runs=1;
Y = zeros(n,n,3,reps*runs);
for run_ind = 1:runs
    i=(1:reps)+((run_ind-1)*reps);
    Y(:,:,:,i)=generate_waterwave(n,reps,init_h,dw_min,dw_max,dh_min, ...
        dh_max,bdd,ts,init_ts,init_d,per_d);
end

waterwave_anim(Y,60,100);
% adding derivatives
Y=ww_derivatives(Y);
%
%
% to avoid effect on borders wtih boundary conditions and to be able to
% plot covariances, only centerd subdomain 1x6x16 will be examined
n=8;
l = (size(Y,1) - n)/2;
i = l+1:l+n;
Z = Y(i,i,:,:);  
nvar = size(Y,4);
var_index = cell(nvar,1);
for var_ind = 1:nvar
    var_index{var_ind} = (1:n*n)+((var_ind-1)*n*n);
end
X = pack_state2d(Y);
% centering 
    
subplot(131);
plot(mean(X(var_index{1},:),1));
%ylim([9 11]);
xlabel('realization');ylabel('water level');title('mean water level');
subplot(132);
plot(mean(X(var_index{2},:),1));
xlabel('realization');ylabel('momentum');title('mean X momentum');
%ylim([-1 1]);
subplot(133);
plot(mean(X(var_index{3},:),1));
%ylim([-1 1]);
xlabel('realization');ylabel('momentum');title('mean Y momentum');
print('img/mean_stat_ww.png','-dpng')

[cov_s,loc_cov_s,cor_s,loc_cor_s] = cov_cor_2d(Z);


%Spatial space

figure();
subplot(121);
surface(squeeze(cov_s(1,1,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariance of water level.');
colorbar();
subplot(122);
surf(squeeze(cov_s(1,1,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariation of water level.');
print('img/cov_height_spatial.png','-dpng')

subplot(121);
surface(squeeze(cov_s(2,2,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariance of X momentum.');
colorbar();
subplot(122);
surf(squeeze(cov_s(2,2,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariance of X momentum.');
print('img/cov_Xmom_spatial.png','-dpng')

figure()
surf(squeeze(cor_s(1,2,:,:)));
colorbar();
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Correlation between water level and X momentum.');
print('img/cov_height_Xmom_spatial.png','-dpng')

figure()
surf(squeeze(cor_s(4,2,:,:)));
colorbar();
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Correlation water level derivation in x direction and X momentum.');
print('img/cor_Xder_Xmom_spatial.png','-dpng')

figure()
surf(squeeze(cor_s(5,2,:,:)));
colorbar();
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Correlation water level derivation in y direction and X momentum.');
print('img/cor_Yder_Xmom_spatial.png','-dpng')

% Fourier space
ZF = transf2d(Z(:,:,1:3,:),@(x) dct_iv(x));
[cov_f,~,cor_f,~] = cov_cor_2d(ZF);

surf(squeeze(cov_f(1,1,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariation of water level in Fourier space.');
print('img/cov_height_fourier.png','-dpng')

surf(squeeze(cov_f(2,2,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariance of X momentum in Fourier space.');
print('img/cov_Xmom_fourier.png','-dpng')

surf(squeeze(cor_f(1,2,:,:)));
colorbar();
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Correlation between water level and X momentum in Fourier space.');
print('img/cov_height_Xmom_fourier.png','-dpng')

% Wavelet space
% wavelet space
qmf=MakeONFilter('Coiflet',2);
L=2;
ZW = transf2d(Z, @(x) FWT_PO(x,L,qmf));
[cov_w,~,cor_w,~] = cov_cor_2d(ZW);

surf(squeeze(cov_w(1,1,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariation of water level in wavelet space.');
print('img/cov_height_wavelet.png','-dpng')

surf(squeeze(cov_w(2,2,:,:)));
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Covariance of X momentum in wavelet space.');
print('img/cov_Xmom_wavelet.png','-dpng')

surf(squeeze(cor_w(1,2,:,:)));
colorbar();
xlim([1 n*n]);ylim([1 n*n]);
xlabel('grid point index');
ylabel('grid point index');
title('Correlation between water level and X momentum in wavelet space.');
print('img/cor_height_Xmom_wavelet.png','-dpng')


