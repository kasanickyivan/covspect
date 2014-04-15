% Experiments with 'true' covariances in SWE model
%
%   definitions, see waterwave2 and generate_waterwave for details:
n=64;
reps=1000;
init_h =10;
dw_min=4;dw_max=8;
dh_min=2;dh_max=4;
ts=1;
init_ts=5000;
init_d=1;
per_d = 1e6;
Y=generate_waterwave(n,reps,init_h,dw_min,dw_max,dh_min,dh_max,ts,init_ts,init_d,per_d);
waterwave_anim(Y,0.02);
% adding derivatives
Y=ww_derivatives(Y);
%
%
% to avoid effect on borders wtih boundary conditions and to be able to
% plot covariances, only centerd subdomain 1x6x16 will be examined
n=16;
l = (size(Y,1) - n)/2;
i = l+1:l+16;
Z = Y(i,i,:,:);  
nvar = size(Y,4);
var_index = cell(nvar,1);
for var_ind = 1:nvar
    var_index{var_ind} = (1:n*n)+((var_ind-1)*n*n);
end
X = pack_state(Z);


subplot(131);
plot(mean(X(var_index{1},:),1));
xlabel('realization');ylabel('water level');title('mean water level');
subplot(132);
plot(mean(X(var_index{2},:),1));
xlabel('realization');ylabel('momentum');title('mean X momentum');
subplot(133);
plot(mean(X(var_index{3},:),1));
xlabel('realization');ylabel('momentum');title('mean Y momentum');


[cov_s,~,cor_s,~] = cov_cor_2d(Z);

surf(cov_s(1,1,:,:));

[cov_2,~,cor_2,~] = cov_cor_2d(X(5:8,5:8,:,:));



[cov_yz,loc_cov_yz,cor_yz,loc_cor_yz] = cov_cor_2d(YZ);
 
ind = 1:6:900;

qmf=MakeONFilter('Coiflet',2);
L=4;
YW = transf2d(Y, @(x) FWT_PO(x,L,qmf));
YF = transf2d(Y, @(x) dct_iv(x));

ZW = transf2d(Z, @(x) FWT_PO(x,L,qmf));
ZF = transf2d(Z, @(x) dct_iv(x));



[cov_f,loc_cov_f,cor_f,loc_cor_f] = cov_cor_2d(ZF);
[cov_w,loc_cov_w,cor_w,loc_cor_w] = cov_cor_2d(ZW);


var1=1;var2=4;x=1;y=1;
subplot(131);
surf_loc(loc_cor_s,var1,var2,x,y);
subplot(132);
surf_loc(loc_cor_f,var1,var2,x,y);
subplot(133);
surf_loc(loc_cor_w,var1,var2,x,y);

plot_ind = 1:8:900;
surf(squeeze(cor_f(1,2,plot_ind,plot_ind)));
surf(squeeze(cov_f(1,2,plot_ind,plot_ind)));


% spatial space
surf(squeeze(cov_s(1,1,plot_ind,plot_ind)));
title('cov between height and height in spatial space');
surf(squeeze(cov_s(1,1,1:128,1:128)));

surf(squeeze(cov_s(1,2,plot_ind,plot_ind)));
title('cov between height and x momentum in spatial space');

surf(squeeze(cov_s(1,2,1:128,1:128)));

surf(squeeze(cov_s(2,3,plot_ind,plot_ind)));
title('cov between x momentum and y momentum in spatial space');
surf(squeeze(cov_s(2,3,1:128,1:128)));

% furier space
surf(squeeze(cov_f(1,1,plot_ind,plot_ind)));
title('cov between height and height in furier space');
surf(squeeze(cov_f(1,1,1:128,1:128)));

surf(squeeze(cov_f(1,2,plot_ind,plot_ind)));
title('cov between height and x momentum in furier space');
surf(squeeze(cov_f(1,2,1:256,1:256)));

surf(squeeze(cov_f(2,3,plot_ind,plot_ind)));
title('cov between x momentum and y momentum in furier space');
surf(squeeze(cor_f(2,3,plot_ind,plot_ind)));
title('cor between x momentum and y momentum in furier space');
surf(squeeze(cor_f(2,3,1:128,1:128)));

% wavelet space
surf(squeeze(cor_w(1,1,plot_ind,plot_ind)));surf(squeeze(cov_w(1,1,plot_ind,plot_ind)));
title('cov between height and height in wavelet space');

surf(squeeze(cov_w(1,1,1:128,1:128)));

surf(squeeze(cov_w(1,2,plot_ind,plot_ind)));
title('cov between height and x momentum in wavelet space');
surf(squeeze(cov_w(1,2,1:128,1:128)));
surf(squeeze(cor_w(1,2,1:128,1:128)));

surf(squeeze(cov_w(2,3,plot_ind,plot_ind)));
title('cov between x momentum and y momentum in wavelet space');
surf(squeeze(cor_w(2,3,plot_ind,plot_ind)));
title('cor between x momentum and y momentum in wavelet space');
surf(squeeze(cor_w(2,3,1:128,1:128)));
