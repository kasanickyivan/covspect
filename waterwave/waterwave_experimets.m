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


[cov_s,loc_cov_s,cor_s,loc_cor_s] = cov_cor_2d(Z);

figure()
surf(squeeze(cov_s(1,1,:,:)));

