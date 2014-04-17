
%% Script to compare ENKF performance on Lorenz96 model
% with different strategies to encode covariance
%

n = 64;
N = 4;
steps = 32;
reps = 32;

cov_diag = { @(XF) cov(XF'), ...
             @(XF) cov_xform(XF,dst_matrix(n), eye(n)), ...
             @(XF) cov_xform(XF,dct_matrix(n), eye(n)), ...
             @(XF) cov_xform(XF,fft_matrix(n), eye(n)), ...
             @(XF) cov_xform(XF,wav_matrix(n,4,'Coiflet',2), eye(n)), ...
             @(XF) cov_xform(XF,wav_matrix(n,4,'Coiflet',4), eye(n)), ...
             @(XF) cov_xform(XF,wav_matrix(n,4,'Beylkin',2), eye(n)) };

cov_names = { 'Sample', 'DST', 'DCT', 'FFT', 'Coi-2', 'Coi-4', 'Bey' };

cov_diag = cov_diag([1 3 4 6]);
cov_names = cov_names([1 3 4 6]);

RMSE = zeros(steps, length(cov_diag),reps);

for i=1:length(cov_diag)
    fprintf('Run %d/%d ', i, length(cov_diag));
    for j=1:reps
        fprintf('.');
        f_cov = cov_diag{i};
        [XA,YA] = enkf1d_lorenz96(n,N,steps,f_cov);
        RMSE(:,i,j) = sum((XA - YA).^2).^0.5;
    end
    fprintf('\n')
end

figure('name', 'ENKF convergence vs. cov encoding');
errorbar(repmat((1:steps)',1,length(cov_diag)), mean(RMSE,3), std(RMSE,[],3)/32^0.5, 'linewidth', 1.5);
title('ENKF convergence vs. cov encoding', 'fontsize', 16);
legend(cov_names);
xlabel('# assimilation steps');