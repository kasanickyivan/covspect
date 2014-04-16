
%%
%  Basic plotting function for 1D EnKF results.
%
%
%   function show_enkf1d(XA,YA)
%
%   XA - mean analysis ensemble 
%   YA - true state
%%



function show_enkf1d(XA,YA)

    figure('name','ENKF diagnostics')
    
    subplot(221);
    imagesc(XA);
    title('Ensemble average');
    xlabel('time');
    subplot(222);
    imagesc(YA);
    xlabel('time');
    title('"True" state');
    subplot(223);
    imagesc(XA-YA);
    xlabel('time');
    title('Kf error');
    subplot(224);
    plot(sum((XA-YA).^2).^0.5);
    title('RMSE');
    xlabel('time');
end