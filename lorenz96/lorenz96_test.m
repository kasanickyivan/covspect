function lorenz96_test
    clear all
    close all
    dt = 0.1;
    n = 100;
    m= 50;
    J = 64;
    F = 8;
    kappa = 1;
    u = zeros(J,1);
    u(20)=10;%1e-1;
    U20 = zeros(n,m);
    for j=1:m
        for i=1:n
            U20(i,j)=u(20);
            u = lorenz96(u,dt,J,F,kappa);

            %figure(1)
            %plot(u)
            %axis([0,40,-10,10])
            %drawnow
        end
%    u = 1e-1*randn(J,1);
%    u(20)=10+1e-1*randn(1);
    end
    hh = figure(1)
    %plot(U20)
    hold on
    plot(U20(:,2:end),'k')
    plot(U20(:,1),'r')
    ylabel('x(20)');
 xlabel('time steps');
% set(hh,'Fontsize',9);
 saveas(hh,'lorenz96_perturbation.fig')
 eval(['print -depsc lorenz96_perturbation.eps']);
 eval(['print -depsc lorenz96_perturbation.png']);
  
    