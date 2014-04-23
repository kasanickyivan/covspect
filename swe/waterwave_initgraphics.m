function [surfplot,top,start,stop] = waterwave_initgraphics(Y)
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,start,stop] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.
% copied from waterwave_orig.m
   [n,~,~,reps] = size(Y);
   h_max = max(reshape(Y(:,:,1,:),n*n*reps,1));
   h_min = min(reshape(Y(:,:,1,:),n*n*reps,1));
   h_mean = mean(reshape(Y(:,:,1,1),n*n,1));
   clf
   shg
   set(gcf,'numbertitle','off','name','Shallow_water')
   x = (0:n-1)/(n-1);
   surfplot = surf(x,x,ones(n,n)*h_mean,zeros(n,n));
   grid off
   axis([0 1 0 1 h_min h_max])
   caxis([-1 1])
   shading faceted
   c = (1:n)'/n;
   cyan = [0*c c c];
   colormap(cyan)
   top = title('Click start');
   start = uicontrol('position',[20 20 80 20],'style','toggle','string','start');
   stop = uicontrol('position',[120 20 80 20],'style','toggle','string','stop');
end