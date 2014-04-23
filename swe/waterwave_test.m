% grid size
d=21;
n=64;
% dt
dt = 0.01;

% initialization
H = ones(n,n);   
U = zeros(n,n); 
V = zeros(n,n);

figure(1);
clf
shg
set(gcf,'numbertitle','off','name','Shallow_water')
%subplot(1,2,1);
x = (0:n-1)/(n-1);
surfplot1 = surf(x,x,ones(n,n),zeros(n,n));
grid off
axis([0 1 0 1 -1 3])
caxis([-1 1])
shading faceted
c = (1:64)'/64;
cyan = [0*c c c];
colormap(cyan)
top = title('Click start');
start = uicontrol('position',[20 20 80 20],'style','toggle','string','start');
stop = uicontrol('position',[120 20 80 20],'style','toggle','string','stop');

% subplot(1,2,2);
% img2 = imagesc(ones(n));

H(1:d,1:d) = H(1:d,1:d) + droplet(1,d);
F = wav_matrix(n,4,'Coiflet',2);
for s = 1:100
    [H,U,V]=waterwave(H,U,V,0.01,1,1,100);
    C = abs(U) + abs(V);  % Color shows momemtum
    t = s*dt;
    tv = norm(C,'fro');
    set(surfplot1,'zdata',H,'cdata',C);
    set(top,'string',sprintf('t = %6.2f,  tv = %6.2f',t,tv))
%    W = FWT2_PO(H,4,MakeONFilter('Coiflet', 2));
%    set(img2,'cdata', F*H*F');
    drawnow
end
   