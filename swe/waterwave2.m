% WATERWAVE   2D Shallow Water Model
%
% Lax-Wendroff finite difference method.
% Reflective boundary conditions.
% Random water drops initiate gravity waves.
% Surface plot displays height colored by momentum.
% Plot title shows t = simulated time and tv = a measure of total variation.
% An exact solution to the conservation law would have constant tv.
% Lax-Wendroff produces nonphysical oscillations and increasing tv.
%
% See:
%    http://en.wikipedia.org/wiki/Shallow_water_equations
%    http://www.amath.washington.edu/~rjl/research/tsunamis
%    http://www.amath.washington.edu/~dgeorge/tsunamimodeling.html
%    http://www.amath.washington.edu/~claw/applications/shallow/www
%
%    originally written by Randy LeVeque
%
%   function [X]=waterwave2(X,dt,dx,dy,nsteps)
%   
%   X(:,:,1)    - water height
%   X(:,:,2)    - X momentum
%   X(:,:,3)    - Y momentum
%   dt          - timestep
%   dx          
%   dy
%   nsteps      number of time steps

function [X]=waterwave2(X,dt,dx,dy,nsteps)
   g = 9.8;                 % gravitational constant
   n = size(X,2);          %grid size, must by n x n
   % extend state variables for boundary conditions
   H = ones(n+2,n+2); U = zeros(n+2,n+2);  V = zeros(n+2,n+2);
   Hx = zeros(n+1,n+1); Ux = zeros(n+1,n+1); Vx = zeros(n+1,n+1);
   Hy = zeros(n+1,n+1); Uy = zeros(n+1,n+1); Vy = zeros(n+1,n+1);
    
   H(2:n+1,2:n+1)=squeeze(X(:,:,1));
   U(2:n+1,2:n+1)=squeeze(X(:,:,2));
   V(2:n+1,2:n+1)=squeeze(X(:,:,3));

   for step = 1:nsteps   
       % Reflective boundary conditions
       H(:,1) = H(:,2);U(:,1) = U(:,2); V(:,1) = -V(:,2);
       H(:,n+2) = H(:,n+1);U(:,n+2) = U(:,n+1);V(:,n+2) = -V(:,n+1);
       H(1,:) = H(2,:);U(1,:) = -U(2,:);V(1,:) = V(2,:);
       H(n+2,:) = H(n+1,:);U(n+2,:) = -U(n+1,:);V(n+2,:) = V(n+1,:);

       % First half step
   
       % x direction
       i = 1:n+1;
       j = 1:n;
   
       % height
       Hx(i,j) = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));
   
       % x momentum
       Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
                 dt/(2*dx)*((U(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                            (U(i,j+1).^2./H(i,j+1) + g/2*H(i,j+1).^2));
   
       % y momentum
       Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
                 dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
                            (U(i,j+1).*V(i,j+1)./H(i,j+1)));
       
       % y direction
       i = 1:n;
       j = 1:n+1;
   
       % height
       Hy(i,j) = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));
   
       % x momentum
       Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
                 dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - ...
                            (V(i+1,j).*U(i+1,j)./H(i+1,j)));
       % y momentum
       Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
                 dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                            (V(i+1,j).^2./H(i+1,j) + g/2*H(i+1,j).^2));
   
       % Second half step
       i = 2:n+1;
       j = 2:n+1;
   
       % height
       H(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
                         (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
       % x momentum
       U(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2) - ...
                         (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx(i-1,j-1).^2)) ...
                       - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./Hy(i-1,j)) - ...
                         (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)));
       % y momentum
       V(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)) - ...
                         (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
                       - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j) + g/2*Hy(i-1,j).^2) - ...
                         (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy(i-1,j-1).^2));
         
       %if all(all(isnan(H)))
       %error('Model unstable.');
       %end
   end
   X(:,:,1)=H(2:n+1,2:n+1);
   X(:,:,2)=U(2:n+1,2:n+1);
   X(:,:,3)=V(2:n+1,2:n+1);
 end