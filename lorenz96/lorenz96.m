function out=lorenz96(u,dt,J,F,kappa)
% solving Lorenz 96 equation with J variables
% dx_j/dt = (1/kappa)*(x_{j-1}(x_{j+1}-x_{j-2})-x_j+F)
%
% in:
% u     : X_{t}
% dt    : time step
% J     : length of the state
% F     : forcing term of the model
% kappa 
% The trajectories are compute using the 
% fourth-order Runge-Kutta method
% it gives the iteration after a time step dt
% out    - updated solution X_{t+dt}

aux     = rk4(u,dt,F,J,kappa);
out   = u(:) + aux;

function deltay = rk4(Xold,dt,F,J,kappa)
Xold = Xold(:);
 k1 = f(Xold,F,J,kappa);
 k2 = f(Xold+1/2.*dt.*k1,F,J,kappa);
 k3 = f(Xold+1/2.*dt.*k2,F,J,kappa);
 k4 = f(Xold+dt.*k3,F,J,kappa);
 deltay= 1/6*dt*(k1+2*k2+2*k3+k4);
 
 
function k = f(X,F,J,kappa)
k=zeros(J,1);
%first the 3 problematic cases: 1,2,J
k(1)=(X(2)-X(J-1))*X(J)-X(1);
k(2)=(X(3)-X(J))*X(1)-X(2);
k(J)=(X(1)-X(J-2))*X(J-1)-X(J);
%then the general case
for j=3:J-1
 k(j)=(X(j+1)-X(j-2)).*X(j-1)-X(j);
end
%add the F    
k= (k+F)/kappa;