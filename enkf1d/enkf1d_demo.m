%%
%   Demo - 1D EnKF asssimilation with diagonal approximation of forecast covarince in spectral
%   space.  
%
%   Model: lorenz 96
%
%
%   
%
%
%
%%
n=128;                   % state dimension
N=4;                     % ensemble size
run_length=65;           
F=fft_matrix(n);         %transformation matrix
M=eye(n);                %mask matrix for covariance estimation in spectral space

[XA,YA] = enkf1d_lorenz96(128,4,65,@(X) cov_xform(X,F,M));
show_enkf(XA,YA);
