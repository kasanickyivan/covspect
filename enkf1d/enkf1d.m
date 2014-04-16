
%%
%   Standart EnKF. 
%
%   Implemenation based on J. Mandel: A  Brief Tutorial on 
%   the Ensemble Kalman Filter 
%   http://arxiv.org/abs/0901.3725
%
%
%   function XA=enkf1d(XF,H,d,r,f_cov_est)
%
%   XF - forecast ensemble
%   H - observation operator
%   d - observed data
%   r - variance of observed data
%   f_cov_est - function for estimating forecast covariance - function
%   handle
%%

function XA=enkf1d(XF,H,d,r,f_cov_est)
    %   n - size of state
    %   N - number of enseble members
    N = size(XF,2);
    %   observed data length
    d_l = length(d);
    R = diag(r);
    %   input data perturbation
    D = repmat(d,1,N) + diag(sqrt(r))*randn(d_l,N);
    
    % forecast covarince
    QF = f_cov_est(XF);
    
    % Kalman gain matrix
    M = (H * QF * H' + R);
    K = QF*H/M;
    
    % Analysis
    XA=XF+K*(D-H*XF);
end