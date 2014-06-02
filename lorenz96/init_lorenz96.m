function X = init_lorenz96(n,dt,F,kappa)
%   function intialize random vector for lorenz 96 model
%
%   in:
%   n       :   lenght of lorenz vector state
%   dt,F,kappa    :  arguments to lorenz96 function
    X = rand(n,1)*0.001-0.0005;
    % initial number of time steps 
    ints = 1000;
    X = lorenz96(X,dt,F,kappa,ints);
end