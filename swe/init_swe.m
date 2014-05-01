function X = init_swe(n,ih,dw,dh,mbd,dt,dx,dy)
%   Function to create initial state to SWE model.
%
%   in:
%
%   n   : grid size
%   ih   : intial height
%   dw  : drop widht
%   dh  : drop height
%   mbd : minimal  distance between boundary and drop 
%   
%   dt,dx,dy arguments to waterwave2

    X = zeros(n,n,3);
    X(:,:,1) = ones(n,n)*ih;
    X(:,:,1) = add_rand_drop(X(:,:,1),dw,dw,dh,dh,mbd);
    X = waterwave2(X,dt,dx,dy,100);
end
