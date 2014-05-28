
% Unpacks the state of the water wave model from one column X
% into the variables H,U,V.
%
% Y = unpack_state(X,n)
%   X - packed state
%   nx - 1st dimension of 2D spatial field
%   ny - 2nd dimension of 2D spatial field
%

function Y = unpack_state(X,nx,ny)
    [m,N] = size(X);
    m = m/nx/ny;
    Y = reshape(X,nx,ny,m,N);
end