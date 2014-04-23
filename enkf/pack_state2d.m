
% Packs the state of the water wave model into one column X.
%
% X = pack_state(Y)
%

function X = pack_state(Y)
    [nx,ny,m,N] = size(Y);
    X=reshape(Y,nx*ny*m,N);
end