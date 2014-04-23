% function counts index gridpoint in packed vector of 2d array
%
%
% input
%   - x - x coordinate
%   - y - y cootdinate
%   - ydim - y dimension of the original array
%
% output
%   same type as x,y - index of (x,y) value in packed vector    

function i=coor_grid2pack(x,y,ydim)
    i=(y-1)*ydim+x;
end