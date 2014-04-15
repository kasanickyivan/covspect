% Function add random drop (random height, random width, random location)
% to 2d array - eater level
%
% inputs:
%   Y - 2d array, water level
%   dw_min - minimal width of drop to be added
%   dw_max - maximal width of drop to be added
%   dw_min - minimal height of drop to be added
%   dw_max - maximal height of drop to be added
%   bdd - boundary drop distance - minimal distance between drop and borders 
%
%
% output:
%   Y - 2d array, water level
function Y=add_rand_drop(Y,dw_min,dw_max,dh_min,dh_max,bdd)
    [nx,ny] = size(Y);
    dw=round(rand*(dw_max-dw_min)+dw_min);
    dh=round(rand*(dh_max-dh_min)+dh_min);
    df_x=round(rand*(nx-dw-2*bdd))+1+bdd; %drop from x coordinates
    dt_x=df_x+dw-1; %drop to x coordinates 
    df_y=round(rand*(ny-dw))+1; %drop from y coordinates
    dt_y=df_y+dw-1; %drop to y coordinates 
    Y(df_x:dt_x,df_y:dt_y)=Y(df_x:dt_x,df_y:dt_y)+droplet(dh,dw);
end