% count central finite diferences for sample generate_waterwave output
%
% input: 
%   Y - 4d array (n x n x 3 x replications)
%
% output:
%   Z - 4d array (n-2 x n-2 x 5 x replications)
%       - (:,:,4,:) - central finite diferences in x direction
%       - (:,:,5,:) - central finite diferences in y direction
function Z = ww_derivatives(Y)
    [n,~,~,reps] = size(Y);
    Z = zeros(n-2,n-2,5,reps);
    Z(:,:,1:3,:) = Y(2:n-1,2:n-1,1:3,:);
    % x direction
    i=1:n-2;
    j=2:n-1;
    Z(:,:,4,:) = (Y(i+2,j,1,:) - Y(i,j,1,:))/2;
    % y direction
    i=2:n-1;
    j=1:n-2;
    Z(:,:,5,:) = (Y(i,j+2,1,:) - Y(i,j,1,:))/2;
end