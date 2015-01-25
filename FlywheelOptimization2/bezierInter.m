function [ b] = bezierInter(points,t)
% Practice making a 3 point bezier curve
% points = [0,0;
%             1,1;
%             -1,3]

% points = [0,0,0;
%     1,1,1;
%     -1,3,5;
%     5,10,-1]

numControlPoints = size(points,1); % get the number of control points
dimensions = size(points, 2); % Get the dimension of the points array
b = zeros(dimensions,size(t,2)); % B will hold the polyomial associated with each point in each dimension
n = numControlPoints-1;

for dim = 1:dimensions
    for i = 0:numControlPoints-1
        % Make the x values.        
       
        % can'use the 0 index, so use the i+1 location to store the values.
        biNomial = nchoosek(n,i);
        b_temp = biNomial*(1-t).^(n-i).*t.^i*points(i+1,dim);
        b(dim,:) =  b(dim,:)+ b_temp;
        %b0x = 1*(1-t).^2.*t.^0*points(1,1);       
    end
end





%  for i = 1:numPoints
%         % Make the x values. 
%         b0x = 1*(1-t).^2.*t.^0*points(1,1);
%         b1x = 2*(1-t).^1.*t.^1*points(2,1);
%         b2x = 1*(1-t).^0.*t.^3*points(3,1);
% 
%         % Make the x\y values. 
%         b0y = 1*(1-t).^2.*t.^0*points(1,2);
%         b1y = 2*(1-t).^1.*t.^1*points(2,2);
%         b2y = 1*(1-t).^0.*t.^3*points(3,2);
%     end
% x = b0x+b1x+b2x;
% y = b0y+b1y+b2y;

