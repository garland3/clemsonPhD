function [ bfinal] = DerivBezierInter(points,t)
% Practice making a 3 point bezier curve
% points = [0,0;
%             1,1;
%             -1,3]

%  https://www.evernote.com/shard/s116/nl/18510164/4c6865dc-8b1b-4bb2-902f-98044e4d3f39
% http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html


n = length(points);
q  = zeros(n-1,2);

if(n>1)
    for i = 1:n-1
        q(i,1)= n*(points(i+1,1)- points(i,1));
        q(i,2) =n*(points(i+1,2)- points(i,2));
    end
end

%points2 = points(n-1,:); % Remove last row to make n-1
bfinal= bezierInter(q,t);


% points = [0,0,0;
%     1,1,1;
%     -1,3,5;
%     5,10,-1]
% delta = 0.001;
% bfinal(2,1)=0;  
% if(numel(t)>1)
%     btemp = bezierInter(points,t);
%     
%     for i  = 2:numel(t)-1
%         bfinal(1,i) = btemp(1,i); % x values are the same location
%         delta = btemp(1,i+1)-btemp(1,i-1);
%         bfinal(2,i) = (btemp(2,i+1)-btemp(2,i-1))/delta;
%     end
%     
% else
%   
%     t_temp = t-delta:delta:t+delta;
%     btemp = bezierInter(points,t_temp);
%     
%         i = 2;
%         bfinal(1,i-1) = btemp(1,i); % x values are the same location
%         delta = btemp(1,i+1)-btemp(1,i-1);
%         bfinal(2,i-1) = (btemp(2,i+1)-btemp(2,i-1))/delta;
%     
% end
% 
% 
% 
% 
% function u= convertToU(x,min,max)
% % Coverts the x value between xMin and xMax, to a u value between 0 and 1
% u = (x-min)/(max-main);
% 
% 
% numControlPoints = size(points,1); % get the number of control points
% dimensions = size(points, 2); % Get the dimension of the points array
% 
% for dim = 1:dim
%     for i = 1:numControlPoints-1
%         pointsNew(i,dim) = (numControlPoints-1)*(points(i+1,dim)-points(i,dim));
%     end
% end
% 
% points = pointsNew; % updated the control poitns to be the new ones
% numControlPoints = numControlPoints -1; % We lost one control point while taking the derivative
% 
% b = zeros(dimensions,size(t,2)); % B will hold the polyomial associated with each point in each dimension
% n = numControlPoints-1;
% 
% for dim = 1:dimensions
%     for i = 0:numControlPoints-1
%         % Make the x values.        
%        
%         % can'use the 0 index, so use the i+1 location to store the values.
%         biNomial = nchoosek(n,i);
%         b_temp = biNomial*(1-t).^(n-i).*t.^i*points(i+1,dim);
%         b(dim,:) =  b(dim,:)+ b_temp;
%         %b0x = 1*(1-t).^2.*t.^0*points(1,1);       
%     end
% end
% 


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
% end
% x = b0x+b1x+b2x;
% y = b0y+b1y+b2y;

