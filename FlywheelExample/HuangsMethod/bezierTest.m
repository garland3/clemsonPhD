function bezierTest(points)
%points are the control points
% t is the 
%  points = [0,0.5;
%              1,0;
%              2,1;
%              3,0.5];
t = 0:0.01:1; % make t go from 0 to 1 in 100 steps
 %b = bezierInter(points,t);
 bd = DerivBezierInter(points,t);
 
dimensions = size(points, 2); % Get the dimension of the points array
if(dimensions ==2)
    plot(b(1,:),b(2,:),'--r')
     hold all
    % plot(bd(1,:),bd(2,:),'--g')
   
    x = transpose(points(:,1));
    y = transpose(points(:,2));
    scatter(x,y,150,'bd','fill');
    
    t= strcat('Volume Composition of Material Along the Length');
    title(t);
    ylabel('Volume Fraction Composition') % x-axis label
    xlabel('Distance along length of Object') % y-axis label
    hold off
else
    plot3(b(1,:),b(2,:),b(3,:))
    hold all
end


