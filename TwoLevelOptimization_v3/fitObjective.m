function diffSquareSum = fitObjective(x,X,Y,Z)
%x is the design vars. 6 by 1 array
%  =  p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2;

% X is data points
% Y is data poitns
% Z the z value at each x,y point

 p00 =x(1);
 p10=x(2);
 p01=x(3) ;
 p20=x(4);
 p11=x(5);
 p02=x(6);

Zexperimental=   p00 + p10*X + p01*Y + p20*X.^2 + p11*X.*Y + p02*Y.^2;

Diff = Zexperimental-Z;
DiffSquared = Diff.^2;
diffSquareSum=sum(DiffSquared);
