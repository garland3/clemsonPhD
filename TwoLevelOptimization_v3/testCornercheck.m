x = rand(35,35)
% d = ones(1,10)
% x = diag(d)
% x=x+flip(x)
x(3,:)=1;

% y=csvread('./completeStucture1.000000_macroIteration_1.csv');
% % sizeToTest=500;
% % startLocation=200;
% % endLocation=startLocation+sizeToTest;
% % x=y(startLocation:endLocation,startLocation:endLocation); % make it smalller
% x=y;
 csvwrite('test.csv',x)

[t, numChanged] = CheckForConerElements(x, size(x,2), size(x,1), 0.3);