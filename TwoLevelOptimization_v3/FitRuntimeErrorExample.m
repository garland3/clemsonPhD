function [] = FitRuntimeErrorExample()

% 
% x = 0:0.1:5;
% x=x';
% y = x.^2;
% 
% 
% sf = fit(x,y,'Poly2');
% plot(sf,x,y)
%  xlabel('ResponseExx');
% ylabel('ResponseEyy');
% zlabel('ResponseDensity');

x = [0 1 2 3 4 5 6 7 8];
t = [0 0.84 0.91 0.14 -0.77 -0.96 -0.28 0.66 0.99];
plot(x,t,'o')


net = feedforwardnet(10);
net = configure(net,x,t);
y1 = net(x)
plot(x,t,'o',x,y1,'x')