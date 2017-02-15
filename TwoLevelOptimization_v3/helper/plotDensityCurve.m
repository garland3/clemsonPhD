x = 0:0.01:1;

% off = -0.45
% y= (x-off)-1./(x-off).^2

off =0
y= (x-off)-1./sin(pi-1+x)
y(y>2)=2
y(y<0)=0
y2 = x;
% dx = 1+sqrt(x)
plot(x,y,x,y2);%,x,dx)
axis square
legend('function','linear');%,'dx')