[x,y] = meshgrid(0:10,0:10);
u = x;
v = zeros(size(y));

close all

subplot(1,2,1)
hold on
quiver(x,y,u,v)

quiver(x,y,v,u)

hold off

subplot(1,2,2)

u = x;
v = zeros(size(y));
quiver(x,y,u,v)
hold on

u = zeros(size(y));
v = ones(size(y))*0.5;
quiver(x,y,u,v)

hold off