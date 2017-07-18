x = 1:1:10;
v =ones(size(x));
v(5:end)=2
xq = 1:0.11:10;

% vq = interpn(x,v,xq,'pchip');
vq = interpn(x,v,xq,'linear');
vq2 = interpn(x,v,xq,'spline');

plot(x,v,'or',xq,vq,'-x',xq,vq2,'-c*');
legend('Data Points','Linear Interploation','Spline Interploation');