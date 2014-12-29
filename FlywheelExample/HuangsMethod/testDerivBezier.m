points = [0,0;
             1,3;
             2,0.2;
              3,0;
             4,0.5]
         
         
t =0:0.01:1;

deriv = DerivBezierInter(points,t)
b =  bezierInter(points,t)


plot(b(1,:),b(2,:),'--r')
hold all
plot(b(1,:),deriv(2,:),'--g')

x =points(:,1);
y = points(:,2);
scatter(x,y,150,'bd','fill');
hold off