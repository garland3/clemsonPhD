function res = exampBV(ya,yb)
res = [ya(1)-10;yb(1)+10]

% y(0) = 10
% y(5) = -10

%bvp4c(odefun,bcfun,solinit) integrates a system of ordinary differential equations of the form

%y? = f(x,y)

%on the interval [a,b] subject to two-point boundary value conditions bc(y(a),y(b)) = 0.

% so the ya is the value at the beginning. It needs to evaluted to 0, so
% rearrange to make an equation that is equal to zero

% yb is the same, but at the other endof the interval. 

% So for y(0) = 0, and y(5) = 0, on the interval from 0<y<5, 
% res = [ya(1), yb(1)] 

% For y(0) = 5, and y(5) -5, then 
% res = [ya(1)-5,yb(1)+5]