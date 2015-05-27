function exampleBVP()
% Solves a simple boundary value problem, where the differential equation
% is dependent on a parameter 'r'. r is changed and the BVP is resolved
% each time. 


clear;clc;
figure(1)
hold on
for count = 0:10:100    
    runEquation(count)
end
 hold off

function dydx = examp_ODE(x,y,r)
dydx = [y(2);-y(2)-50+r*x];

function res = examp_BV(ya,yb)
res = [ya(1);yb(1)];

function runEquation(r)

% http://stackoverflow.com/questions/2256229/matlab-how-do-i-pass-a-parameter-to-a-function
solinit = bvpinit(linspace(0,4,5),[1 0]);
sol = bvp4c(@(x,y) examp_ODE(x,y,r),@examp_BV,solinit);

x = linspace(0,4);
y = deval(sol,x);
plot(x,y(1,:))
% t = 0:0.1:5
% [t,y] = ode45(@(x,y) examp_ODE(x,y,r),t,[0,0]);
% plot(t,y)
% hold off
