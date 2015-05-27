% min surface area of a "bottle" (extruded circle)
% s.t. volume is less than 10 units

clear; clc; close all;
% x(1) = radius
% x(2) = height

A = zeros(2,2);
b = zeros(2,1);
Aeq = A;
beq = b;
lb = ones(2,1)*-1000;
ub = ones(2,1)*1000;

x0 = [10,10]

options =   optimoptions('fmincon','MaxFunEvals',1000)
[x,fval] = fmincon(@surfacearea,x0,A,b,Aeq,beq,lb,ub,@volumeconstraint, options);