% flywheel optimizatoin using fmincon


clear; clc; close all;
% x(1) = radius
% x(2) = height

numberInputs = 10;
% x = numberInputs/2 control heights + numberInputs/2 vol fract heights


A = zeros(numberInputs,numberInputs)
b = zeros(numberInputs,1)
Aeq = []
beq = []
lb = zeros(numberInputs,1)
ub = ones(numberInputs,1)*1

x0 = 0.5*ones(1,numberInputs)

options =   optimoptions('fmincon','MaxFunEvals',1000);
[x,fval] = fmincon(@callCaculateEnergy,x0,A,b,Aeq,beq,lb,ub,@nonlinearConstraints, options);
x

% Result 1
%   0.7569    0.1432    0.4097    0.4089    0.4268    0.4668    0.3171    0.4557    0.4593    0.4733