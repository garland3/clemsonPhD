% This script is helping me learn fmin con. 

% 0 ? x1 + 2x2 + 2x3 ? 72.
clear; clc; close all;
% rewrite as 
% –x1–2x2–2x3 ? 0
% x1 + 2x2 + 2x3? 72

% A = [-1 -2 -3; ...
%     1 2 2];
% b = [0;72];
A = zeros(3,3);
b = zeros(3,1);
Aeq = A;
beq = b;
lb = ones(3,1)*-1000;
ub = ones(3,1)*1000;

x0 = [10,10,10]

options =   optimoptions('fmincon','MaxFunEvals',1000)
[x,fval] = fmincon(@objective,x0,A,b,Aeq,beq,lb,ub,@mycon, options);







