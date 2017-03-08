% rho = 0:0.01:1;
% 
% rhoSqred = rho.^(1/3);
% k = rho.^3;
% 
% 
% plot(rho, rhoSqred,rho,k)
% legend('rho^2', 'k')
clear
clc
close all
Dtarget=[ 20648       17593       17194;
       17593       96639       17194;
       17194       17194       20525];
   
%    Dtarget=[ 20648       5000       5000;
%        5000       21000       5000;
%        5000       5000       10000];
   
   
   e1=[0:0.01:1];
   [s1, s2]=size(e1);
   e2 = 0.3*ones(s1,s2);
   e3 = 0.3*ones(s1,s2);
   strain=([e1; e2; e3]);
   
   for i = 1:s2
            strainLocal = strain(:,i);
            strainEnergy(i) = transpose(strainLocal)*Dtarget*strainLocal;
   end
   
    plot(e1,strainEnergy);
   
   
   fun = @(x)-[x(1) x(2) x(3)]*Dtarget*[x(1);x(2);x(3)];
% A = [];
% b = [];
A = [1 1 1];
%     0 0 0;
%     0 0 0];
b = [1];
%         0 ;
%         0];
Aeq=[];
beq=[];
   
lb = [0 0 0];
ub = [1 1 1];
%lb=-ub;

x0 = ones(1,3)/3;
%x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
%x = fmincon(fun,x0,A,b)
x = fmincon(fun,x0,A,b,[],[],lb,ub)
x(1)+x(2)+x(3)

%constraint = @(x)(x(1)+x(2)+x(3)-1
x2 = ga(fun,3,A,b,[],[],lb,ub)