% rho = 0:0.01:1;
% 
% rhoSqred = rho.^(1/3);
% k = rho.^3;
% 
% 
% plot(rho, rhoSqred,rho,k)
% legend('rho^2', 'k')

Dtarget=[ 20648       17593       17194;
       17593       96639       17194;
       17194       17194       20525]
   
   e1=[0:0.01:1];
   e2 = 0.3;
   e3 = 0.3;
   strain=[e1;e2;e3];
   
   strainEnergy = transpose(strain)*Dtarget*strain