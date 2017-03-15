function [sumDiff]= ObjectiveCalculateEffectiveVars(x,DmatrixIN, matProp,config)
Exx=x(1);
Eyy=x(2);
Theta=x(3);
topDensity=1;
material1Fraction=[];
D_out =  matProp.getDmatMatrixTopExxYyyRotVars(config,topDensity,Exx, Eyy,Theta,material1Fraction);

diff = abs(D_out-DmatrixIN);
sumDiff = sum(sum(diff));