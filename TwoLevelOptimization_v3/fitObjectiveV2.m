function diffSquareSum= fitObjectiveV2(x,X,Y,Z,R,o,config,matProp)
%x is the design vars. 10 by 1
%  = p_000  + p_100 E_xx  + p_010 E_yy+p_001 ?+p_200 E_xx^2  + p_020 E_yy^2+p_002 ?^2+p_110 E_xx E_yy  + p_011 E_yy ?+p_101 E_xx ?

% X is data points
% Y is data poitns
% Z is the theta values
% R is the rho data

% p_000 =x(1);
% p_100=x(2);
% p_010=x(3) ;
% p_001=x(4);
% p_200=x(5);
% p_020=x(6);
% p_002=x(7);
% p_110=x(8);
% p_011=x(9);
% p_101=x(10);

%    scaleUp = matProp.E_material1;

E_xx=X;%/scaleUp;
E_yy=Y;%/scaleUp;
theta=Z;
 

Coefficents=x;

    [~, ~,rhoExperimental] = o.CalculateDensitySensitivityandRho(E_xx,E_yy,theta,Coefficents,config,matProp);

% rhoExperimental= p_000  + p_100* E_xx  + p_010* E_yy+p_001 *theta+p_200 *E_xx.^2  + p_020* E_yy.^2+p_002*theta.^2+p_110 *E_xx.* E_yy  + p_011 *E_yy .*theta+p_101* E_xx.*theta;
% rhoExperimental= x(1)  + x(2)* exp(E_xx)  + x(3)* exp(E_yy)+x(4) *exp(theta) +x(5)*E_xx  + x(6)* E_yy +x(7)*theta+ x(8)*E_xx.*E_yy; %+p_110 *E_xx.* E_yy  + p_011 *E_yy .*theta+p_101* E_xx.*theta;
%    rhoExperimental= x(1)  + x(2)* E_xx  + x(3)* E_yy+x(4) *theta+x(5) *E_xx.^2  + ...
%              x(6)* E_yy.^2+x(7)*theta.^2+x(8) *E_xx.* E_yy  + x(9) *E_yy .*theta+x(10)* E_xx.*theta;

Diff = rhoExperimental-R;
DiffSquared = Diff.^2;
diffSquareSum=sum(DiffSquared);
