function vol = massDifferential(t,hpoints,VolFractpoints, innerR, outerR)

% Coverts the t value to an r value. 
dr = outerR-innerR;
r = dr*t+innerR;
 
z= bezierInter(hpoints,t);
z = z(2,:);

volFrac =  bezierInter(VolFractpoints,t);
volFrac = volFrac(2,:);
density = calculateDensity(volFrac);

vol =  2*pi*density.*z.*r.*dr;


