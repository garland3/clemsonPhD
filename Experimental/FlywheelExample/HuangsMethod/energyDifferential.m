function energy = energyDifferential(t,hpoints,VolFractpoints,w, innerR, outerR)


% Coverts the t value to an r value. 
dr = outerR-innerR;
r = dr*t+innerR;
 
z= bezierInter(hpoints,t);
z = z(2,:);

volFrac =  bezierInter(VolFractpoints,t);
volFrac = volFrac(2,:);
density = caculateDensity(volFrac);

energy =  w^2*pi*density.*z.*(r.^3).*dr;


