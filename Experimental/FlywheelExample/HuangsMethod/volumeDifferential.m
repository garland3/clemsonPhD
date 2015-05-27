function vol = volumeDifferential(t,hpoints,innerR, outerR)
% Coverts the t value to an r value. 
dr = outerR-innerR;
r = dr*t+innerR;
 
z= bezierInter(hpoints,t);
z = z(2,:);

vol =  2*pi*z.*r.*dr;