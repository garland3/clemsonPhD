function vol = volumeDifferential3(r,theta, z, hpoints,innerR, outerR)


t= convertToT(r,innerR,outerR)
dr = DerivBezierInter(hpoints,t)
dr = dr(2,:);
 






vol =  2*pi*z.*r.*dr;