function mass = CalculateMass(hpoints, VolFractpoints,innerR, outerR)

mass = integral(@(t) massDifferential(t,hpoints,VolFractpoints,innerR, outerR),0,1);


