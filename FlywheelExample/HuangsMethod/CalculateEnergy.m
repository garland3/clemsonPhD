function energy = CalculateEnergy(hpoints, VolFractpoints,w,innerR, outerR)

energy = integral(@(t) energyDifferential(t,hpoints,VolFractpoints,w,innerR, outerR),0,1);


disp(energy);