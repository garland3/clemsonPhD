function vol = CalculateVolume(hpoints, innerR, outerR)


vol = integral(@(t) volumeDifferential(t,hpoints,innerR, outerR),0,1);


disp(vol);

