nely = 10;
nelx = 20;
structure = ones(nely, nelx); 
volFraction = structure*0; % initialize the volfraction composition
volFraction(1:3:nely, 1:3:nelx) = 1;
doPlot = 1;
alpha = 1;

  [U, g1_local_square, volFracV1, volFracV2]  =  FEALevelSet_2D(structure,volFraction,  doPlot, alpha);