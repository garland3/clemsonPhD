function levelSetSlice()

% generate 4x4 inche grid
% put a z = 1 point at (1,2) and (3,2)
% make a distance function from them

[X,Y] = meshgrid(0:.2:4, 0:.2:4);
Z = sin(X)+sin(2*Y)


% normalize
maxZ = max(max(Z))
minZ = min(min(Z))
d = maxZ - minZ
Z = Z-minZ;
Z = Z/d;

surf(X,Y,Z)

% sanitiy check, ok!
maxZ = max(max(Z))
minZ = min(min(Z))


la
