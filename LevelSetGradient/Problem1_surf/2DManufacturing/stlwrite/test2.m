tmpvol = zeros(20,20,20);       % Empty voxel volume
tmpvol(8:12,8:12,5:15) = 1;     % Turn some voxels on
tmpvol(13:14,13:15,5:10) = 1;     % Turn some voxels on
fv = isosurface(tmpvol, 0.99);  % Create the patch object
stlwrite('test2.stl',fv, 'mode','ascii')         % Save to binary .stl