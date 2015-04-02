
[xpoints,ypoints] = meshgrid(0:1:10,0:1:10);
zpoints = -(xpoints-5).^2+-(ypoints-5).^2;






stepX = 0.1;
stepY = 0.1;
splineXX = 0:stepX:10; % columns

splineYY = 0:stepY:10; % rows

[splineXX_v2,splineYY_v2]=meshgrid(splineXX,splineYY);
splineZZ_v2 =  interp2(xpoints,ypoints,zpoints,splineXX_v2,splineYY_v2,'cubic'); % spline value at each x column

minz = min(min(splineZZ_v2));
maxZ = max(max(splineZZ_v2));
diffZ = maxZ - minz;

%splineZZ_v2 = splineZZ_v2+abs(minz)-10; % make the min -1
scale = 2/diffZ;
splineZZ_v2 = scale*splineZZ_v2+1; % make it range from -1 to 1

sz = size(splineYY_v2);
planeZ = zeros(sz);
szColor = ones(sz)


colormap('default')

doublesurf(splineZZ_v2,winter(64),planeZ,copper(64))

%surf(splineXX_v2,splineYY_v2,splineZZ_v2,'edgecolor','none')
%hold on


%surf(splineXX_v2,splineYY_v2,planeZ,'edgecolor','none')
%hold off