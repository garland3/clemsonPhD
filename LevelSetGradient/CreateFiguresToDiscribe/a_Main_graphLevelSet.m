% ---------------------------
% Creates figures used in journal paper. 
% ---------------------------



% ---------------------------
% 1D level set, 2d phi function. 
% ---------------------------
clear
clc
close all
x = 0:0.1:10
y = 1-x.^2/60
y2 = x*0


y3 = vertcat(y,y2);
 createfigureLevelSet(x, y3); % actually does the plotting

% ---------------------------
% 2D level set, 3d phi function. 
% ---------------------------

x2 = 0:0.1:6
y2 = 0:0.1:2

[xMesh, yMesh]=meshgrid(x2,y2);
% Z = np.sin(R)-0.4+X/4+Y/4
% R = np.sqrt(X**2 + Y**2)
R = sqrt(xMesh.^2+yMesh.^2)
zMesh = sin(R)-0.4+xMesh/4+yMesh/4;

levelSet = zMesh*0
colorSet = levelSet-2

[xl,yl] = size(zMesh);

for i = 1:xl
    for j = 1:yl
        zz = zMesh(i,j);
        if(zz>1)
             zMesh(i,j)=1;
        end
    end
end


    figure(2)

h = surfc(xMesh,yMesh,zMesh)
colormap winter
set(h, 'edgecolor','none')
hold on
h = surf(xMesh,yMesh,levelSet, colorSet )

set(h, 'edgecolor','none')
hold on
[C,hh] = contour3(xMesh,yMesh,zMesh,[0,0],'k')
set(hh,'LineWidth',5)


xlabel('X','FontSize',10, 'FontName','Arial')
ylabel('Y','FontSize',10, 'FontName','Arial')

%createfigureLevelSetFunction3D(xMesh, yMesh, levelSet, zMesh, colorSet)

[az,el] = view % get the current view perspective, so that we can reuse it 

for i = 1:xl
    for j = 1:yl
        zz = zMesh(i,j);
        if(zz>1)
             zMesh(i,j)=1;
        end
        
        if(zz<0)
             zMesh(i,j)=-1;
        end
    end
end

% Set the values in in ches and give mm label
% X axis
ax = gca;
mmPerInch = 25.4;
step = 0.590551; % inches (15 mm)
inchTicks = 0:step:6;
set(ax,'XTick',inchTicks ); % 50 mm = 1.9685 inches
mmTicks = inchTicks*mmPerInch;
set(ax,'XTickLabel',mmTicks );

% Set the values in in ches and give mm label
% Y axis
ax = gca;
mmPerInch = 25.4;
step = 0.590551; % inches (15 mm)
inchTicks = 0:step:2;
set(ax,'YTick',inchTicks ); % 50 mm = 1.9685 inches
mmTicks = inchTicks*mmPerInch;
set(ax,'YTickLabel',mmTicks );


lsfProjection_v2(xMesh,yMesh,zMesh)
%contourf(xMesh,yMesh,zMesh,30)
%colormap winter

%view(az,el)

%hold off

