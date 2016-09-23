function  macroElementProps = GetMacroElementPropertiesFromCSV(settings,e,macro_meso_iteration)

macroElementProps = macroElementProp;
macroElementProps.elementNumber = e;

folderNum = settings.iterationNum;

% GEt element->node mapping 
outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,macro_meso_iteration);IEN = csvread(outname);

% % Get displacement field
outname = sprintf('./out%i/displacement%i.csv',folderNum,macro_meso_iteration);

U =  csvread(outname);  
nodes1=  IEN(e,:);
macroElementProps.elementNodes=nodes1;
xNodes = nodes1*2-1;
yNodes = nodes1*2;
dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
u_local =   U(dofNumbers);
offsetX = mean(u_local([1 3 5 7]));
offsetY = mean(u_local([2 4 6 8]));
% offsetX = u_local(7);
% offsetY = u_local(8);
u_local([1 3 5 7]) = u_local([1 3 5 7])-offsetX;
u_local([2 4 6 8]) = u_local([2 4 6 8])-offsetY;
macroElementProps.disp  = u_local;
% macroElementProps.disp  = transpose(macroElementProps.disp );


% Get element to XY position map (needed for x and w vars retrival)
%
%
% NodeToXYArrayMap is the node to XY map, not the element.
%
% outname =
% sprintf('./out%i/NodeToXYArrayMap%i.csv',folderNum,macro_meso_iteration);

% NodeToXYArrayMap = csvread(outname);
% [result] =  NodeToXYArrayMap(macroElementProps.elementNumber,:);

macroElementProps.elementNumber

% Save element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
elementXYposition=csvread(outname);

results = elementXYposition(macroElementProps.elementNumber,:);

macroElementProps.yPosition = results(1);
macroElementProps.xPosition = results(2);

% Get the density field
outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
x = csvread(outname);
macroElementProps.density = x(macroElementProps.yPosition,macroElementProps.xPosition );


% % Get the volume fraction field
outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
w = csvread(outname);
macroElementProps.material1Fraction = w(macroElementProps.yPosition,macroElementProps.xPosition );

% if(macro_meso_iteration>1)
outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
D = csvread(outname);
macroElementProps.D_given =D;
% end