function  macroElementProps = GetMacroElementPropertiesFromCSV(settings,e,macro_meso_iteration)

macroElementProps = macroElementProp;
macroElementProps.elementNumber = e;

folderNum = settings.iterationNum;

% GEt element->node mapping 
outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,macro_meso_iteration);
IEN = csvread(outname);

% % Get displacement field
outname = sprintf('./out%i/displacement%i.csv',folderNum,macro_meso_iteration);
U =  csvread(outname);  
nodes1=  IEN(e,:);
macroElementProps.elementNodes=nodes1;
xNodes = nodes1*2-1;
yNodes = nodes1*2;
dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
u_local =   U(dofNumbers);
macroElementProps.disp  = u_local;
macroElementProps.disp  = transpose(macroElementProps.disp );


% Get element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/NodeToXYArrayMap%i.csv',folderNum,macro_meso_iteration);
NodeToXYArrayMap = csvread(outname);
[result] =  NodeToXYArrayMap(macroElementProps.elementNumber,:);
macroElementProps.yPosition = result(1);
macroElementProps.xPosition = result(2);

% Get the density field
outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
x = csvread(outname);
macroElementProps.density = x(macroElementProps.yPosition,macroElementProps.xPosition );


% % Get the volume fraction field
outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
w = csvread(outname);
macroElementProps.material1Fraction = w(macroElementProps.yPosition,macroElementProps.xPosition );