function IEN = readElementCordinateMap(refineMesh)
% Read the map betweeen elements and nodes
%    IEN(count,:)=[rowMultiplier*elementsInRow+j,
%                      rowMultiplier*elementsInRow+j+1,
%                      (rowMultiplier +1)*elementsInRow+j+1,
%                       (rowMultiplier +1)*elementsInRow+j];
if(refineMesh ==1)
    IEN = importdata('elementNodeMapV2.csv');
else
      IEN = importdata('elementNodeMap.csv');
end
% Remove the first column which is the element number. Instead we use the
% index as the identifier for the element number
IEN = IEN(:,[2:5]);



