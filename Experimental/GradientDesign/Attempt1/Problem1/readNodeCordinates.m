function M = readNodeCordinates(refineMesh)
% Read the node globalposition file
%  globalPosition(count,:) = [xloc yloc];
if(refineMesh ==1)
    M = csvread('nodeNumbersAndLocationsV2.csv');
else
    M = csvread('nodeNumbersAndLocations.csv');
end
% Remove the first column which is the node number. I just use the index as
% the node number
M = M(:,2:3);


