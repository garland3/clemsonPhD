function [x]= GetMesoUnitCellDesignFromCSV(folderNum,macro_meso_iteration,elementNumber)

% save the density field
outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
x=csvread(outname);

% % save the density field
% outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
% csvwrite(outname,designVarsMeso.x);