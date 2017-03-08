function [x]= GetMesoUnitCellDesignFromCSV(settings,elementNumber)

% save the density field
folderNum= settings.iterationNum;
outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,settings.macro_meso_iteration,elementNumber);

if exist(outname, 'file') == 2
x=csvread(outname);
else
    fprintf('Error Getting Elemeent %i\n',elementNumber);
    x = zeros(settings.nelyMeso,settings.nelxMeso);
end

% % save the density field
% outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
% csvwrite(outname,designVarsMeso.x);