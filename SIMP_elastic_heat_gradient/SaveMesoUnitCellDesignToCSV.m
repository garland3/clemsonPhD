function []= SaveMesoUnitCellDesignToCSV(designVarsMeso,macroElemProps,folderNum,macro_meso_iteration,elementNumber)

% save the density field
outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
csvwrite(outname,designVarsMeso.x);

% save the density field
outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
csvwrite(outname,designVarsMeso.x);