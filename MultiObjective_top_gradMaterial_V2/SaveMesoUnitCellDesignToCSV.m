function []= SaveMesoUnitCellDesignToCSV(designVarsMeso,macroElemProps,folderNum,macro_meso_iteration,elementNumber,newDesign)

if(newDesign ==1)
    % save the density field
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
    csvwrite(outname,designVarsMeso.x);

    % save the density field
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
    csvwrite(outname,designVarsMeso.x);
    
    D = macroElemProps.D_homog;
else
     D = macroElemProps.D_given;
end

 outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,elementNumber);
 csvwrite(outname,D);