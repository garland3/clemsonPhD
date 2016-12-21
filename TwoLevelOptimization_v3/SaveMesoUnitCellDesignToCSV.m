function []= SaveMesoUnitCellDesignToCSV(designVarsMeso,macroElemProps,folderNum,macro_meso_iteration,elementNumber,newDesign)

if(newDesign ==1)
    % save the density field
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
    x = designVarsMeso.x;
    csvwrite(outname,x);

    % save the sensitivity field
%      outname = sprintf('./out%i/sensitivity%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
%      csvwrite(outname,designVarsMeso.temp1);
    
    D = macroElemProps.D_homog;
else
     D = macroElemProps.D_given;
end

 outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,elementNumber);
 csvwrite(outname,D);