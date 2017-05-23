function []= SaveMesoUnitCellDesignToCSV(DVmeso,macroElementProperties,configMeso,newDesign)

folderNum=configMeso.iterationNum;
macro_meso_iteration=configMeso.macro_meso_iteration;
elementNumber=macroElementProperties.elementNumber;


if(configMeso.multiscaleMethodCompare~=1)
    if(newDesign ==1)
        % save the density field
        outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
        x = DVmeso.x;
        csvwrite(outname,x);
        
        % save the psuedo strain values
        outname = sprintf('./out%i/psuedostrain_Ite%i_forElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
        p =  macroElementProperties.psuedoStrain;
        csvwrite(outname,p);
        
        % save the final volume
        outname = sprintf('./out%i/volumeUsed_Ite%i_forElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
        v =    configMeso.totalVolume;
        csvwrite(outname,v);
        
        
        % save the sensitivity field
        %      outname = sprintf('./out%i/sensitivity%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
        %      csvwrite(outname,designVarsMeso.temp1);
    end
else
    % save the density field
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
    x = DVmeso.x;
    csvwrite(outname,x);
    
end

   
D = macroElementProperties.D_subSys;
    

    

 outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,elementNumber);
 csvwrite(outname,D);