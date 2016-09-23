function[ designVars]= ReadXMacroFromCSV(macro_meso_iteration, settings,designVars)

if(macro_meso_iteration>1)
    folderNum = settings.iterationNum;
    previousIterationNum = macro_meso_iteration-1;
    outname = sprintf('./out%i/densityfield%i.csv',folderNum,previousIterationNum);
    designVars.x = csvread(outname);
else
    message = 'First iteration, no x matrix to read';
end
