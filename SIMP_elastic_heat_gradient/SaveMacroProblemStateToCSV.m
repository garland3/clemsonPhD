function []= SaveMacroProblemStateToCSV(settings,designVars,macro_meso_iteration)

folderNum = settings.iterationNum;

% Save element->node mapping 
outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,macro_meso_iteration);
csvwrite(outname,designVars.IEN);

% Save element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/NodeToXYArrayMap%i.csv',folderNum,macro_meso_iteration);
csvwrite(outname,designVars.NodeToXYArrayMap);


% Save displacement field
outname = sprintf('./out%i/displacement%i.csv',folderNum,macro_meso_iteration);
csvwrite(outname,full(designVars.U));

% save the density field
outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
csvwrite(outname,designVars.x);


% save the volume fraction field
outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
csvwrite(outname,designVars.w);