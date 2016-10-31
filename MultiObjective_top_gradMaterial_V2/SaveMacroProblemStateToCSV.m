function []= SaveMacroProblemStateToCSV(settings,designVars,matProp)

folderNum = settings.iterationNum;
mm_iteration = settings.macro_meso_iteration;

% Save element->node mapping 
outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.IEN);

% % Save element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.elementXYposition);




% Save displacement field
outname = sprintf('./out%i/displacement%i.csv',folderNum,mm_iteration);
uout = full(designVars.U);
csvwrite(outname,uout);

% save the density field
outname = sprintf('./out%i/densityfield%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.x);


% save the volume fraction field
outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.w);

% save the lambda value
outname = sprintf('./out%i/lambda%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.lambda1);

% Write the D constituitive matrix to .csv files. 
ne = settings.nelx*settings.nely;
for e = 1:ne
    outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,mm_iteration,e);
    D_flat =  matProp.SavedDmatrix(e,:);
    D = reshape(D_flat,3,3);
    csvwrite(outname,D);
end