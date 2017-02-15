function[ designVars]= ReadXMacroFromCSV( settings,designVars)

if(settings.macro_meso_iteration>1)
    folderNum = settings.iterationNum;
    previousIterationNum = settings.macro_meso_iteration-1;
    
      % get the topology densities. 
    outname = sprintf('./out%i/densityfield%i.csv',folderNum,previousIterationNum);
    designVars.x = csvread(outname);
    
      % get the volume fraction optimization vars
    outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,previousIterationNum);
    designVars.w=csvread(outname);
    
    % get the lambda value
    outname = sprintf('./out%i/lambda%i.csv',folderNum,previousIterationNum); 
     designVars.lambda1=csvread(outname);
else
    message = 'First iteration, no x matrix to read';
end
