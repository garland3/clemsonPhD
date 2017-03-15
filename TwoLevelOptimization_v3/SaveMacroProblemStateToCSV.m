function []= SaveMacroProblemStateToCSV(config,designVars,matProp)

folderNum = config.iterationNum;
mm_iteration = config.macro_meso_iteration;

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
outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.x);

% save the Exx field
outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.Exx);

% save the Eyy field
outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.Eyy);

% save the Theta field
outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.t);


% save the volume fraction field (only if using the old methodd)
if(config.useExxEyy~=1)
    outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,mm_iteration);
    csvwrite(outname,designVars.w);
end

% save the lambda value
outname = sprintf('./out%i/lambda%i.csv',folderNum,mm_iteration);
csvwrite(outname,designVars.lambda1);   

% for the new meso design method using the Exx and Eyy and the consistency
% constraints I need
% 1. Topology var. Tells me if I need to make a new design or not. DONE
% 2. XY position map, DONE
% 3. D_system matrix for each element in the macro design. 
% 4. Displacement field. So that I can calculate the strain on each. DONE
% element. 
% 5. The target density of the meso problem. 





ne = config.nelx*config.nely;
targetMesoDensity = zeros(ne,1);


% Write the D constituitive matrix to .csv files.
% Calculuate the target density 
for e = 1:ne
    outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    %  D_flat =  matProp.SavedDmatrix(e,:);
    
    
    results = designVars.elementXYposition(e,:);
    yPosition = results(1);
    xPosition = results(2);
    
    
   x=   designVars.x(yPosition,xPosition);
%      x=   1;
   Exx=   designVars.Exx(yPosition,xPosition);
   Eyy=   designVars.Eyy(yPosition,xPosition);
   w=   designVars.w(yPosition,xPosition);
   rotation=   designVars.t(yPosition,xPosition);
    
   D_out= matProp.getDmatMatrixTopExxYyyRotVars(config,x,Exx, Eyy,rotation,w);
    
    csvwrite(outname,D_out);
    
    targetMesoDensity(e) = matProp.CalculateDensityTargetforMeso(w,x,Exx,Eyy,config);
end

%write the target density
outname = sprintf('./out%i/TargetMesoDensities%i.csv',folderNum,mm_iteration);
csvwrite(outname,targetMesoDensity);

