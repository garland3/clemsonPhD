function [designVarsMeso meso_settings] = GenerateDesignVarsForMesoProblem(settings,e,macro_meso_iteration)

% --------------------------------------
% %% meso_settings
% --------------------------------------------
meso_settings = Configuration;
meso_settings.elasticMaterialInterpMethod = 2; % Hashin–Shtrikam law (average of upper and lower boundary)
meso_settings.heatMaterialInterpMethod = 5; % Hashin–Shtrikam law (average of upper and lower boundary)

% target volumes of material 1 and 2
meso_settings.v1 = 0.2; 
meso_settings.v2 = 0.2;


    
  
meso_settings.w1 = 1; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
meso_settings.iterationNum = 0;
meso_settings.doSaveDesignVarsToCSVFile = 0;
meso_settings.doPlotFinal = 0;
meso_settings.terminationCriteria =0.1; % 10%

% if meso structure designing, then make a smaller initial mesh

meso_settings.nelx = 10;
meso_settings.nely =10;
    
meso_settings= meso_settings.UpdateVolTargetsAndObjectiveWeights();
%meso_settings


% ---------------------------------
% Initialization of varriables
% ---------------------------------
designVars = DesignVars(meso_settings);

% Reuse the existing X matrix if it exists. 
if(macro_meso_iteration>1)
    oldDesignNumber = macro_meso_iteration-1
    folderNum = settings.iterationNum;
     outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,oldDesignNumber,e);
     if exist(outname, 'file') == 2
        designVars.x = csvread(outname);
     else
          designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = meso_settings.totalVolume; % artificial density of the elements
     end
    
else
    designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = meso_settings.totalVolume; % artificial density of the elements
end
designVars.w(1:meso_settings.nely,1:meso_settings.nelx)  = 1; % actual volume fraction composition of each element
fractionCurrent_V1Local =1;

designVars.temp1(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.temp2(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.complianceSensitivity(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.totalStress(1:meso_settings.nely,1:meso_settings.nelx) = 0;

% designVars.g1elastic(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.g1heat(1:meso_settings.nely,1:meso_settings.nelx) = 0;

 designVars = designVars.CalcIENmatrix(meso_settings);
 designVars =  designVars.CalcElementLocation(meso_settings);
designVars = designVars.PreCalculateXYmapToNodeNumber(meso_settings);



% if doing meso optimization, setup optimization configurations
%if ( meso_settings.mode == 4) 

% designVars = designVars.CalcElementNodeMapmatrixWithPeriodicXandY(meso_settings);
% designVars =  designVars.CalcNodeLocationMeso(meso_settings);
%end

designVarsMeso=designVars;

