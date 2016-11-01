function [designVarsMeso, meso_settings] = GenerateDesignVarsForMesoProblem(settings,e)

% --------------------------------------
% %% meso_settings
% --------------------------------------------
meso_settings = settings;
meso_settings.macro_meso_iteration = settings.macro_meso_iteration;
% copy a bunch of var set earlier
meso_settings.numXElmPerDV=settings.numXElmPerDV;
meso_settings.numYElmPerDV=    settings.numYElmPerDV;
meso_settings.doPlotAppliedStrain =   settings.doPlotAppliedStrain ;
meso_settings.loadingCase =   settings.loadingCase ;
meso_settings.doUseMultiElePerDV= settings.doUseMultiElePerDV; % Tell the meso settings about the mode.

meso_settings.elasticMaterialInterpMethod = 2; % Hashin�Shtrikam law (average of upper and lower boundary)
meso_settings.heatMaterialInterpMethod = 5; % Hashin�Shtrikam law (average of upper and lower boundary)

% target volumes of material 1 and 2
meso_settings.v1 = 0.2; 
meso_settings.v2 = 0.2;


    
  
meso_settings.w1 = 1; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
meso_settings.iterationNum = 0;
meso_settings.doSaveDesignVarsToCSVFile = 0;
meso_settings.doPlotFinal = 0;
meso_settings.terminationCriteria =0.1; % 10%

% if meso structure designing, then make a smaller initial mesh

meso_settings.nelx = settings.nelxMeso;
meso_settings.nely =settings.nelyMeso;
    
meso_settings= meso_settings.UpdateVolTargetsAndObjectiveWeights();
%meso_settings


% ---------------------------------
% Initialization of varriables
% ---------------------------------
designVars = DesignVars(meso_settings);

% Reuse the existing X matrix if it exists. 
if(settings.macro_meso_iteration>1)
    oldDesignNumber = settings.macro_meso_iteration-1;
    folderNum = settings.iterationNum;
     outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,oldDesignNumber,e);
     if exist(outname, 'file') == 2
        designVars.x = csvread(outname);
     else
          designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = meso_settings.totalVolume; % artificial density of the elements
     end
    
else
    %designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = meso_settings.totalVolume; % artificial density of the elements
    designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = randi([0, meso_settings.totalVolume*100],meso_settings.nely,meso_settings.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere. 
end
designVars.w(1:meso_settings.nely,1:meso_settings.nelx)  = 1; % actual volume fraction composition of each element
% fractionCurrent_V1Local =1;

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

