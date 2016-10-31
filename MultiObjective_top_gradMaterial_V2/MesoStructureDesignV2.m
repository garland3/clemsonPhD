function [D_homog,designVarsMeso,macroElemProps]=MesoStructureDesignV2(matProp,mesoSettings,designVarsMeso,masterloop,~,macroElemProps,dcGiven)
% 	------------------ nope 1. Tile the meso design domain
% 	2. Apply the strain 
% 	3. Calcualute sensitive of every locaiton. 
% 	4. Then, sum the tiled sensitives together. 
% 		1. Cut out the top and right edge tiles, they dont' seen to acurately reflect the design at all.             % 
% 	5. Update the X var of the meso unit cell. 
% 		1. Do this 3 times, (no need to rerun the applystrain method).             % 
% 	6. Get homogenous properties

doPlot = 1; % For debugging allow plotting of some information. 

% Calcualte the strain, epsilon = B*d
% get the B matrix. 
[macroElemProps.K ,~,macroElemProps.B] = matProp.effectiveElasticKEmatrix( macroElemProps.material1Fraction, mesoSettings,[]);
macroElemProps.strain = macroElemProps.B* transpose(macroElemProps.disp); % transpose the disp to be vertical

if(doPlot ==1)
    p = plotResults;
    figure(1)
    subplot(2,2,3);
end

% give periodic boundary condition. 

designVarsMeso = designVarsMeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoSettings);
designVarsMeso =  designVarsMeso.CalcNodeLocationMeso(mesoSettings);



% --------------------------------------------
%    CALCULATE SENSITIVITY
% 
%    UPDATE THE DESGIN OF THE UNIT CELL
% --------------------------------------------
% Loop 2 times. Calculatint the sensitivity and changing the design var X
for mesoLoop = 1:20
    designVarsMeso = anthonyHomgenization(designVarsMeso, mesoSettings, matProp, macroElemProps,mesoLoop);
%     designVarsMeso = designVarsMeso.CalculateSensitiviesMesoStructure( mesoSettings, matProp, masterloop,macroElemProps,U);
     % designVarsMeso = designVarsMeso.CalculateSensitiviesMesoStructure_Tile( mesoSettings, matProp, masterloop,macroElemProps,U);
%      if(dcGiven~=1)
% %         designVarsMeso = designVarsMeso.CalculateSensitiviesMesoStructureNoPeriodic( mesoSettings, matProp, masterloop,macroElemProps,U);
%     
%         designVarsMeso.dc=-designVarsMeso.temp1;
%          designVarsMeso.dc(isnan( designVarsMeso.dc))=0.1;
%      end
%      designVarsMeso.dc = designVarsMeso.temp1;
    % FILTERING OF SENSITIVITIES
    [designVarsMeso.dc]   = check(mesoSettings.nelx,mesoSettings.nely,mesoSettings.rmin,designVarsMeso.x,designVarsMeso.dc);
    
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [designVarsMeso.x] = OC(mesoSettings.nelx,mesoSettings.nely,designVarsMeso.x,mesoSettings.totalVolume,designVarsMeso.dc, designVarsMeso, mesoSettings);
    
    if(doPlot ==1)
        figure(1)
        subplot(2,2,3)
        p.PlotArrayGeneric(designVarsMeso.x,'meso design -> topology var'); % plot the results.
          subplot(2,2,4)
        p.PlotArrayGeneric(designVarsMeso.dc,'meso design -> sensitivity'); % plot the results.
        drawnow
%         figure(3)
    end
end

% --------------------------------------------
%    CALCULATE Effective constitutive matrix of the meso structure. This is
%    need for the macro optimization
% --------------------------------------------

% give periodic boundary condition. 

% designVarsMeso = designVarsMeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoSettings);
% designVarsMeso =  designVarsMeso.CalcNodeLocationMeso(mesoSettings);

macroElemProps = designVarsMeso.GetHomogenizedProperties(mesoSettings,mesoSettings, matProp, masterloop,macroElemProps);
D_homog =  macroElemProps.D_homog;