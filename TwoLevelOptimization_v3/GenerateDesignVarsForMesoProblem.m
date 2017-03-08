function [DVmeso, mesoConfig] = GenerateDesignVarsForMesoProblem(mesoConfig,e,macroElementProperties)

% --------------------------------------
%  mesoConfig
% --------------------------------------------

% target volumes of material 1 and 2
mesoConfig.v1=macroElementProperties.targetDensity;
mesoConfig.v2=0;
mesoConfig.totalVolume= mesoConfig.v1+0;

mesoConfig.nelx = mesoConfig.nelxMeso;
mesoConfig.nely =mesoConfig.nelyMeso;

% ---------------------------------
% Initialization of varriables
% ---------------------------------
DVmeso = DesignVars(mesoConfig);

    % % Reuse the existing X matrix if it exists.
    % if(mesoConfig.macro_meso_iteration>1)
    %     DesignNumber = mesoConfig.macro_meso_iteration-1;
    %     folderNum = mesoConfig.iterationNum;
    %     outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,e);
    %     if exist(outname, 'file') == 2
    %         DVmeso.x = csvread(outname);
    %     else
    %         %DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = mesoConfig.totalVolume; % artificial density of the elements
    %         %  DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([0, mesoConfig.totalVolume*100],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
    %         %          DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = zeros(mesoConfig.nely,mesoConfig.nelx);       
    %     end
    % else
    DVmeso = DVmeso.PreCalculateXYmapToNodeNumber(mesoConfig);
    % DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) =randi([0, 100],mesoConfig.nely,mesoConfig.nelx)/100*mesoConfig.totalVolume; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
    
    
% ---------------------------------
% MOved initialization of the x matrix to the designVars object
% ---------------------------------
    
   



