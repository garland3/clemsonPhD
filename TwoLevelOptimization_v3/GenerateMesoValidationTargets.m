function [DV] = GenerateMesoValidationTargets(DV,config,matProp,step)

temp=config.validationGridSizeNelx-1;
scale = 1;
ExxVector =0:scale*matProp.E_material1/temp:matProp.E_material1*scale;
EyyVector =0:scale*matProp.E_material1/temp:matProp.E_material1*scale;
thetaVector = 0:(pi/2)/temp:pi/2;
[ExxValues, EyyValues, ThetaValues] = meshgrid(ExxVector,EyyVector,thetaVector);

if(step==1)
    fprintf('Generate target values for meso validation mode\n');
    Usize = size(config.loadingCase,2);
    nnodes =(config.nelx+1)*(config.nely+1)*2;
    DV.U=ones(Usize,nnodes);% 1; % so it is not NULL. Just assign a placeholder value
    DV.sensitivityElastic=1;
    DV.sensitivityElasticPart2=1;
    DV.w=1;
    DV.lambda1=1;
    
    % DV.x = ones(config.nely, config.nelx);
    
    % totalValidationProblems=config.nely*config.nelx;
    % numSegments = floor(totalValidationProblems^(1/3));
    % numSegmentsExx = numSegments;
    % numSegmentsTheta = floor(totalValidationProblems/(numSegmentsExx^2));
    %
    % numSegmentsExx=numSegmentsExx-1;
    % % numSegmentsTheta=numSegmentsTheta-1;
    
    
    
    
    [ t1 t2 t3]=size(ExxValues);
    ne = t1*t2*t3;
    simpDensity = ones(ne,1);
    
    ExxValues = reshape(ExxValues,[],1);
    EyyValues = reshape(EyyValues,[],1);
    ThetaValues = reshape(ThetaValues,[],1);
    
    % xValues
    % % save the density field
    % outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
    % csvwrite(outname,DV.x);
    folderNum=0;
    % save the Exx field
    mm_iteration=1;
    outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
    
    csvwrite(outname,ExxValues);
    
    % save the Eyy field
    outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
    csvwrite(outname,EyyValues);
    
    % save the Theta field
    outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
    csvwrite(outname,ThetaValues);
    
    % Save a placeholder file.
    outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,mm_iteration);
    csvwrite(outname,[1 2]);
    
    outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
    csvwrite(outname,simpDensity);
elseif(step==2)
    % -------------------------
    % Show all the different designs in a complete structure
    % -------------------------
    
    
    %     postProcess = 1;
    %
    %     close all
    %     p = plotResults;
    %
    %
    %     temp = config.mesoAddAdjcentCellBoundaries;
    %     config.mesoAddAdjcentCellBoundaries=0;
    %     macroElementProps = macroElementProp;
    %     macroElementProps.targetDensity=0.5; % make up some value for now.
    %     [DVMeso, mesoconfig] = GenerateDesignVarsForMesoProblem(config,1,macroElementProps);
    %     config.mesoAddAdjcentCellBoundaries=temp;
    %
    %     mesoconfig.doUseMultiElePerDV =config.doUseMultiElePerDV;
    %     numTilesX=config.numTilesX;
    %     numTilesY = config.numTilesY;
    
    % Generate huge area
    %      temp=config.validationGridSizeNelx-1;
    %     totalX=config.nelx*mesoconfig.nelx*numTilesX
    %     totalY=config.nely*mesoconfig.nely*numTilesY
    
    
    %     ne = config.nelx*config.nely; % number of elements
    %    folderNum=0;
    %    mm_iteration=1;
    %   outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
    %   ExxValues2= csvread(outname);
    %
    % t1= size( ExxVector,2) ;% =0:scale*matProp.E_material1/temp:matProp.E_material1*scale;
    % t2= size( EyyVector,2) ;
    % t3= size( thetaVector,2);
    %
    %
    %
    %   ExxVector
    %   ExxValues3=reshape(ExxValues2,[t1 t2 t3]);
    
    
    %--------------------------------------------
    % Get the density field
    %     %--------------------------------------------
    %    macro_meso_iteration = config.macro_meso_iteration;
    %     %macroElementProps = macroElementProp;
    %     % macroElementProps.elementNumber = e;
    %     folderNum = config.iterationNum;
    % GET the saved element to XY position map (needed for x and w vars retrival)
    %     outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
    %     elementXYposition=csvread(outname);
    % Get the density field
    %     outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,macro_meso_iteration);
    %     xxx = csvread(outname);
    %
    
    % -----------------------------
    % -- PLOT the designs
    % --------------------------
    
    makevide=1;
    if(makevide ==1)
        config.recvid=1;
        video = VideoManager;
        [vidObj, framedNumber] = video.InitializeVideo( config,'./mesoValidation.avi');
        F=getframe();
    end
    
    config.   numTilesX = 1;
    config.     numTilesY = 1;
    
    config.nelx=config.validationGridSizeNelx;
    config.nely=config.validationGridSizeNelx;
    
    configMeso = config;
    configMeso.nelx=configMeso.nelxMeso;
    configMeso.nely=configMeso.nelxMeso;
    
    p = plotResults;
    % plot each theta on a different plot, and make a video.
    elementCount=1;
    for i = 1:size(thetaVector,2)
        % make a blank huge array
        totalY = config.validationGridSizeNelx* config.nelyMeso
        completeStruct = zeros(totalY,totalY);
        
        thetaValue = thetaVector(i);
        xxx=ones(config.validationGridSizeNelx,config.validationGridSizeNelx);
        for j = 1:(config.validationGridSizeNelx)
            for k = 1:(config.validationGridSizeNelx)
                fprintf('element %i  with theta %f\n',elementCount,thetaValue);
                macroElementProps.elementNumber=elementCount;
                %                 results = elementXYposition(macroElementProps.elementNumber,:);
                macroElementProps.yPos = j;
                macroElementProps.xPos = k;
                
                x=GetMesoUnitCellDesignFromCSV(config,elementCount);
                DV.x = x;
                
                step = 1;
                completeStruct= TileMesoStructureV2(configMeso,config, DV,macroElementProps,xxx,completeStruct,step);
                elementCount=elementCount+1;
            end
        end
        completeStruct( completeStruct>1)=1;
        
        completeStruct(completeStruct>config.voidMaterialDensityCutOff)=1;
        completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;
        
        size(completeStruct)
        
        plotname = sprintf('meso designs with theta %f',thetaValue);
        p.PlotArrayGeneric( completeStruct, plotname)
          xlabel('Eyy');
                ylabel('Exx');
                draw now
        
        %     p.PlotArrayGeneric( completeStruct, plotname)
        %     rgbSteps = 100;  caxis([0,1]);
        %     map = colormap; % current colormap
        %     middPoint = floor(rgbSteps/4);
        %     map(1:middPoint,:) = [ones(middPoint,1),ones(middPoint,1),ones(middPoint,1)];
        %     for zz =    middPoint:rgbSteps
        %         map(zz,:) = [0,               1- zz/rgbSteps, 0.5];
        %     end
        %     colormap(map)
        %     %     colorbar
        %     freezeColors
        nameGraph = sprintf('./mesoDesignsWithThetaIndex%i', i);
        %         print(nameGraph,'-dpng', '-r1200')
        print(nameGraph,'-dpng')
        
        
        if(makevide ==1)
            [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
        end
    end
    
    if(makevide ==1)
        video.CloseVideo( config, F,vidObj)
    end
    
    %     if(postProcess==1)
    %         for e = 1:ne
    %             fprintf('step 2element %i of %i\n',e,ne);
    %             macroElementProps.elementNumber=e;
    %             results = elementXYposition(macroElementProps.elementNumber,:);
    %             macroElementProps.yPos = results(1);
    %             macroElementProps.xPos = results(2);
    %             macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
    %
    %             % Check if void
    %             %if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff)
    %             step = 2;
    %             completeStruct= TileMesoStructureV2(mesoconfig,config, DVMeso,macroElementProps,xxx,completeStruct,step);
    %             %end
    %         end
    %     end
    
    
    
    %if(config.multiscaleMethodCompare~=1)
    % set the max value to be 1
    
    %     if (config.generateCompleteStructureCSV==1)
    %         outname = sprintf('./completeStucture%f_macroIteration_%i.csv', config.w1,config.macro_meso_iteration);
    %         csvwrite(outname,completeStruct);
    %     end
    
end

%ExxValues=padarray(ExxValues,
% ExxValues=reshape(ExxValues,1,[]);
% [t1 t2]=size(ExxValues);
% ExxValues=padarray(ExxValues,[0 totalValidationProblems-t2],'post');
% ExxValues=reshape(ExxValues,config.nely,config.nelx);
%
% EyyValues=reshape(EyyValues,1,[]);
% EyyValues=padarray(EyyValues,[0 totalValidationProblems-t2],'post');
% EyyValues=reshape(EyyValues,config.nely,config.nelx);
%
% ThetaValues=reshape(ThetaValues,1,[]);
% ThetaValues=padarray(ThetaValues,[0 totalValidationProblems-t2],'post');
% ThetaValues=reshape(ThetaValues,config.nely,config.nelx);
%
% xSimp = ones(t1, t2);
% xSimp=reshape(xSimp,1,[]);
% xSimp=padarray(xSimp,[0 totalValidationProblems-t2],'post');
% xSimp=reshape(xSimp,config.nely,config.nelx);
%
% DV.Exx = ExxValues;
% DV.Eyy = EyyValues;
% DV.t = ThetaValues;
% DV.x = xSimp;

% SaveMacroProblemStateToCSV(config,DV,matProp)

end
