clear
clc
close all
% material properties Object
matProp = MaterialProperties;
% settings object
settings = Configuration;
% plotting tool
plotter = plotResults;
settings.doPlotVolFractionDesignVar = 1;
settings.doPlotTopologyDesignVar = 1;
settings.doPlotHeat = 0;
settings.doPlotHeatSensitivityTopology = 0;
settings.doPlotStress = 0;
settings.doPlotFinal = 1;
settings.doSaveDesignVarsToCSVFile = 0; % set to 1 to write plotFinal csv file instead

% set the design var object.
designVars = DesignVars(settings);
recvid = 1;

lsitOfFolders = ls( 'out*');
count = 0;
for folderS = lsitOfFolders'
    folderS = folderS';
    folder = folderS(1,:);
    folder = strtrim(folder)    
    %  totalStringLength = numel(folder);
    iterationNum = str2num(folder(4:end));    
    
    if recvid==1
        videoOut = strcat(folder,'/resultsOuts.avi');
        vidObj = VideoWriter(videoOut);    %Prepare the new file for video
        vidObj.FrameRate = 5;
        vidObj.Quality = 100;
        open(vidObj);
        vid=1;
    end
    
    settings.w1 = iterationNum/10;
    %count = count+1;    
    folderNum = iterationNum;
    outname = sprintf('./out%i/storeOptimizationVar.csv',folderNum);
    designVars.storeOptimizationVar = csvread(outname);
  
    
    for i = 1:256
        nameTopology = sprintf('./%s/topDensity%i.csv',folder, i);
        nameVolFractionGrad = sprintf('./%s/volFractionVar%i.csv',folder, i);
        
        % if the file does not exist, then save the final graph, and break.
        if exist(nameTopology, 'file') == 0
            %  nameGraph = sprintf('./%s/gradTopOptimization%i',folder, i);
            nameGraph = sprintf('./gradTopOptimization%f.png', settings.w1);
            print(nameGraph,'-dpng')
            break;
        end
        [nameTopology nameVolFractionGrad]
        
        %-----------------------
        % Read the actual files
        % ---------------------
        designVars.x = csvread(nameTopology);
        designVars.w = csvread(nameVolFractionGrad);        
        FEACalls = i;
        plotter.plotTopAndFraction(designVars,  settings, matProp, FEACalls); % plot the results.
        
        if recvid==1
            drawnow
            F(vid) = getframe(figure(1)); % %Get frame of the topology in each iteration
            writeVideo(vidObj,F(vid)); %Save the topology in the video
            vid=vid+1;
        end        
    end
    
    if recvid==1
        close(vidObj);  %close video
    end
end


% designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];