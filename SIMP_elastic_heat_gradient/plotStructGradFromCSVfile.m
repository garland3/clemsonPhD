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

plotJustFinalResults = 1; % 1 for yes, 0 for no. If Yes, the dont' make a video, and only plot the final design and design metrics

recvid = 1;

if plotJustFinalResults ==1
    settings.doPlotVolFractionDesignVar = 0;
    settings.doPlotTopologyDesignVar = 0;
    settings.doPlotHeat = 0;
    settings.doPlotHeatSensitivityTopology = 0;
    settings.doPlotStress = 0;
    settings.doPlotFinal = 1;
    settings.doSaveDesignVarsToCSVFile = 0; % set to 1 to write plotFinal csv file instead
    recvid = 0;
end

% set the design var object.
designVars = DesignVars(settings);
plotter.CountPlots(settings);

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
    
    if plotJustFinalResults ==1
        % find the last iteration in this folder. Then plot it.
        finaliterationNumber = 1;
        for i = 1:6000
             nameTopology = sprintf('./%s/topDensity%i.csv',folder, i);
              if exist(nameTopology, 'file') == 0
                  finaliterationNumber = i-1;
                  break;
              end
        end
        status = plotter.plotParticularIterationNumInFolder( folder, finaliterationNumber,designVars, settings,matProp);
        nameGraph = sprintf('./gradTopOptimization%f.png', settings.w1);
        print(nameGraph,'-dpng')
        response = fig2plotly()
        plotly_url = response.url;
        
    else
  
         % Loop over the iteration data saved as .csv files and plot it. 
        for i = 1:6000


            %-----------------------
            % Read the actual files and plot
            % ---------------------
            status = plotter.plotParticularIterationNumInFolder( folder, i,designVars, settings,matProp);
            if status ==-1
                break;
            end

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
end


% designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];