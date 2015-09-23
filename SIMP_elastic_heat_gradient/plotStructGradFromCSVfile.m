
clear
clc
close all

% material properties Object
matProp = MaterialProperties;

% settings object
settings = Configuration;


% plotting tool
plotter = plotResults;
settings.plotFinal = 1; % set to final plotting mode. 
settings.plotToCSVFile = 0; % do not replot to the file. 

% set the design var object. 
designVars = DesignVars(settings);


recvid = 1;




lsitOfFolders = ls( 'out*');
count = 0;
for folderS = lsitOfFolders'
    folderS = folderS';
    folder = folderS(1,:)
    
    
    if recvid==1
        videoOut = strcat(folder,'/resultsOuts.avi');
        vidObj = VideoWriter(videoOut);    %Prepare the new file for video
        vidObj.FrameRate = 5;
        vidObj.Quality = 100;
        open(vidObj);
        vid=1;
    end
    
    settings.w1 = count/10;
    count = count+1;

    for i = 1:256
        name = sprintf('./%s/gradAndStuct%i.csv',folder, i);

        % if the file does not exist, then save the final graph, and break. 
        if exist(name, 'file') == 0
            
           %  nameGraph = sprintf('./%s/gradTopOptimization%i',folder, i);
            nameGraph = sprintf('./gradTopOptimization%f.png', settings.w1);
            print(nameGraph,'-dpng')
            break;
        end
        name

        structGradArray = csvread(name);
        figure(1)

        
        plotter. ActualPlotStructGradArray(structGradArray, settings,matProp,i)
        %designVars.


             if recvid==1
                 drawnow
                F(vid) = getframe(figure(1)); %#ok<AGROW> %Get frame of the topology in each iteration
                writeVideo(vidObj,F(vid)); %Save the topology in the video
                vid=vid+1;
            end

    end


    if recvid==1
        close(vidObj);  %close video
    end
end