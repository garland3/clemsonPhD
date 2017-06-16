classdef VideoManager
    %VideoManager Manages Video Functions
    
    properties
    end
    
    methods
        % --------------------
        % Call 1 time at beginning. 
        % --------------------
        function [vidObj, framedNumber] = InitializeVideo(obj, config,name)
            vidObj = 0;
            framedNumber= 1;
            if config.recvid==1
                videoOut =name;
                vidObj = VideoWriter(videoOut);    %Prepare the new file for video
                vidObj.FrameRate = 5;
                vidObj.Quality = 100;
                open(vidObj);            
            end
        end
        
        % --------------------
        % Call each time you want to add a frame
        % --------------------
        function [framedNumber, F]  = RecordFrame(obj, config,framedNumber, F,vidObj)
            if config.recvid==1
                drawnow
                F(framedNumber) = getframe(figure(1)); % %Get frame of the topology in each iteration
                 writeVideo(vidObj,F(framedNumber)); %Save the topology in the video
                framedNumber=framedNumber+1;
            end
        end
        
         % --------------------
        % Call each time you want to add a frame
        % --------------------
        function  CloseVideo(obj, config, F,vidObj)
            
            if config.recvid==1                
                close(vidObj);  %close video
            end
            
        end
    end
    
end

