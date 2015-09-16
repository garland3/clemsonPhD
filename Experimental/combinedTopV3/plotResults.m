classdef plotResults
    methods
        function obj = plotResults
        end
        
        function plotTopAndFraction(obj,designVars)
         
           %  ----------------------------
           % Plot Topology Optimization DENSITIES  
           %  ----------------------------
           subplot(2,2,1);
           colormap(winter); 
           imagesc(designVars.x); 
           set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
           % image(x); 
           axis equal; axis tight; axis off;
            
           %  ----------------------------
           % Plot Volume Fraction Optimization
           %  ----------------------------
            
           subplot(2,2,2);
           colormap(winter); 
           imagesc(designVars.w); 
           set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
           % image(x); 
           axis equal; axis tight; axis off;pause(1e-6);
           %    F(vid) = getframe; %Get frame of the topology in each iteration
           %      writeVideo(vidObj,F(vid)); %Save the topology in the video
           %      vid=vid+1;
            
        end
    end
end