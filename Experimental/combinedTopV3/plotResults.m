classdef plotResults
    methods
        function obj = plotResults
        end
        
        function plotTopAndFraction(obj,designVars, settings, matProp)
         
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
            colorbar;
           set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
           % image(x); 
           axis equal; axis tight; axis off;pause(1e-6);
          
           %    F(vid) = getframe; %Get frame of the topology in each iteration
           %      writeVideo(vidObj,F(vid)); %Save the topology in the video
           %      vid=vid+1;
           
           if settings.mode ==3
               obj.PlotStrucAndGrad(designVars, settings, matProp)
           end
            
        end
        
        % Plot a graph showing the topology optimization results and the
        % gradient optimization results. 
        function PlotStrucAndGrad(obj, designVars, settings, matProp)
            
            % ----------------------------------
             % Generate the matrix used for ploitting. 
             % ----------------------------------
            structGradArray(1:settings.nely,1:settings.nelx)  = 0; % Initialize
            
         
             for i = 1:settings.nelx
                for j = 1:settings.nely
                    
                    x_local = designVars.x(j,i);
                    if(x_local <= settings.voidMaterialDensityCutOff) % if void region
                       % E_atElement(j,i) = E_empty;
                       % K_atElement(i,j) = K_empty;
                       structGradArray(j,i) = matProp.E_material2-matProp.E_material2*0.25; % make the void region 25 less than the least strong material for plotting purposes
                    else % if a filled region
                       volFraclocal = designVars.w(j,i);
%                        volume1 = volume1 +volFraclocal; % sum up the total use of material 1 
%                        volume2 = volume2 + (1- volFraclocal); % sum up the total use of material 2 

                      % K_atElement(i,j) = KheatPLA*volFraclocal+(1-volFraclocal)*KheatNylon;
                       E_atElement= matProp.E_material1*volFraclocal+(1-volFraclocal)*matProp.E_material2;  % simple mixture ratio
                       structGradArray(j,i) = E_atElement;
                    end
                end
            end
            
            subplot(2,2,3);
            imagesc(structGradArray); axis equal; axis tight; axis off;
            set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
            % caxis([-1 1 ]);
            title('Structure and Elastic Mod Gradient')
            colormap winter
            %  cmap = colormap;
            rgbSteps = matProp.E_material1- matProp.E_material2 +0.1 ; % plus 1 for 1 below the Enylon for void
            % [cmin,cmax] = caxis;
            caxis([matProp.E_material2-0.5,matProp.E_material1])
            map = colormap; % current colormap
            %map = [colormap(1,1):1/rgbSteps:colormap(1:-1)
            for zz =    1:rgbSteps
                map(1,:) = [1,1,1];
                map(zz,:) = [0,       zz*7/(8*rgbSteps)+1/8,          0.5];
            end    
            colormap(map)   
            colorbar
            
        end
    end
end