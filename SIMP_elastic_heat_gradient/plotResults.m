classdef plotResults
    methods
        function obj = plotResults
        end
        
        function plotTopAndFraction(obj,designVars, settings, matProp, loopNumb)
         
           %  ----------------------------
           % Plot Topology Optimization DENSITIES  
           %  ----------------------------
           if(settings.plotFinal ~=1)
                 figure(1)
               subplot(2,2,1);

               %colormap(winter); 
               imagesc(designVars.x); 
                 title('Topology Opt');
               set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
               % image(x); 
               axis equal; axis tight; axis off;


               %  ----------------------------
               % Plot Volume Fraction Optimization
               %  ----------------------------
    %              subplot(2,2,1);
    %            imagesc(designVars.w);
    %             title('w');
    %             colorbar;
    %            set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
    %            % image(x); 
    %            axis equal; axis tight; axis off;pause(1e-6);

                  figure(1)
               subplot(2,2,2);

             %  colormap(winter); 
               imagesc(designVars.w*matProp.E_material1+(1-designVars.w)*matProp.E_material2);
                title('Effective Elastic Fraction Opt');
                colorbar;
               set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
               % image(x); 
               axis equal; axis tight; axis off;pause(1e-6);


               %    F(vid) = getframe; %Get frame of the topology in each iteration
               %      writeVideo(vidObj,F(vid)); %Save the topology in the video
               %      vid=vid+1;
           end
           
%            if settings.mode ==1 || 2||3
              obj.PlotStrucAndGrad(designVars, settings, matProp,loopNumb)
%            end
            
        end
        
        % Plot a graph showing the topology optimization results and the
        % gradient optimization results. 
        function PlotStrucAndGrad(obj, designVars, settings, matProp,loopNumb)
            
            % ----------------------------------
             % Generate the matrix used for ploitting. 
             % ----------------------------------
            structGradArrayElastic(1:settings.nely,1:settings.nelx)  = 0; % Initialize
           %  structGradArrayHeat(1:settings.nely,1:settings.nelx)  = 0; % Initialize
            
            
            
         
             for i = 1:settings.nelx
                for j = 1:settings.nely
                    
                    x_local = designVars.x(j,i);
                    if(x_local <= settings.voidMaterialDensityCutOff) % if void region
                       % E_atElement(j,i) = E_empty;
                       % K_atElement(i,j) = K_empty;
                       structGradArrayElastic(j,i) = matProp.E_material2-0.01*matProp.E_material2; % make the void region 25 less than the least strong material for plotting purposes
                       structGradArrayHeat(j,i) = matProp.E_material1-1;
                    else % if a filled region
                       volFraclocal = designVars.w(j,i);
%                        volume1 = volume1 +volFraclocal; % sum up the total use of material 1 
%                        volume2 = volume2 + (1- volFraclocal); % sum up the total use of material 2 

                      % K_atElement(i,j) = KheatPLA*volFraclocal+(1-volFraclocal)*KheatNylon;
                      
                      
                       E_heat_atElement = matProp.effectiveHeatProperties(volFraclocal);
                       E_atElement=  matProp.effectiveElasticProperties(volFraclocal);  % simple mixture ratio
                       structGradArrayElastic(j,i) = E_atElement;
                       structGradArrayHeat(j,i) = E_heat_atElement;
                    end
                end
             end
             
             
             TcontourMatrix = 1;

 
  
            if settings.doPlotHeat ==1
                numNodesInRow = settings.nelx +1;
                numNodesInColumn = settings.nely+1;
                
                TcontourMatrix = zeros(numNodesInRow,numNodesInColumn);
                for j = 1:numNodesInColumn % y
                      rowMultiplier = j-1;
                     for i = 1:numNodesInRow % x
                         nodeNumber = i+numNodesInRow*rowMultiplier;
                         TcontourMatrix(i,j) = designVars.U_heatColumn(nodeNumber);

                     end
                end
            end
            
             
             ActualPlotStructGradArray(obj,structGradArrayElastic,TcontourMatrix, settings,matProp,designVars, loopNumb)
             
            
            
%              figure(1)
%             subplot(2,2,4);
%             imagesc(structGradArrayHeat); axis equal; axis tight; axis off;
%             set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
%             % caxis([-1 1 ]);
%             title('Structure and heat condunction Gradient')
%             %colormap winter
%             %  cmap = colormap;
%             rgbSteps = matProp.K_material2- matProp.K_material1;  % plus 1 for 1 below the Enylon for void
%             rgbSteps = rgbSteps*50;
%             % [cmin,cmax] = caxis;
%             caxis([matProp.E_material2+1,matProp.E_material1])
%             map = colormap; % current colormap
%             %map = [colormap(1,1):1/rgbSteps:colormap(1:-1)
%             map(rgbSteps,:) = [1,1,1];
%             map(rgbSteps-1,:) = [1,1,1];
%             for zz =    1:rgbSteps-1
%                 
%                 map(zz,:) = [0,       zz*7/(8*rgbSteps)+1/8,          0.5];
%             end    
%             colormap(map)   
%             colorbar
            
        end
        
        
        function ActualPlotStructGradArray(obj,structGradArrayElastic, temperatureField, settings,matProp,designVars, loopNum)
            
            if(settings.plotToCSVFile ==1)                
                 
                 % -------------------------------
                 % Plot to CSV file
                 % ------------------------------
                % loopNumb
                folderNum = settings.iterationNum;
                  name = sprintf('./out%i/gradAndStuct%i.csv',folderNum, loopNum);
                    csvwrite(name,structGradArrayElastic);
                    
            else
                
                % Plot normally. 
              figure(1)
             if(settings.doPlotHeat ==1)
                    subplot(1,2,1);
                    maxH = max(max(temperatureField));
                    minH = min(min(temperatureField));
                    contour(designVars.XLocations,designVars.YLocations,temperatureField,'ShowText','on');
                    % imagesc(temperatureField); axis equal; axis tight; axis off;
                    set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/an
                    
                     titleText = sprintf('Heat,\n max=%f, min = %f',maxH,minH);
                  
                    title(titleText);
                    
                     subplot(1,2,2);

             else
                   subplot(1,1,1);
                 %  subplot(1,1,1);
             end
            imagesc(structGradArrayElastic); axis equal; axis tight; axis off;
            set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
            % caxis([-1 1 ]);
            titleText = sprintf('Structure and Elastic Mod Gradient,\n w1=%f, iter = %i',settings.w1,loopNum);
            title(titleText)
            %colormap winter
            %  cmap = colormap;
            rgbSteps = (matProp.E_material1- matProp.E_material2);  % plus 1 for 1 below the Enylon for void
            rgbSteps = 100;
            %rgbSteps = rgbSteps*50;
            % [cmin,cmax] = caxis;
            caxis([matProp.E_material2-0.01*matProp.E_material2,matProp.E_material1])
            map = colormap; % current colormap
            %map = [colormap(1,1):1/rgbSteps:colormap(1:-1)
            map(1,:) = [1,1,1];
            for zz =    2:rgbSteps
                
                map(zz,:) = [0,       zz*7/(8*rgbSteps)+1/8,          0.5];
            end    
            colormap(map)   
            colorbar
            
            
            drawnow
            
                 
             end
        end
    end
end