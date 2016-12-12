classdef plotResults
    methods
        function obj = plotResults
        end
        
        function plotTopAndFraction(obj,designVars, settings, matProp, loopNumb)
            
            % ------------------------------------------------------
            % Count how many things to plot.
            % Geneate the subplot number
            % Use if statements for each potential plot
            % ------------------------------------------------------
            
            numberOfPlots=obj.CountPlots(settings);
            bestSquareSize = ceil(sqrt(numberOfPlots));
            plotDim2 = ceil(numberOfPlots/sqrt(numberOfPlots));
            plotDim1=bestSquareSize;
            clf
            plotcount = 1;
            figure(1)
            if numberOfPlots == 2
                plotDim2 =1;  plotDim1=2;
            end
            
            %  ----------------------------
            % Plot Topology Optimization DENSITIES design vars
            %  ----------------------------
            if(settings.doPlotTopologyDesignVar ==1)
                titleText = 'Topology Opt density';
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotArrayGeneric(designVars.x,titleText)
            end
            
            %  ----------------------------
            % Plot Topology Optimization volume fraction design vars
            %  ----------------------------
            if(settings.doPlotVolFractionDesignVar ==1)
                titleText = 'Vol Fraction Design Var';
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotArrayGeneric(designVars.w,titleText)
            end
            
            %  ----------------------------
            % Plot Heat
            %  ----------------------------
            if(settings.doPlotHeat ==1)
                figure(1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotTemperatures(settings, designVars)
                colorbar
            end
            
            %  ----------------------------
            % Plot Stress
            %  ----------------------------
            if(settings.doPlotStress ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = sprintf('Von Mises Stress, maxF = %f, maxU=%f', designVars. maxF, designVars.maxU );
                obj.PlotArrayGeneric(designVars.totalStress,titleText)
            end
            
            %  ----------------------------
            % Plot Final
            %  ----------------------------
            if(settings.doPlotFinal == 1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotStrucAndGrad(designVars, settings, matProp,loopNumb)
            end
            
            
            %  ----------------------------
            % Plot Heat topology sensitivity, temp2
            %  ----------------------------
            if(settings.doPlotHeatSensitivityTopology == 1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = sprintf('Heat Top Sensit');
                obj.PlotArrayGeneric(designVars.temp2,titleText)
            end
                   
            if(settings.doSaveDesignVarsToCSVFile ==1)                
                % -------------------------------
                % Plot to CSV file
                % ------------------------------
                % loopNumb
                folderNum = settings.iterationNum;
                name = sprintf('./out%i/topDensity%i.csv',folderNum, loopNumb);
                csvwrite(name,designVars.x);
                
                name = sprintf('./out%i/volFractionVar%i.csv',folderNum, loopNumb);
                csvwrite(name,designVars.w);                
            end
            
            %  ----------------------------
            % Plot Design Update Metrics
            %  ----------------------------
            if(settings.doPlotMetrics == 1)
                 subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotDesignMetrics(designVars, settings, matProp, loopNumb);
            end
            
              drawnow
        end
        
        function PlotDesignMetrics(obj, designVars, settings, matProp, loopNumb)
            %  designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            x = 1:loopNumb;
           % y1 = designVars.storeOptimizationVar(1:loopNumb,1)';
            y2 = designVars.storeOptimizationVar(1:loopNumb,2)'; % Elastic Compliance
            y3 = designVars.storeOptimizationVar(1:loopNumb,3)'; % Heat Compliance         
            y4 = designVars.storeOptimizationVar(1:loopNumb,4)'; % volume fraction material 1
            y5 = designVars.storeOptimizationVar(1:loopNumb,5)'; % volume fraction material 2
            
            y2 = y2/max(y2); % normalize to make plotting nice
            y3 = y3/max(y3); % normalize to make plotting nice
            plot(x, y4,'y', x, y5, 'm', x, y2, 'c', x, y3, 'r')            
            legend('vol1', 'vol2','Elast Obj','Heat Obj')
        end
        
        
        % Plots an array
        function PlotArrayGeneric(obj, array, titleText)
                imagesc(array); axis equal; axis tight; axis off;
            % colormap winter
            set(gca,'YDir','normal');
            title(titleText);
            colorbar
        end        
        
        % --------------------------------------------
        % Plots temperature countours
        % --------------------------------------------
        function PlotTemperatures(obj,settings, designVars)
            % Set up temperature array
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
            
            % Plot Contours
            averageTemp = mean2(TcontourMatrix);
            contour(designVars.XLocations,designVars.YLocations,TcontourMatrix,'ShowText','on');
            set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/an
            titleText = sprintf('Heat,\n average T=%f',averageTemp);
            title(titleText);
            freezeColors
            
        end
        
        
        
        % ==============================================
        % -----------------------------------------------
        % Counts the number of Plots needed based on the configuration
        % settings
        % -----------------------------------------------
        function [numberOfPlots]=CountPlots(obj,settings)
            
            numberOfPlots = 0;
            
            % Plotting information
            %         doPlotVolFractionDesignVar = 0;
            %         doPlotTopologyDesignVar = 1;
            %         doPlotHeat = 1;
            %         doPlotHeatSensitivityTopology = 1;
            %         doPlotStress = 0;
            %         doPlotFinal = 1;
            
            if(settings.doPlotHeat ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
            if(settings.doPlotStress ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
            if(settings.doPlotFinal ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
            if(settings.doPlotTopologyDesignVar ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
            if(settings.doPlotVolFractionDesignVar ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
            if(settings.doPlotHeatSensitivityTopology ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
            if(settings.doPlotMetrics ==1)
                numberOfPlots = numberOfPlots+1;
            end
            
        end
        
        
        
        
        % -------------------------------------
        % Plot a graph showing the topology optimization results and the
        % gradient optimization results.
        % -------------------------------------
        function PlotStrucAndGrad(obj, designVars, settings, matProp,loopNumb)
            
            % ----------------------------------
            % Generate the matrix used for ploitting.
            % ----------------------------------
            % structGradArrayElastic(1:settings.nely,1:settings.nelx)  = 0; % Initialize
            structGradDesignVarArray(1:settings.nely,1:settings.nelx)  = 0; % Initialize
            %  structGradArrayHeat(1:settings.nely,1:settings.nelx)  = 0; % Initialize
            
            %smallerE = min( matProp.E_material2, matProp.E_material1);
            %largerE = max( matProp.E_material2, matProp.E_material1);
            %diffE = largerE-smallerE;
            
            for i = 1:settings.nelx
                for j = 1:settings.nely
                    x_local = designVars.x(j,i);
                    if(x_local <= settings.voidMaterialDensityCutOff) % if void region
                        % structGradArrayElastic(j,i) =0; % make the void region 10 less than the least strong material for plotting purposes
                        structGradDesignVarArray(j,i) =-0.5;
                    else % if a filled region
                        volFraclocal = designVars.w(j,i);
                        % E_atElement=  matProp.effectiveElasticProperties(volFraclocal,settings);  % simple mixture ratio
                        %structGradArrayElastic(j,i) = E_atElement;
                        structGradDesignVarArray(j,i) = volFraclocal;
                    end
                end
            end
            ActualPlotStructGradArray(obj,structGradDesignVarArray, settings,matProp,designVars, loopNumb)
        end
        
        % -------------------------------------
        % ActualPlotStructGradArray
        % -------------------------------------
        function ActualPlotStructGradArray(obj,structGradDesignVarArray,  settings,matProp,designVars, loopNum)
            
            %             if(settings.doPlotFinalToCSVFile == 1)
            %                 % -------------------------------
            %                 % Plot to CSV file
            %                 % ------------------------------
            %                 % loopNumb
            %                 folderNum = settings.iterationNum;
            %                 name = sprintf('./out%i/gradAndStuct%i.csv',folderNum, loopNum);
            %                 csvwrite(name,structGradDesignVarArray);
            %             else
            % Plot normally.
            %  figure(1)
            %                 if(settings.doPlotHeat ==1 && settings.doPlotSensitivityComparison~=1)
            %
            %                 elseif(settings.doPlotSensitivityComparison ==1 && settings.doPlotHeat ==1)
            %                     % PLot heat
            %                     subplot(2,2,1);
            %                     averageTemp = mean2(temperatureField);
            %                     contour(designVars.XLocations,designVars.YLocations,temperatureField,'ShowText','on');
            %                     set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/an
            %                     titleText = sprintf('Heat,\n average T=%f',averageTemp);
            %                     title(titleText);
            %                     freezeColors
            %                     % PLot comparison
            %                     subplot(2,2,3);
            %                     tStress= max(max(abs(designVars.totalStress)));
            %                     fprintf('max Stress=%f\n', tStress);
            %                 else
            %                     subplot(1,1,1);
            %                     %  subplot(1,1,1);
            %                 end
            
            imagesc(structGradDesignVarArray); axis equal; axis tight; axis off;
            set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
            % caxis([-1 1 ]);
            titleText = sprintf('Structure and vol fraction,\n w1=%f, iter = %i',settings.w1,loopNum);
            title(titleText)
            rgbSteps = 100;
            %smallerE = min( matProp.E_material2, matProp.E_material1);
            %largerE = max( matProp.E_material2, matProp.E_material1);
            
            caxis([-0.02,1])
            
            map = colormap; % current colormap
            %map = [colormap(1,1):1/rgbSteps:colormap(1:-1)
            map(1,:) = [1,1,1];
            for zz =    2:rgbSteps
                map(zz,:) = [0,       zz*7/(8*rgbSteps)+1/8,          0.5];
            end
            colormap(map)
            colorbar
            freezeColors
            %             end
        end
        
        function status = plotParticularIterationNumInFolder(obj, folder, i, designVars, settings,matProp)
            status = 1;
            nameTopology = sprintf('./%s/topDensity%i.csv',folder, i);
            nameVolFractionGrad = sprintf('./%s/volFractionVar%i.csv',folder, i);

            % if the file does not exist, then save the final graph, and break.
            if exist(nameTopology, 'file') == 0
                %  nameGraph = sprintf('./%s/gradTopOptimization%i',folder, i);
                nameGraph = sprintf('./gradTopOptimization%f.png', settings.w1);
                print(nameGraph,'-dpng')
               %  p = fig2plotly; needs more work, this looks horrible. 
               status = -1;
                return;
            end
            [nameTopology nameVolFractionGrad]
            %-----------------------
            % Read the actual files
            % ---------------------
            designVars.x = csvread(nameTopology);
            designVars.w = csvread(nameVolFractionGrad);      
            [settings.nelx]   = size(  designVars.x,2);
              [   settings.nely] = size(  designVars.x,1);

            FEACalls = i;
            obj.plotTopAndFraction(designVars,  settings, matProp, FEACalls); % plot the results.

            
            
        end
        
        
       
    end
end