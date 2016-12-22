classdef plotResults
    methods
        function obj = plotResults
        end
        
        function plotTopAndFraction(obj,DV, config, matProp, loopNumb)
            
            % ------------------------------------------------------
            % Count how many things to plot.
            % Geneate the subplot number
            % Use if statements for each potential plot
            % ------------------------------------------------------
            
            numberOfPlots=obj.CountPlots(config);
            bestSquareSize = ceil(sqrt(numberOfPlots));
            plotDim2 = ceil(numberOfPlots/sqrt(numberOfPlots));
            plotDim1=bestSquareSize;
            clf
            plotcount = 1;
            figure(1)
            set(gcf, 'Position', get(0, 'Screensize'));
            %             figure('Name', '1', 'units','normalized','outerposition',[0 0 1 1])
            if numberOfPlots == 2
                plotDim2 =1;  plotDim1=2;
            end
            
            %  ----------------------------
            % Plot Topology Optimization DENSITIES design vars
            %  ----------------------------
            if(config.doPlotTopologyDesignVar ==1)
                titleText = 'Topology Opt density';
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotArrayGeneric(DV.x,titleText)
            end
            
            %  ----------------------------
            % Plot volume fraction design vars
            %  ----------------------------
            if(config.doPlotVolFractionDesignVar ==1)
                titleText = 'Vol Fraction Design Var';
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotArrayGeneric(DV.w,titleText)
            end
            
            %  ----------------------------
            % Plot Heat
            %  ----------------------------
            if(config.doPlotHeat ==1)
                figure(1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotTemperatures(config, DV)
                colorbar
            end
            
            %  ----------------------------
            % Plot Stress
            %  ----------------------------
            if(config.doPlotStress ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = sprintf('Von Mises Stress, maxF = %f, maxU=%f', DV. maxF, DV.maxU );
                obj.PlotArrayGeneric(DV.totalStress,titleText)
            end
            
            %  ----------------------------
            % Plot Orthogonal Distribution Var
            %  ----------------------------
            %             if(config.doPlotOrthDistributionVar ==1)
            %                 subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
            %
            %                 obj.PlotOrthDistributionAsArrows(DV, config, loopNumb);
            %             end
            
            %  ----------------------------
            % Plot Exx and Eyy as Arrows
            %  ----------------------------
            if(config.doPlotEyyExxArrows==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotExxEyyVarAsArrows(DV, config, loopNumb);
            end
            
             
            %  ----------------------------
            % Plot Rotation Var Value
            %  ----------------------------
            if(config.doPlotRotationValue ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotRotationVarAsArrows(DV, config, loopNumb);
            end
            
            %  ----------------------------
            % PLOT COMBINED ROTATION AND ORTHOGONAL DIRECTION GRAPH
            %  ----------------------------
            if(config.doPlotCombinedExxEyyAndRotation ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotOrthAndRotationTogether(DV, config, loopNumb);
            end
            
            %  ----------------------------
            % PLOT Exx
            %  ----------------------------
            if(config.doPlotExx ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = sprintf(' Exx ' );
                obj.PlotArrayGeneric(DV.Exx,titleText)
            end
            
            %  ----------------------------
            % PLOT Eyy
            %  ----------------------------
            if(config.doPlotEyy ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = sprintf(' Eyy ' );
                obj.PlotArrayGeneric(DV.Eyy,titleText)
            end
            
            
            
            
            %  ----------------------------
            % Plot ElasticSensitivity
            %  ----------------------------
            if(config.doPlotElasticSensitivity ==1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = 'doPlotElasticSensitivity';
                obj.PlotArrayGeneric(DV.sensitivityElastic,titleText)
            end
            
            %  ----------------------------
            % Plot Final
            %  ----------------------------
            if(config.doPlotFinal == 1)
                graphichandle = subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotStrucAndGrad(DV, config, matProp,loopNumb,graphichandle)
            end
            
            
            %  ----------------------------
            % Plot Heat topology sensitivity, temp2
            %  ----------------------------
            if(config.doPlotHeatSensitivityTopology == 1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                titleText = sprintf('Heat Top Sensit');
                obj.PlotArrayGeneric(DV.temp2,titleText)
            end
            
            if(config.doSaveDesignVarsToCSVFile ==1)
                % -------------------------------
                % Plot to CSV file
                % ------------------------------
                % loopNumb
                folderNum = config.iterationNum;
                name = sprintf('./out%i/topDensity%i.csv',folderNum, loopNumb);
                csvwrite(name,DV.x);
                
                name = sprintf('./out%i/volFractionVar%i.csv',folderNum, loopNumb);
                csvwrite(name,DV.w);
            end
            
            %  ----------------------------
            % Plot Design Update Metrics
            %  ----------------------------
            if(config.doPlotMetrics == 1)
                subplot(plotDim1,plotDim2,plotcount); plotcount = plotcount + 1;
                obj.PlotDesignMetrics(DV, config, matProp, loopNumb);
            end
            
            drawnow
        end
        
        function PlotDesignMetrics(obj, DV, config, matProp, loopNumb)
            %  DV.c, DV.cCompliance, DV.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            x = 1:loopNumb;
            % y1 = DV.storeOptimizationVar(1:loopNumb,1)';
            y2 = DV.storeOptimizationVar(1:loopNumb,2)'; % Elastic Compliance
            y3 = DV.storeOptimizationVar(1:loopNumb,3)'; % Heat Compliance
            y4 = DV.storeOptimizationVar(1:loopNumb,4)'; % volume fraction material 1
            y5 = DV.storeOptimizationVar(1:loopNumb,5)'; % volume fraction material 2
            
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
        function PlotTemperatures(obj,config, DV)
            % Set up temperature array
            TcontourMatrix = 1;
            if config.doPlotHeat ==1
                numNodesInRow = config.nelx +1;
                numNodesInColumn = config.nely+1;
                TcontourMatrix = zeros(numNodesInRow,numNodesInColumn);
                for j = 1:numNodesInColumn % y
                    rowMultiplier = j-1;
                    for i = 1:numNodesInRow % x
                        nodeNumber = i+numNodesInRow*rowMultiplier;
                        TcontourMatrix(i,j) = DV.U_heatColumn(nodeNumber);
                    end
                end
            end
            
            % Plot Contours
            averageTemp = mean2(TcontourMatrix);
            contour(DV.XLocations,DV.YLocations,TcontourMatrix,'ShowText','on');
            set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/an
            titleText = sprintf('Heat,\n average T=%f',averageTemp);
            title(titleText);
            freezeColors
            
        end
        
        
        
        % ==============================================
        % -----------------------------------------------
        % Counts the number of Plots needed based on the configuration
        % config
        % -----------------------------------------------
        function [numberOfPlots]=CountPlots(obj,config)
            
            %             numberOfPlots = 0;
            
            numberOfPlots=config.doPlotVolFractionDesignVar+ ...
                config.doPlotTopologyDesignVar + ...
                config.doPlotHeat+...
                config.doPlotHeatSensitivityTopology +...
                config.doPlotStress+...
                config.doPlotFinal +...
                config.doPlotMetrics +...
                config.doSaveDesignVarsToCSVFile +...
                config.doPlotAppliedStrain+...
                config.doPlotExx+...
                config.doPlotEyy+...
                config.doPlotElasticSensitivity +...
                config.doPlotRotationValue +...
                config.doPlotCombinedExxEyyAndRotation +...
                config.doPlotEyyExxArrows;
            
            
            
        end
        
        % -------------------------------------
        % Plot the orth distribution Var as arrows
        % -------------------------------------
        function PlotOrthDistributionAsArrows(obj, DV, config, loopNumb)
            titleText = 'Orth Distribution Var Values' ;
            
            [X,Y] = meshgrid(1:config.nelx,1:config.nely);
            dx = DV.d;
            dy =  1-DV.d;
            quiver(X,Y,dx,dy);
            %             hold on
            %             imagesc(DV.d); axis equal; axis tight; axis off;
            title(titleText);
        end
        
        % -------------------------------------
        % Plot the ROTATION VAR AS ARROWS
        % -------------------------------------
        function PlotRotationVarAsArrows(obj, DV, config, loopNumb)
            titleText = 'Rotation Var Values' ;
            [X,Y] = meshgrid(1:config.nelx,1:config.nely);
            dx = cos(DV.t);
            dy = sin(DV.t);
            quiver(X,Y,dx,dy);
            title(titleText);
        end
        
         % -------------------------------------
        % Plot the Exx and Eyy as arrows
        % -------------------------------------
        function   PlotExxEyyVarAsArrows(obj, DV, config, loopNumb)
            titleText = 'PlotExxEyyVarAsArrows' ;
            [X,Y] = meshgrid(1:config.nelx,1:config.nely);
            dx = DV.Exx;
            dy = DV.Eyy;
            quiver(X,Y,dx,dy);
             axis equal
            title(titleText);
        end
      
        
        
        % -------------------------------------
        % Plot the ROTATION VAR and the orthogona vars together AS ARROWS
        % -------------------------------------
        function PlotOrthAndRotationTogether(obj, DV, config, loopNumb)
            
            % Big arrow
            theta1(1:config.nely,1:config.nelx) = zeros(config.nely,config.nelx);
            magnitude1(1:config.nely,1:config.nelx) = zeros(config.nely,config.nelx);
            
            % Little arrow
            theta2(1:config.nely,1:config.nelx) = zeros(config.nely,config.nelx);
            magnitude2(1:config.nely,1:config.nelx) = zeros(config.nely,config.nelx);
            
            for ely = 1:config.nely
                for elx = 1:config.nelx
                    d_local = DV.d(ely,elx);
                    %                          material1Fraction  = DV.w(ely,elx);
                    t_local =  DV.t(ely,elx);
                    x_local = DV.x(ely,elx);
                    
                    if(x_local<config.voidMaterialDensityCutOff)
                        theta1(ely,elx)=0;
                        magnitude1(ely,elx)=0;
                        theta2(ely,elx)=0;
                        magnitude2(ely,elx)=0;
                    else
                        
                        theta1(ely,elx)=t_local;
                        magnitude1(ely,elx)=d_local;%*material1Fraction*x_local;
                        
                        theta2(ely,elx)=t_local+pi/2;
                        magnitude2(ely,elx)=(1- d_local);%;*material1Fraction*x_local;
                        
                    end
                end
            end
            
            titleText = 'Distribution and Rotation of Material, Red = E_xx,Green = E_yy' ;
            [X,Y] = meshgrid(1:config.nelx,1:config.nely);
            
            % RED MAIN DIRECTION ARROW
            % GREEN SECONDARY
            % PLOT green first, so that red will be on top.
            dx2 = cos(theta2).*magnitude2;
            dy2 = sin(theta2).*magnitude2;
            q = quiver(X,Y,dx2,dy2,'Autoscale','off');
            q.Color = 'green';
            
            
            %                 set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/an
            hold on
            
            dx = cos(theta1).*magnitude1;
            dy = sin(theta1).*magnitude1;
            q = quiver(X,Y,dx,dy,'Autoscale','off');
            q.Color = 'red';
            
            axis equal
            %                set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/an
            
            
            
            title(titleText);
            
            
            
            
            
        end
        
        
        
        
        
        
        
        
        % -------------------------------------
        % Plot a graph showing the topology optimization results and the
        % gradient optimization results.
        % -------------------------------------
        function PlotStrucAndGrad(obj, DV, config, matProp,loopNumb, graphicHandle)
            
                   
            structGradDesignVarArray(1:config.nely,1:config.nelx)  = 0; % Initialize
            
             
              structGradDesignVarArray(DV.x<config.voidMaterialDensityCutOff)=-0.5;
                structGradDesignVarArray(DV.x>config.voidMaterialDensityCutOff)=DV.w(DV.x>config.voidMaterialDensityCutOff);
             
            
%             for i = 1:config.nelx
%                 for j = 1:config.nely
%                     x_local = DV.x(j,i);
%                     if(x_local <= config.voidMaterialDensityCutOff) % if void region
%                         % structGradArrayElastic(j,i) =0; % make the void region 10 less than the least strong material for plotting purposes
%                         structGradDesignVarArray(j,i) =-0.5;
%                     else % if a filled region
%                         volFraclocal = DV.w(j,i);
%                         % E_atElement=  matProp.effectiveElasticProperties(volFraclocal,config);  % simple mixture ratio
%                         %structGradArrayElastic(j,i) = E_atElement;
%                         structGradDesignVarArray(j,i) = volFraclocal;
%                     end
%                 end
%             end
            ActualPlotStructGradArray(obj,structGradDesignVarArray, config,matProp,DV, loopNumb, graphicHandle)
        end
        
        % -------------------------------------
        % ActualPlotStructGradArray
        % -------------------------------------
        function ActualPlotStructGradArray(obj,structGradDesignVarArray,  config,matProp,DV, loopNum,graphicHandle)
            
           
            
            imagesc(structGradDesignVarArray); axis equal; axis tight; axis off;
            set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
            % caxis([-1 1 ]);
            titleText = sprintf('Structure and vol fraction,\n w1=%f, iter = %i',config.w1,loopNum);
            title(titleText)
            rgbSteps = 100;
            %smallerE = min( matProp.E_material2, matProp.E_material1);
            %largerE = max( matProp.E_material2, matProp.E_material1);
            
            caxis([-0.02,1])
            
            map = colormap(); % current colormap
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
        
        function status = plotParticularIterationNumInFolder(obj, folder, i, DV, config,matProp)
            status = 1;
            nameTopology = sprintf('./%s/topDensity%i.csv',folder, i);
            nameVolFractionGrad = sprintf('./%s/volFractionVar%i.csv',folder, i);
            
            % if the file does not exist, then save the final graph, and break.
            if exist(nameTopology, 'file') == 0
                %  nameGraph = sprintf('./%s/gradTopOptimization%i',folder, i);
                nameGraph = sprintf('./gradTopOptimization%f.png', config.w1);
                print(nameGraph,'-dpng')
                %  p = fig2plotly; needs more work, this looks horrible.
                status = -1;
                return;
            end
            [nameTopology nameVolFractionGrad]
            %-----------------------
            % Read the actual files
            % ---------------------
            DV.x = csvread(nameTopology);
            DV.w = csvread(nameVolFractionGrad);
            [config.nelx]   = size(  DV.x,2);
            [   config.nely] = size(  DV.x,1);
            
            FEACalls = i;
            obj.plotTopAndFraction(DV,  config, matProp, FEACalls); % plot the results.
            
            
            
        end
        
        
        %------------------------------------------------------------------
        % plot strain field after complete set of iterations.
        %------------------------------------------------------------------
        function  []= plotStrainField(obj,config,DV,folderNum,loadcaseIndex)
            
            ne = config.nelx*config.nely; % number of elements
            U2 = (DV.U(loadcaseIndex,:));
            maxU = max(max(abs(DV.U)));
            multiplierScale=1/maxU;
            
            for e = 1:ne
                %     if(DV.x>config.voidMaterialDensityCutOff
                % loop over local node numbers to get their node global node numbers
                % for
                j = 1:4;
                % Get the node number
                coordNodeNumber = DV.IEN(e,j);
                %   arrayCoordNumber(j) = coordNodeNumber;
                % get the global X,Y position of each node and put in array
                coord(j,:) = DV.globalPosition(coordNodeNumber,:);
                %     coord = coord+0.5; % because each node is a 1 by 1, so to make it line up with the fgm and top plot
                %  end
                arrayCoordNumber = coordNodeNumber;
                
                
                % plot the element outline
                %     if(doplotDisplacement ==1)
                hold on
                coordD = zeros(5,1);
                %           for
                temp = 1:4;
                %nodeNumber = IEN(e,j);
                %arrayCoordNumber(j) = coordNodeNumber(temp;
                coordD(temp,1) =  coord(temp,1) + transpose(multiplierScale*U2(2*arrayCoordNumber(temp)-1)); % X value
                coordD(temp,2) =  coord(temp,2) + transpose(multiplierScale*U2(2*arrayCoordNumber(temp))); % Y value
                %               coordD(temp,1) =  coord(temp,1)+ multiplierScale*U_displaced(arrayCoordNumber(temp),1); % X value
                %              coordD(temp,2) =  coord(temp,2)+ multiplierScale*U_displaced(arrayCoordNumber(temp),2); % Y value
                %           end
                
                
                coord2 = coord;
                coordD(5,:) = coordD(1,:) ;
                coord2(5,:) = coord2(1,:);
                
                coordD = coordD+0.5;
                coord2 = coord2+0.5;
                
                
                if(config.doUseMultiElePerDV ==1)
                    coord2(:,1)=coord2(:,1)/config.numXElmPerDV;
                    coordD(:,1)=coordD(:,1)/config.numXElmPerDV;
                    
                    coord2(:,2)=coord2(:,2)/config.numYElmPerDV;
                    coordD(:,2)=coordD(:,2)/config.numYElmPerDV;
                    
                    
                end
                
                plot(coord2(:,1),coord2(:,2),'-g');
                plot(coordD(:,1),coordD(:,2), '-b');
                %     end
            end
            hold off
            nameGraph = sprintf('Iter: %i, meso-macro iter %i, load case %i',folderNum,config.macro_meso_iteration,loadcaseIndex);
            title(nameGraph);
        end
        
        function  [] = PlotObjectiveFunctionAndConstraintsOverSevearlIterations(obj,maxiterations)
            
            close all
            
            numloops=maxiterations;
            
            folderNum = 0;
            all = [];
            iterationDiff = [];
            for i =1:numloops
                macro_meso_iteration = i;
                outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,macro_meso_iteration);
                storeOptimizationVar=  csvread(outname);
                all=[all;storeOptimizationVar];
                [rr,cc] = size (storeOptimizationVar);
                iterationDiff= [iterationDiff zeros(1,rr-1)];
                iterationDiff = [iterationDiff 1];
            end
            config = Configuration;
            DV = DesignVars(config);
            DV.storeOptimizationVar=all;
            p = plotResults;
            [r,c] = size (all);
            loopNumb=r;
            
            % Get size of the grid
            previousIterationNum=1;
            outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,previousIterationNum);
            [gridrows, gridcolumns]=size(csvread(outname));
            totalvolume = gridrows*gridcolumns;
            
            % normalize the volumes
            %mtemp = max(DV.storeOptimizationVar(1:loopNumb,4));
            % DV.storeOptimizationVar(1:loopNumb,4) = DV.storeOptimizationVar(1:loopNumb,4)/totalvolume;
            % %mtemp2 = max(DV.storeOptimizationVar(1:loopNumb,5));
            % DV.storeOptimizationVar(1:loopNumb,5) = DV.storeOptimizationVar(1:loopNumb,5)/totalvolume;
            
            
            
            config = 0;
            matProp = 0;
            hold on
            p.PlotDesignMetrics( DV, config, matProp, loopNumb)
            xvalues = 1:loopNumb;
            stairs(xvalues,iterationDiff);
            hold off
            
            w1 = 0;
            nameGraph = sprintf('./optiParaViaIterations%i.png', w1);
            print(nameGraph,'-dpng');
        end
        
        
        
    end
end