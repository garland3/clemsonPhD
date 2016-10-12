function    [U2, maxF,maxU] = AppliedStrainNoPeriodicNoTile(designVars, settings, matProp,macroElemProps)


doplotDisplacement = settings.doPlotAppliedStrain;
 multiplierScale = 5;
 
if doplotDisplacement ==1
    hi=  figure(2);
    cla(hi);
end


% ------------------------------------------------------
% For each displacement field case
% ------------------------------------------------------
[~, temp3] = size(settings.loadingCase);    


if( settings.averageMultiElementStrain==1)
    elementsPerDesignVar = settings.numXElmPerDV*settings.numYElmPerDV;
    temp3 = temp3*elementsPerDesignVar;
end
countElementsX = 1;
countElementsY = 1;
macroLoadCase = 1;

nn = (settings.nelx+1)*(settings.nely+1); % number of nodes
U2 = zeros(temp3,nn*2);

for loadcaseIndex = 1:temp3
     
    % ------------------------------------------------------
    % Single element per design var.
    % ------------------------------------------------------
    if(settings.doUseMultiElePerDV==0 )
        [Y,X] = ndgrid(0:1,0:1);  
        X = X*(  settings.nelx);
        Y = Y*(settings.nely);  

        macrodisplacementvector = macroElemProps.disp(loadcaseIndex,:);
        XD = [macrodisplacementvector(1)  macrodisplacementvector(3);
            macrodisplacementvector(7)  macrodisplacementvector(5)];

        YD = [macrodisplacementvector(2)  macrodisplacementvector(4);
            macrodisplacementvector(8)  macrodisplacementvector(6)];   

        Fx = griddedInterpolant(Y,X,XD,'linear');
        Fy = griddedInterpolant(Y,X,YD,'linear');
    elseif(settings.doUseMultiElePerDV==1 && settings.averageMultiElementStrain==1)
        % ------------------------------------------------------
        % Multiple element per design var. But, average the strain fields. 
        % ------------------------------------------------------
          [Y,X] = ndgrid(0:1,0:1);  
          X = X*(  settings.nelx);
          Y = Y*(settings.nely); 
          [t1,t2] = size( macroElemProps.mesoXnodelocations);
          
        %  macroLoadCase = floor((loadcaseIndex-1)/elementsPerDesignVar)+1;
        xd = macroElemProps.xDisplacements(macroLoadCase,:)';
        XD = reshape(xd,t2,t1)';

        yd = macroElemProps.yDisplacements(macroLoadCase,:)';
        YD = reshape(yd,t2,t1)';
        
        % then get the xd, yd for this element only
       elementsInRow = t1; % think this is actually "nodes in a row. "
        
        j = countElementsX;
        i = countElementsY;
         rowMultiplier = i-1;
         
         % not quite the same as IEN
        nodes = [rowMultiplier*elementsInRow+j, ...
                    rowMultiplier*elementsInRow+j+1, ...
%                     (rowMultiplier +1)*elementsInRow+j+1,...
                    (rowMultiplier +1)*elementsInRow+j,...
                      (rowMultiplier +1)*elementsInRow+j+1];
        
        [ macroLoadCase countElementsX countElementsY loadcaseIndex]
        countElementsX =  countElementsX +1;
        if(countElementsX>(t2-1)) % t2 is nodes, but we want elements, so minus 1
             countElementsY =  countElementsY + 1;
             countElementsX=1;
        end
        if(countElementsY>(t1-1)) % t2 is nodes, but we want elements, so minus 1
             countElementsY =  1; % reset to loading case. All elements finished. 
             countElementsX=1;
             macroLoadCase = macroLoadCase+1;
        end
   

        XD = XD(nodes);
        YD = YD(nodes);
        
       
        
        Fx = griddedInterpolant(Y,X,XD,'linear');
        Fy = griddedInterpolant(Y,X,YD,'linear');
        
        
    elseif(settings.doUseMultiElePerDV==1&& settings.averageMultiElementStrain~=1)
        % ------------------------------------------------------
        % Multiple element per design var.
        % ------------------------------------------------------
        X = macroElemProps.mesoXnodelocations;
        Y = macroElemProps.mesoYnodelocations;
        [t1,temp3] = size(X);
        X = X*(  settings.nelx)/(temp3-1);
        Y = Y*(settings.nely)/(t1-1); 

        xd = macroElemProps.xDisplacements(loadcaseIndex,:)';
        XD = reshape(xd,temp3,t1)';

        yd = macroElemProps.yDisplacements(loadcaseIndex,:)';
        YD = reshape(yd,temp3,t1)';


        Fx = griddedInterpolant(Y,X,XD,'linear');
        Fy = griddedInterpolant(Y,X,YD,'linear');

    end
   

    % Fx = scatteredInterpolant(displacementTargetX(:,1),displacementTargetX(:,2),displacementTargetX(:,3));
    % Fy = scatteredInterpolant(displacementTargetY(:,1),displacementTargetY(:,2),displacementTargetY(:,3));
    % [X,Y] = ndgrid(0:1,0:1);

    % Fx = griddedInterpolant(X,Y,XD,'linear');
    % Fy = griddedInterpolant(X,Y,YD,'linear');


    for i = 1:nn
        utemp = designVars.globalPosition(i,:);
        %     deltaX2 = Fx(utemp(1),utemp(2));
        %     deltaY2 = Fy(utemp(1),utemp(2));

        deltaX2 = Fx(utemp(2),utemp(1));
        deltaY2 = Fy(utemp(2),utemp(1));

        U2 (loadcaseIndex,2*i-1) = deltaX2; % X value
        U2 (loadcaseIndex,2*i) = deltaY2; % Y value

    end


    % loop over the elements
    if(doplotDisplacement ==1)
         subplot(1,temp3,loadcaseIndex)
        ne = settings.nelx*settings.nely; % number of elements
        for e = 1:ne
            % loop over local node numbers to get their node global node numbers
            % for
            j = 1:4;
            % Get the node number
            coordNodeNumber = designVars.IEN(e,j);
            %   arrayCoordNumber(j) = coordNodeNumber;
            % get the global X,Y position of each node and put in array
            coord(j,:) = designVars.globalPosition(coordNodeNumber,:);
            %  end
            arrayCoordNumber = coordNodeNumber;


            % plot the element outline
            if(doplotDisplacement ==1)
                hold on
                coordD = zeros(5,1);
                %           for
                temp = 1:4;
                %nodeNumber = IEN(e,j);
                %arrayCoordNumber(j) = coordNodeNumber(temp;
                coordD(temp,1) =  coord(temp,1) + transpose(multiplierScale*U2(loadcaseIndex,2*arrayCoordNumber(temp)-1)); % X value
                coordD(temp,2) =  coord(temp,2) + transpose(multiplierScale*U2(loadcaseIndex,2*arrayCoordNumber(temp))); % Y value
                %               coordD(temp,1) =  coord(temp,1)+ multiplierScale*U_displaced(arrayCoordNumber(temp),1); % X value
                %              coordD(temp,2) =  coord(temp,2)+ multiplierScale*U_displaced(arrayCoordNumber(temp),2); % Y value
                %           end

                coord2 = coord;
                coordD(5,:) = coordD(1,:) ;
                coord2(5,:) = coord2(1,:);
                %             plot(coord2(:,1),coord2(:,2),'-g');
                plot(coordD(:,1),coordD(:,2), '-b');
                 axis square
            end
            hold off
        end
    end

    plotoutline= doplotDisplacement;

    if(plotoutline==1 && settings.doUseMultiElePerDV==0)
        coord(:,1) = [0 1 1 0]*( settings.nelx);
        coord(:,2)  = [0 0 1 1]*( settings.nely);

        U3 = macroElemProps.disp(loadcaseIndex,:)*multiplierScale;
        coordD = zeros(5,2);
        for temp = 1:4
            coordD(temp,1) =  coord(temp,1)+ U3(2*temp-1); % X value
            coordD(temp,2) =  coord(temp,2)+ U3(2*temp); % Y value
        end
        coord2 = coord;
        coordD(5,:) = coordD(1,:) ;
        coord2(5,:) = coord2(1,:);
        %             coord3 = coord2;

        %     subplot(2,2,1);
        plot(coordD(:,1),coordD(:,2), '-r',coord2(:,1),coord2(:,2), '-g');
        %     axis([-0.3 1.3 -0.3 1.3])
        axis square

    end


    if(doplotDisplacement ==1)
        axis equal
        tti= strcat('Displacement of the elements shown in Blue');
        title(tti);
        hold off
        drawnow
    end




    maxF=1;
    maxU=1;
  
end


U2 = transpose(U2);