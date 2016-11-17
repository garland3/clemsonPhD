function  []= plotStrainField(settings,designVars,folderNum,loadcaseIndex)

ne = settings.nelx*settings.nely; % number of elements
U2 = (designVars.U(loadcaseIndex,:));
 maxU = max(max(abs(designVars.U)));
multiplierScale=1/maxU;

for e = 1:ne
%     if(designVars.x>settings.voidMaterialDensityCutOff
    % loop over local node numbers to get their node global node numbers
    % for
    j = 1:4;
    % Get the node number
    coordNodeNumber = designVars.IEN(e,j);
    %   arrayCoordNumber(j) = coordNodeNumber;
    % get the global X,Y position of each node and put in array
    coord(j,:) = designVars.globalPosition(coordNodeNumber,:);
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
        
        
        if(settings.doUseMultiElePerDV ==1)
            coord2(:,1)=coord2(:,1)/settings.numXElmPerDV;
            coordD(:,1)=coordD(:,1)/settings.numXElmPerDV;
            
            coord2(:,2)=coord2(:,2)/settings.numYElmPerDV;
            coordD(:,2)=coordD(:,2)/settings.numYElmPerDV;
            
            
        end
        
        plot(coord2(:,1),coord2(:,2),'-g');
        plot(coordD(:,1),coordD(:,2), '-b');
%     end
end
hold off
  nameGraph = sprintf('Iter: %i, meso-macro iter %i, load case %i',folderNum,settings.macro_meso_iteration,loadcaseIndex);
       title(nameGraph);

%   nameGraph = sprintf('./out%i/strainGraph%i.png',folderNum,macro_meso_iteration);
% print(nameGraph,'-dpng')