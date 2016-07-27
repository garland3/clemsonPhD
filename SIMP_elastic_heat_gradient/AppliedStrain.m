function    [U2, maxF,maxU] = AppliedStrain(designVars, settings, matProp,macroElemProps)

% count = 1;
%     elementsInRow = settings.nelx+1;
%     nn = (settings.nelx+1)*(settings.nely+1); % number of nodes
%     U = zeros(nn,4);
%     % Each row, so nely # of row
%     for i = 1:settings.nely
%         rowMultiplier = i-1;
%         % Each column, so nelx # of row
%         for j= 1:settings.nelx
%             U(count,:)=[rowMultiplier*elementsInRow+j, ...
%                 rowMultiplier*elementsInRow+j+1, ...
%                 (rowMultiplier +1)*elementsInRow+j+1,...
%                 (rowMultiplier +1)*elementsInRow+j];
%             count = count+1;
%         end
%     end

% displacementTargetX = [Xcoord, Ycoord, Xdisplace] %each row is a new sample point
% displacementTargetY = [Xcoord, Ycoord, Ydisplace] %each row is a new sample point

% first column is the X value,secnod is the Y value, 3rd is the
% displacement at tht particlar XY position in the X direction. 
macrodisplacementvector = macroElemProps.disp;
displacementTargetX = [ 0 0 macrodisplacementvector(1);
                        settings.nelx 0 macrodisplacementvector(3);
                        0 settings.nely macrodisplacementvector(5);
                        settings.nelx settings.nely macrodisplacementvector(7)];
                    
displacementTargetY = [ 0 0 macrodisplacementvector(2);
                        settings.nelx 0 macrodisplacementvector(4);
                        0 settings.nely macrodisplacementvector(6);
                        settings.nelx settings.nely macrodisplacementvector(8)];

% diplacementTopRightX = 0.6;
% diplacementTopRightY = 0.6;
multiplierScale = 1;


doplotDisplacement = 1;
if doplotDisplacement ==1
    figure(2)
end
nn = (settings.nelx)*(settings.nely); % number of nodes

% 
% slopeX = diplacementTopRightX/settings.nelx;
% slopeY = diplacementTopRightY/settings.nely;
U2 = zeros(1,nn*2);

%U_displaced = zeros(nn,2); % delta X an Delta Y at each location
for i = 1:nn
    utemp = designVars.globalPosition(i,:);
    
%     deltaX = slopeX*utemp(1);
%     deltaY = slopeY*utemp(2);
    deltaX = griddata(displacementTargetX(:,1),displacementTargetX(:,2), displacementTargetX(:,3),utemp(1),utemp(2));
    deltaY = griddata(displacementTargetY(:,1),displacementTargetY(:,2), displacementTargetY(:,3),utemp(1),utemp(2));
    
%     newX = utemp(1)+deltaX;
%     newY = utemp(2) + deltaY;

    U2 (2*i-1) = deltaX; % X value
    U2 (2*i) = deltaY; % Y value
    
    
   % U_displaced(i,:) = [deltaX deltaY];
end
%arrayCoordNumber= zeros(1,4);
ne = settings.nelx*settings.nely; % number of elements

% loop over the elements
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
          for temp = 1:4
             %nodeNumber = IEN(e,j);
             %arrayCoordNumber(j) = coordNodeNumber(temp;
             coordD(temp,1) =  coord(temp,1)+ multiplierScale*U2(2*arrayCoordNumber(temp)-1); % X value
             coordD(temp,2) =  coord(temp,2)+ multiplierScale*U2(2*arrayCoordNumber(temp)); % Y value
%               coordD(temp,1) =  coord(temp,1)+ multiplierScale*U_displaced(arrayCoordNumber(temp),1); % X value
%              coordD(temp,2) =  coord(temp,2)+ multiplierScale*U_displaced(arrayCoordNumber(temp),2); % Y value
          end    

         coord2 = coord;
          coordD(5,:) = coordD(1,:) ;
         coord2(5,:) = coord2(1,:); 
        % plot(coord2(:,1),coord2(:,2),'-g');  
          plot(coordD(:,1),coordD(:,2), '-b');   
       end    
end

if(doplotDisplacement ==1)
    axis equal
    tti= strcat('Displacement of the elements shown in Blue');
    title(tti);
    hold off 
      drawnow
end

hold off
m = 'done';

maxF=1;
maxU=1;
U2 = transpose(U2);