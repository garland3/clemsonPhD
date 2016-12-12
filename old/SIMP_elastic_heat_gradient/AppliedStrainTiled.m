function    [U2, maxF,maxU] = AppliedStrainTiled(designVars, settings, matProp,macroElemProps)

% first column is the X value,secnod is the Y value, 3rd is the
% displacement at tht particlar XY position in the X direction.
macrodisplacementvector = macroElemProps.disp;
displacementTargetX = [ 0 0 macrodisplacementvector(1);
    designVars.nelxTile 0 macrodisplacementvector(3);
    designVars.nelxTile designVars.nelyTile   macrodisplacementvector(5);
    0 designVars.nelyTile   macrodisplacementvector(7)];

displacementTargetY = [ 0 0 macrodisplacementvector(2);
    designVars.nelxTile 0 macrodisplacementvector(4);
    designVars.nelxTile designVars.nelyTile      macrodisplacementvector(6);
    0 designVars.nelyTile   macrodisplacementvector(8)];


multiplierScale = 10;
doplotDisplacement = 0;
if doplotDisplacement ==1
    %     figure(2)
end
nn = designVars.nelyTile *designVars.nelxTile ; % number of nodes = number elements since it wraps around

U2 = zeros(1,nn*2);

%U_displaced = zeros(nn,2); % delta X an Delta Y at each location
for i = 1:nn
    utemp = designVars.globalPositionTile(i,:); 
    
    deltaX = griddata(displacementTargetX(:,1),displacementTargetX(:,2), displacementTargetX(:,3),utemp(1),utemp(2));
    deltaY = griddata(displacementTargetY(:,1),displacementTargetY(:,2), displacementTargetY(:,3),utemp(1),utemp(2));
    
    U2 (2*i-1) = deltaX; % X value
    U2 (2*i) = deltaY; % Y value
    
end

ne =nn;

% ------------------------------
%    PLOT THE DISPLACEMENTS
%    (if plotting is on)
% ------------------------------
if doplotDisplacement ==1
    % loop over the elements
    for e = 1:ne
        % loop over local node numbers to get their node global node numbers
        % for
        j = 1:4;
        % Get the node number
        coordNodeNumber = designVars.IENTile(e,j);
        %   arrayCoordNumber(j) = coordNodeNumber;
        % get the global X,Y position of each node and put in array
        coord(j,:) = designVars.globalPositionTile(coordNodeNumber,:);
        %  end
        arrayCoordNumber = coordNodeNumber;
        
        % plot the element outline
        
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
        plot(coord2(:,1),coord2(:,2),'-g');
        plot(coordD(:,1),coordD(:,2), '-b');
    end
    
    
    axis equal
    tti= strcat('Displacement of the elements shown in Blue');
    title(tti);
    hold off
    % drawnow, uncommenting forced it to plot one at a time.
    hold off
end
m = 'done';

maxF=1;
maxU=1;
U2 = transpose(U2);