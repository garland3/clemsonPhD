function    [U2, maxF,maxU] = AppliedStrainNoPeriodicNoTile(designVars, settings, matProp,macroElemProps)

% first column is the X value,secnod is the Y value, 3rd is the
% displacement at tht particlar XY position in the X direction.
macrodisplacementvector = macroElemProps.disp;
% displacementTargetX = [ 0 0 macrodisplacementvector(1);
%                         settings.nelx+1 0 macrodisplacementvector(3);
%                         settings.nelx+1 settings.nely+1    macrodisplacementvector(5);
%                         0 settings.nely+1   macrodisplacementvector(7)];
% 
% displacementTargetY = [ 0 0 macrodisplacementvector(2);
%                         settings.nelx+1 0 macrodisplacementvector(4);
%                         settings.nelx+1 settings.nely+1     macrodisplacementvector(6);
%                         0 settings.nely+1  macrodisplacementvector(8)];

 [X,Y] = ndgrid(0:1,0:1);
 
%  X(2,:) =   settings.nelx+1;
%  Y(:,2) = settings.nely+1;
 
 X = X*(  settings.nelx+1);
 Y = Y*(settings.nely+1);
 
 XD = [macrodisplacementvector(1)  macrodisplacementvector(7);
      macrodisplacementvector(3)  macrodisplacementvector(5)];
  
 YD = [macrodisplacementvector(2)  macrodisplacementvector(8);
      macrodisplacementvector(4)  macrodisplacementvector(6)];
            

multiplierScale = 50;
doplotDisplacement = settings.doPlotAppliedStrain;
if doplotDisplacement ==1
        hi=  figure(2);
        cla(hi);
end
nn = (settings.nelx+1)*(settings.nely+1); % number of nodes

U2 = zeros(1,nn*2);

% Fx = scatteredInterpolant(displacementTargetX(:,1),displacementTargetX(:,2),displacementTargetX(:,3));
% Fy = scatteredInterpolant(displacementTargetY(:,1),displacementTargetY(:,2),displacementTargetY(:,3));
% [X,Y] = ndgrid(0:1,0:1);

Fx = griddedInterpolant(X,Y,XD,'linear');
Fy = griddedInterpolant(X,Y,YD,'linear');
for i = 1:nn
    utemp = designVars.globalPosition(i,:);    
    deltaX2 = Fx(utemp(1),utemp(2));
    deltaY2 = Fy(utemp(1),utemp(2));    
    
    U2 (2*i-1) = deltaX2; % X value
    U2 (2*i) = deltaY2; % Y value
    
end


% loop over the elements
if(doplotDisplacement ==1)
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
            coordD(temp,1) =  coord(temp,1) + transpose(multiplierScale*U2(2*arrayCoordNumber(temp)-1)); % X value
            coordD(temp,2) =  coord(temp,2) + transpose(multiplierScale*U2(2*arrayCoordNumber(temp))); % Y value
            %               coordD(temp,1) =  coord(temp,1)+ multiplierScale*U_displaced(arrayCoordNumber(temp),1); % X value
            %              coordD(temp,2) =  coord(temp,2)+ multiplierScale*U_displaced(arrayCoordNumber(temp),2); % Y value
            %           end

            coord2 = coord;
            coordD(5,:) = coordD(1,:) ;
            coord2(5,:) = coord2(1,:);
%             plot(coord2(:,1),coord2(:,2),'-g');
            plot(coordD(:,1),coordD(:,2), '-b');
        end
    end
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
U2 = transpose(U2);