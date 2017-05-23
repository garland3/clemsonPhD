% ------------------------------------------------------------
%
% Tiles the meso repeating unit cell. Need the following things
%   1. tiled element locations
%   2. the node and element numbers need to be updated
%   3. the density and material volume fraciton need to be updated
%
% ------------------------------------------------------------
function completeStruct = TileMesoStructureV2(mesoSettings, macroSettings,DVMeso,macroElementProps,densityFieldMacro,completeStruct, step)


xx=densityFieldMacro;

% DVMeso.nelxTile=nelxTile;
% DVMeso.nelyTile=nelyTile;


numTilesX=mesoSettings.numTilesX;
numTilesY = mesoSettings.numTilesY;




% preform the logic tests  to see if we need to add material or remove
% material for the post processing.
cut_topLeft = 0;
cut_topRight = 0;
cut_bottomLeft = 0;
cut_bottomRight = 0;

add_top=0;
add_bottom=0;
add_left=0;
add_right=0;

xCurrent=macroElementProps.xPos;
xRight=xCurrent+1;
xLeft=xCurrent-1;

yCurrent=macroElementProps.yPos;
yUp=yCurrent+1;
yDown=yCurrent-1;

 % densityCurrentElement = xx(yCurrent, xCurrent);

% densityUp = xx(yUp, xCurrent);
% densityDown = xx(yDown, xCurrent);
%
% densityxLeft= xx(yCurrent, xLeft);
% densityxRight= xx(yCurrent, xRight);

if(macroSettings.multiscaleMethodCompare==1)
    nelxTile = mesoSettings.nelx *(mesoSettings.numTilesX+add_left+add_right );
    nelyTile = mesoSettings.nely *(mesoSettings.numTilesY+add_top+add_bottom );
    yShift = (macroElementProps.yPos-1)*mesoSettings.nely*numTilesY+1-(mesoSettings.nely*add_bottom);
    xShift = (macroElementProps.xPos-1)*mesoSettings.nelx*numTilesX+1-(mesoSettings.nelx*add_left);
    
    % --------------------------------------------------
    % Generate the actual tiled matrix
    % --------------------------------------------------
    tiledDensity = zeros(nelyTile,nelxTile);
    % ZeroRegions = ones(nelyTile,nelxTile);
    
    % count = 1;
    for i = 1:nelyTile  % y
        for j= 1:nelxTile % x
            
            microY = mod(i-1,mesoSettings.nely)+1;
            microX = mod(j-1,mesoSettings.nelx)+1;
            tiledDensity(i,j) = DVMeso.x(microY,microX);
        end
    end
    
    
    ybounds = yShift:(yShift+nelyTile-1);
    xbounds = xShift:(xShift+nelxTile-1);
    completeStruct(ybounds,xbounds)= completeStruct(ybounds,xbounds)+tiledDensity;
    return;
end


% -----------------------------
% Check Right
% -------------------
if(xRight<=macroSettings.nelx)
    
    density = xx(yCurrent, xRight);
    if(density>macroSettings.voidMaterialDensityCutOff)
        add_right=1;
    end
end


% -----------------------------
% Check Left
% -------------------
if(xLeft>0)
    
    density = xx(yCurrent, xLeft);
    if(density>macroSettings.voidMaterialDensityCutOff)
        add_left=1;
    end
end


% -----------------------------
% Check down
% -------------------
if(yDown>0)
    
    density = xx(yDown, xCurrent);
    if(density>macroSettings.voidMaterialDensityCutOff)
        add_bottom=1;
    else
        % Check Right, then Check left; checking for empty
        if(xRight<macroSettings.nelx)
            densityxRight= xx(yCurrent, xRight);
            if(densityxRight<macroSettings.voidMaterialDensityCutOff)
                cut_bottomRight=1;
            end
        end
        
        if(xLeft>0)
            densityxLeft= xx(yCurrent, xLeft);
            if(densityxLeft<macroSettings.voidMaterialDensityCutOff)
                cut_bottomLeft=1;
            end
        end
    end
end


% -----------------------------
% Check up
% -------------------
if(yUp<=macroSettings.nely)
    
    densityUp = xx(yUp, xCurrent);
    if(densityUp>macroSettings.voidMaterialDensityCutOff)
        add_top=1;
    else
        % Check Right, then Check left; checking for empty
        if(xRight<macroSettings.nelx)
            densityxRight= xx(yCurrent, xRight);
            if(densityxRight<macroSettings.voidMaterialDensityCutOff)
                cut_topRight=1;
            end
        end
        
        if(xLeft>0)
            densityxLeft= xx(yCurrent, xLeft);
            if(densityxLeft<macroSettings.voidMaterialDensityCutOff)
                cut_topLeft=1;
            end
        end
    end
end

add_leftV2 = add_left;
add_topV2 = add_top;
add_rightV2 = add_right;
add_bottomV2 = add_bottom;

% changed my mind, only do bottom and right overlaps
add_left=0;
add_top=0;

% changed my mind again. Do none, only do bottom and right overlaps
add_right=0;
add_bottom=0;

nelxTile = mesoSettings.nelx *(mesoSettings.numTilesX+add_left+add_right );
nelyTile = mesoSettings.nely *(mesoSettings.numTilesY+add_top+add_bottom );
yShift = (macroElementProps.yPos-1)*mesoSettings.nely*numTilesY+1-(mesoSettings.nely*add_bottom);
xShift = (macroElementProps.xPos-1)*mesoSettings.nelx*numTilesX+1-(mesoSettings.nelx*add_left);


ybounds = yShift:(yShift+nelyTile-1);
xbounds = xShift:(xShift+nelxTile-1);

if(step ==1)
    
    
    % --------------------------------------------------
    % Generate the actual tiled matrix
    % --------------------------------------------------
    tiledDensity = zeros(nelyTile,nelxTile);
    % ZeroRegions = ones(nelyTile,nelxTile);
    
    % count = 1;
    for i = 1:nelyTile  % y
        for j= 1:nelxTile % x
            
            microY = mod(i-1,mesoSettings.nely)+1;
            microX = mod(j-1,mesoSettings.nelx)+1;
            tiledDensity(i,j) = DVMeso.x(microY,microX);
        end
    end
    
    
    % --------------------------------------------------------------
    %
    % Add to existing complete structure!
    %
    % Add the design, rather than replace it.
    % --------------------------------------------------------------
    completeStruct(ybounds,xbounds)= completeStruct(ybounds,xbounds)+tiledDensity;
    
    
elseif(step==2)
    
    % --------------------------------------------------------------
    % Add the Cuts!
    % --------------------------------------------------------------
%     if(cut_topRight ==1 || cut_topLeft ==1 || cut_bottomLeft==1  || cut_bottomRight==1 )
        temp = completeStruct(ybounds,xbounds);
        %tiledDensity(flipud(triu(tiledDensity)==tiledDensity)) = 0;
        if(cut_topRight ==1)
            temp(flipud(triu(temp)==temp)) = 0;
        end
        
        if(cut_topLeft ==1)
            temp(tril(temp)==temp) = 0;
        end
        
        if(cut_bottomRight ==1)
            temp(triu(temp)==temp) = 0;
        end
        
        if(cut_bottomLeft ==1)
            temp(flipud(tril(temp)==temp)) = 0;
        end
        
        if(macroElementProps.densitySIMP<macroSettings.voidMaterialDensityCutOff)
            temp=zeros(size(temp));
        end
        
        
        % ----------------------------------------
        % Add external Boundary
        % ----------------------------------------
        numrows =9;
        % -----------------------------------
        if(macroElementProps.densitySIMP>macroSettings.voidMaterialDensityCutOff)
            % if solid, and on the edge of design domain, then add material
            % -----------------------------------
            if(yCurrent==1) % then on the bottom row
                temp(1:numrows,:)=1;
            end
            
            if(yCurrent==macroSettings.nely  ) % then on the bottom row
                temp(end-numrows:end,:)=1;
            end
            
            if(xCurrent==1 ) % then on the bottom row
                temp(:,1:numrows)=1;
            end
            
            if(xCurrent==macroSettings.nelx) % then on the bottom row
                temp(:,end-numrows:end)=1;
            end
            
            % ----------------------------------------------
            % if  solid, and AddXXX is false then add material in the XX
            % direction
            % ----------------------------------
            if( add_topV2==0 && (cut_topLeft==0 && cut_topRight ==0))
                 temp(end-numrows:end,:)=1;
            end
            
            if( add_bottomV2==0 && (cut_bottomLeft==0 && cut_bottomRight ==0))
                temp(1:numrows,:)=1;
            end
            
            if( add_leftV2==0 && (cut_topLeft==0 && cut_bottomLeft ==0))
                temp(:,1:numrows)=1;
            end
            
            if( add_rightV2==0 && (cut_topRight==0 && cut_bottomRight ==0))
                temp(:,end-numrows:end)=1;
            end
            
            % ----------------------------
            % Add diagianols
            % ----------------------------
            [t1 t2]=size(temp);
            tempOnes = ones(size(temp));
            
           
            if(cut_topLeft ==1)
                 for i = 0:numrows-1
                    DiagMiddle =diag(ones(1,t1-i),i);
                    temp(DiagMiddle==1) = tempOnes(DiagMiddle==1);
                 end
            end
            
            if(cut_topRight ==1)
                for i = 0:numrows-1
                    DiagMiddle =fliplr(diag(ones(1,t1-i),i));
                    temp(DiagMiddle==1) = tempOnes(DiagMiddle==1);
                 end
            end
            
             if(cut_bottomRight ==1)
                for i = 0:numrows-1
                    DiagMiddle =rot90(diag(ones(1,t1-i),i),2);
                    temp(DiagMiddle==1) = tempOnes(DiagMiddle==1);
                 end
             end
            
            if(cut_bottomLeft ==1)
                for i = 0:numrows-1
                    DiagMiddle =flipud(diag(ones(1,t1-i),i));
                    temp(DiagMiddle==1) = tempOnes(DiagMiddle==1);
                 end
            end
            

            
        end
        
        completeStruct(ybounds,xbounds)= temp;
        %     end


    
end

end