% ------------------------------------------------------------
%
% Tiles the meso repeating unit cell. Need the following things
%   1. tiled element locations
%   2. the node and element numbers need to be updated
%   3. the density and material volume fraciton need to be updated
%
% ------------------------------------------------------------
function designVarsMeso = TileMesoStructure(mesoSettings, designVarsMeso)



nelxTile = mesoSettings.nelx *mesoSettings.numTilesX;
nelyTile = mesoSettings.nely *mesoSettings.numTilesY ;
designVarsMeso.nelxTile=nelxTile;
designVarsMeso.nelyTile=nelyTile;

designVarsMeso.xTile(1:nelyTile,1:nelxTile) = mesoSettings.totalVolume; % artificial density of the elements
designVarsMeso.wTile(1:nelyTile,1:nelxTile)  = 1; % actual volume fraction composition of each element


designVarsMeso.XLocationsTile=zeros(nelyTile,nelxTile);
designVarsMeso.YLocationsTile=zeros(nelyTile,nelxTile);

designVarsMeso = designVarsMeso.CalcElementNodeMapmatrixWithPeriodicXandY_Tile(mesoSettings);
designVarsMeso =  designVarsMeso.CalcNodeLocationMeso(mesoSettings);


count = 1;
for i = 1:nelyTile  % y
    for j= 1:nelxTile % x
        %obj.globalPosition(count,:) = [j-1 i-1];
       
        %obj.XLocations(j,i) = j-1;
        %obj.YLocations(j,i) = i-1;
        %obj.
        microY = mod(i-1,mesoSettings.nely)+1;
        microX = mod(j-1,mesoSettings.nelx)+1;
        designVarsMeso.xTile(i,j) = designVarsMeso.x(microY,microX);
        designVarsMeso.wTile(i,j) = designVarsMeso.w(microY,microX);
        
        designVarsMeso.globalPositionTile(count,:) = [j-1 i-1];
        count = count+1;
        
    end
end
end