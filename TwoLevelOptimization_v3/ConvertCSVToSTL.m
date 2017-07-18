function [] = ConvertCSVToSTL(filename)

testing =1;

if(testing==1)
    completeStructure = csvread('./completeStucture1.000000_macroIteration_2.csv');
    disp('csv read');
elseif(testing==2)
    Z = peaks(20); % Create grid height
    completeStructure=Z;
    completeStructure(completeStructure>0)=1;
    completeStructure(completeStructure<0)=0;
else
    completeStructure = csvread(filename);
    disp('csv read');
end

[ysize,xsize]=size(completeStructure);

config = csvToStlConfig;



config.binary = 1;

config.targetXDimension = 220;
config.targetYDimension = 110;
config.targetZDimension = 20;

config.xscale=config.targetXDimension/xsize;
config.yscale=config.targetYDimension/ysize;
config.zscale=config.targetZDimension/1;

tricount=0;


config.mode = 2;




if(config.binary==0)
    config.outFile=fopen('structure3.stl','wt+');
    fprintf(config.outFile,'solid structure');
else
    config.mode = 1; % find number of triangles we will need
    tricount=LoopOverArray(completeStructure,config);
    
    config. outFile=fopen('structureBinary3.stl','wb');
    fprintf(config.outFile, '%-80s', 'structure');             % Title
    fwrite(config.outFile, tricount, 'uint32');           % Number of facets 1256
    tricount=0; % reset triangle count to 0
   config. mode = 2; % find number of triangles we will need
end

tricount=LoopOverArray(completeStructure,config); % write the data this time.


% Finish up the file format and close the file
if(config.binary==0)
    fprintf(config.outFile,'endsolid structure');
    fclose(config.outFile);
else
    fclose(config.outFile);
end
tricount



% yz, plane is for top and bottom edges
function tricount= addRect_xzPlane(x1,y1,x2Orginal,y2Original,z1,nx,ny,nz,config,tricount)

% point 1
x1_new=x1*config.xscale;
y1_new=y1*config.yscale;
z1_new = z1*config.zscale;

% point 2
x2=x2Orginal*config.xscale;
y2=y1*config.yscale;
z2 = z1_new;

% point 3
x3=(x1)*config.xscale;
y3=(y1)*config.yscale;
z3 = 0;

% point 4
x4=x2Orginal*config.xscale;
y4=(y1)*config.yscale;
z4 = 0;

if(ny>0)
    tricount=  addRectangle123_324(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount);
else
    tricount=   addRectangle321_423(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount);
end



% yz, plane is for left and right edges
function tricount= addRect_yzPlane(x1,y1,x2Orginal,y2Orginal,z1,nx,ny,nz,config,tricount)


% point 1
x1_new=x1*config.xscale;
y1_new=y1*config.yscale;
z1_new = z1*config.zscale;

% point 2
x2=(x1)*config.xscale;
y2=(y2Orginal)*config.yscale;
z2 = z1_new;

% point 3
x3=(x1)*config.xscale;
y3=(y1)*config.yscale;
z3 = 0;

% point 4
x4=(x1)*config.xscale;
y4=(y2Orginal)*config.yscale;
z4 = 0;

if(nx>0)
    tricount=addRectangle123_324(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount);
else
    tricount=addRectangle321_423(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount);
end


function tricount= addRect_XYPlane(x1,y1,x2orginal,y2orginal,z1,nx,ny,nz,config,tricount)


% point 1
x1_new=x1*config.xscale;
y1_new=y1*config.yscale;
z1_new = z1*config.zscale;

% point 2
x2=x2orginal*config.xscale;
y2=y1*config.yscale;
z2 = z1_new;

% point 3
x3=(x1)*config.xscale;
y3=y2orginal*config.yscale;
z3 = z1_new;

% point 4
x4=x2orginal*config.xscale;
y4=y2orginal*config.yscale;
z4 = z1_new;

if(nz>0)
    tricount= addRectangle123_324(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount);
else
    tricount=  addRectangle321_423(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount);
end


function [tricount]= addRectangle123_324(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount)

tricount=addTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz,config,tricount);
% addTriangle(x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz);
tricount=addTriangle(x3,y3,z3,x2,y2,z2,x4,y4,z4,nx,ny,nz,config,tricount);

function [tricount]=  addRectangle321_423(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz,config,tricount)

tricount=addTriangle(x3,y3,z3,x2,y2,z2,x1,y1,z1,nx,ny,nz,config,tricount);
% addTriangle(x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz);
tricount=addTriangle(x4,y4,z4,x2,y2,z2,x3,y3,z3,nx,ny,nz,config,tricount);







function[tricount]= addTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz,config,tricount)



if(config.mode ==2)
    % binary=1;
    if(config.binary ==1)
        fwrite(config.outFile,[nx ny nz],'float');
        fwrite(config.outFile,[x1 y1 z1],'float');
        fwrite(config.outFile,[x2 y2 z2],'float');
        fwrite(config.outFile,[x3 y3 z3],'float');
        fwrite(config.outFile,0,'ushort');
        %     fprintf(outFile,'  ');%uint16(0)
        
    else
        fprintf(config.outFile,'facet normal %f %f %f \n',nx,ny,nz);
        fprintf(config.outFile,'outer loop \n');
        fprintf(config.outFile,'vertex %f %f %f \n',x1,y1,z1);
        fprintf(config.outFile,'vertex %f %f %f \n',x2,y2,z2);
        fprintf(config.outFile,'vertex %f %f %f \n',x3,y3,z3);
        fprintf(config.outFile,'endloop \n');
        fprintf(config.outFile,'endfacet \n');
    end
end
tricount=tricount+1;

function [topedge,bottomedge,leftedge,rightedge,alreadyInArray,xlimit,ylimit,hasMaterial]=  GetEdgeInformationOfElement(x,y,completeStructure,alreadySearchedArray,config)
topedge = 0;
bottomedge = 0;
leftedge =0;
rightedge =0;
alreadyInArray=0;
hasMaterial=0;
% material here, add some top square, bottom square

% -------------------------------------
% check bottom
% -------------------------------------
i=x;
j=y;
xlimit=0;
ylimit=0;

if(x>config.xsize)
    xlimit=1;
    return;
end

if(y>config.ysize)
    ylimit=1;
    return;
end


alreadyInArray=alreadySearchedArray(j,i);
hasMaterial=completeStructure(j,i);

% check if on bottom edge
if(j==1)
    bottomedge=1;
else
    % check cell below's is empty
    value = completeStructure(j-1,i);
    if(value==0)
        bottomedge=1;
    end
end


% -------------------------------------
% check top
% -------------------------------------

% check if on top edge
if(j==config.ysize)
    topedge=1;
else
    % check cell below's is empty
    value = completeStructure(j+1,i);
    if(value==0)
        topedge=1;
    end
end


% -------------------------------------
% check left
% -------------------------------------

% check if on top edge
if(i==1)
    leftedge=1;
else
    % check cell left is empty
    value = completeStructure(j,i-1);
    if(value==0)
        leftedge=1;
    end
end

% -------------------------------------
% check Right
% -------------------------------------

% check if on top edge
if(i==config.xsize)
    rightedge=1;
else
    % check cell right is empty
    value = completeStructure(j,i+1);
    if(value==0)
        rightedge=1;
    end
end

% fprintf('Check Element x: %i y: %i Top: %i Bottom: %i Left: %i Right: %i AlreadyChecked: %i\n',x,y,topedge,bottomedge,leftedge,rightedge,alreadyInArray);

function tricount= LoopOverArray(completeStructure,config)

[ysize,xsize]=size(completeStructure);

config.xsize=xsize;
config.ysize = ysize;
total = ysize*xsize;
count =1;
tricount=0;

% search going right then up
% if element = 1 and alreadySearchedArray=0
% search right until, change in edge type or elem=0 or alreadySearchedArray=1
% -- then search up while, all eleemtns in row have same edge type and elem=1
% -- save the rectangle
% - record in the alreadySearchedArray all the elemements as 1

alreadySearchedArray=completeStructure*0;

updateOutputFrequency = 10000;
updatePlotFrequnecy = updateOutputFrequency*100;


% updateOutputFrequency = 1;
% updatePlotFrequnecy = updateOutputFrequency*1;

numOfColors=10;
colorByRegionArray=1:numOfColors;
regionCount=1;

for j = 1:ysize
    for i = 1:xsize
        
        alreadySerached = alreadySearchedArray(j,i);
        value = completeStructure(j,i);
        
        if(mod(count,updateOutputFrequency)==1)
            fprintf('\n\ncell %d  of %d, percent complete %0.01f with value %d alreadyCheck=%i tricount=%i\n',count,total,100*count/total,value,alreadySerached,tricount);
        end
        
        count=count+1;
        
        
        
        if(value==1 &&alreadySerached==0 )
            regionCount=regionCount+1;
            elementsRight=0;
            [topedge,bottomedge,leftedge,rightedge,alreadyInArray,xlimit,ylimit,hasMaterial]=  GetEdgeInformationOfElement(i,j,completeStructure,alreadySearchedArray,config);
            
            startRectX = i;
            startRectY=j;
            %             newtopedge=1;
            
            %             newtopedge=topedge;
            %             newbottomedge=
            
            elementsRight =elementsRight+1;
            [newtopedge,newbottomedge,newleftedge,newrightedge,newalreadyInArray,xlimit,ylimit,hasMaterial]=  GetEdgeInformationOfElement(i+elementsRight,j,completeStructure,alreadySearchedArray,config);
            
            % Search going right
            while(topedge==newtopedge && ...
                    bottomedge==newbottomedge &&...
                    leftedge== newleftedge && ...
                    rightedge== newrightedge&& ...
                    alreadyInArray== newalreadyInArray && ...
                    xlimit==0 && ...
                    ylimit==0  && ...
                    hasMaterial==1)
                
                
                elementsRight =elementsRight+1;
                [newtopedge,newbottomedge,newleftedge,newrightedge,newalreadyInArray,xlimit,ylimit,hasMaterial]=  GetEdgeInformationOfElement(i+elementsRight,j,completeStructure,alreadySearchedArray,config);
            end
            elementsRight=elementsRight-1; % the latest try was unsucessuful
            
            elementUP=0;
            % Search going Up
            rowOk = 1;
            while(rowOk==1)
                
                for ii = startRectX:startRectX+elementsRight
                    %                     fprintf('row %i',j+elementUP);
                    [newtopedge,newbottomedge,newleftedge,newrightedge,newalreadyInArray,xlimit,ylimit,hasMaterial]=  GetEdgeInformationOfElement(ii,j+elementUP,completeStructure,alreadySearchedArray,config);
                    if (topedge==newtopedge && ...
                            bottomedge==newbottomedge &&...
                            leftedge== newleftedge && ...
                            rightedge== newrightedge&& ...
                            alreadyInArray== newalreadyInArray && ...
                            xlimit==0 && ...
                            ylimit==0  && ...
                            hasMaterial==1)
                        % do nothing. All ok.
                        
                    else
                        rowOk=0;
                        
                        break
                    end
                    
                end
                if(rowOk==1)
                    elementUP =elementUP+1;
                end
                
                
            end
            elementUP =elementUP-1; % becaue the current row was not valid.
            
            endRectX = i+elementsRight;
            endRectY=j+elementUP;
            
            
            % Change the already seached array to show these are completed
            for jjj = startRectX:endRectX
                for iii=startRectY:endRectY
                    colorValue= 5+mod(regionCount,numOfColors);
                    alreadySearchedArray(iii,jjj)=colorValue;
                end
            end
            
            
            
            % -------------------------------------
            % Add material
            % -------------------------------------
            %
            %             tricount=   addRectOneUnit_XYPlane(i,j,1,0,0,1,config,tricount); % add top, up,1 unit
            %             tricount=   addRectOneUnit_XYPlane(i,j,0,0,0,-1,config,tricount); % add top, 0 plane, pointing down
            
            tricount=   addRect_XYPlane(startRectX,startRectY,endRectX+1,endRectY+1,1,0,0,1,config,tricount); % add top, up,1 unit
            tricount=   addRect_XYPlane(startRectX,startRectY,endRectX+1,endRectY+1,0,0,0,-1,config,tricount); % add top, 0 plane, pointing down
            
            if(topedge==1)
                
                tricount=    addRect_xzPlane(startRectX,endRectY+1,endRectX+1,'na',1,0,1,0,config,tricount); % add the top,plus 1 y direction
            end
            
            if(bottomedge==1)
                %                 tricount=    addRect_xzPlane(i,j,1,0,-1,0,config,tricount); % add the bottom, edge
                tricount=    addRect_xzPlane(startRectX,startRectY,endRectX+1,'na',1,0,-1,0,config,tricount); % add the top,plus 1 y direction
            end
            
            if(leftedge==1)
                %                 tricount=     addRect_yzPlane(i,j,1,1,0,0,config,tricount); % add left edge
                tricount=     addRect_yzPlane(startRectX,startRectY,'na',endRectY+1,1,1,0,0,config,tricount); % add left edge
            end
            
            if(rightedge==1)
                tricount=     addRect_yzPlane(endRectX+1,startRectY,'na',endRectY+1,1,-1,0,0,config,tricount); % add right edge , x plus 1
            end
            
           
            
        end % void regions
        
         if(mod(count,updatePlotFrequnecy)==0)
%                 fprintf('------------\nUpdating Plots\n ---------\n');
                p = plotResults;
                subplot(2,2,1)
                p.PlotArrayGeneric(alreadySearchedArray,'alreadySearchedArray');
                  map = colormap; % current colormap
                    map(1,:) = [1,  1,1];
                    colormap(map)
                subplot(2,2,2)
                p.PlotArrayGeneric(completeStructure,'completeStructure');
%                 subplot(2,2,3)
%                 p.PlotArrayGeneric(completeStructure-alreadySearchedArray,'diff');
                drawnow
            end
    end
end

p = plotResults;
subplot(1,1,1)
p.PlotArrayGeneric(alreadySearchedArray,'alreadySearchedArray');
map = colormap; % current colormap
map(1,:) = [1,  1,1];
colormap(map)
% subplot(2,2,2)
% p.PlotArrayGeneric(completeStructure,'completeStructure');
%                 subplot(2,2,3)
%                 p.PlotArrayGeneric(completeStructure-alreadySearchedArray,'diff');
drawnow
