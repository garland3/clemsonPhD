function [] = ConvertCSVToSTL_old(filename)

testing =1;

if(testing~=1)
    completeStructure = csvread('./completeStucture1.000000_macroIteration_1.csv');
    disp('csv read');
else    
    Z = peaks(20); % Create grid height
    completeStructure=Z;
    completeStructure(completeStructure>0)=1;
    completeStructure(completeStructure<0)=0;   
end

[ysize,xsize]=size(completeStructure);



global targetXDimension;
global targetYDimension;% = 100 ; % 100 mm
global targetZDimension; % = 100 ; % 100 mm
global xscale;
global yscale;
global zscale;
global outFile;
global tricount;

global binary;
binary = 1;

targetXDimension = 200;
targetYDimension = 200;
targetZDimension = 20;

xscale=targetXDimension/xsize;
yscale=targetYDimension/ysize;
zscale=targetZDimension/1;

tricount=0;

global mode;
mode = 2;




if(binary~=1)
    outFile=fopen('structure.stl','wt+');
    fprintf(outFile,'solid structure');
else
    mode = 1; % find number of triangles we will need
    LoopOverArray(completeStructure)
    outFile=fopen('structure.stl','wb');    
    fprintf(outFile, '%-80s', 'structure');             % Title
    fwrite(outFile, tricount, 'uint32');           % Number of facets 1256
    tricount=0; % reset triangle count to 0
    mode = 2; % find number of triangles we will need
end

  LoopOverArray(completeStructure); % write the data this time. 


% Finish up the file format and close the file
if(binary~=1)
    fprintf(outFile,'endsolid structure');
    fclose(outFile);
else
    fclose(outFile);   
end
tricount



% yz, plane is for top and bottom edges
function addRectOneUnit_xzPlane(x1,y1,z1,nx,ny,nz)
global xscale;
global yscale;
global zscale;

% point 1
x1_new=x1*xscale;
y1_new=y1*yscale;
z1_new = z1*zscale;

% point 2
x2=(x1+1)*xscale;
y2=y1*yscale;
z2 = z1_new;

% point 3
x3=(x1)*xscale;
y3=(y1)*yscale;
z3 = 0;

% point 4
x4=(x1+1)*xscale;
y4=(y1)*yscale;
z4 = 0;

if(ny>0)
    addRectangle123_324(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)
else
    addRectangle321_423(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)
end



% yz, plane is for left and right edges
function addRectOneUnit_yzPlane(x1,y1,z1,nx,ny,nz)
global xscale;
global yscale;
global zscale;

% point 1
x1_new=x1*xscale;
y1_new=y1*yscale;
z1_new = z1*zscale;

% point 2
x2=(x1)*xscale;
y2=(y1+1)*yscale;
z2 = z1_new;

% point 3
x3=(x1)*xscale;
y3=(y1)*yscale;
z3 = 0;

% point 4
x4=(x1)*xscale;
y4=(y1+1)*yscale;
z4 = 0;

if(nx>0)
    addRectangle123_324(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)
else
    addRectangle321_423(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)
end


function addRectOneUnit_XYPlane(x1,y1,z1,nx,ny,nz)
global xscale;
global yscale;
global zscale;

% point 1
x1_new=x1*xscale;
y1_new=y1*yscale;
z1_new = z1*zscale;

% point 2
x2=(x1+1)*xscale;
y2=y1*yscale;
z2 = z1_new;

% point 3
x3=(x1)*xscale;
y3=(y1+1)*yscale;
z3 = z1_new;

% point 4
x4=(x1+1)*xscale;
y4=(y1+1)*yscale;
z4 = z1_new;

if(nz>0)
    addRectangle123_324(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)
else
    addRectangle321_423(x1_new,y1_new,z1_new,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)
end


function addRectangle123_324(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)

addTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz);
% addTriangle(x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz);
addTriangle(x3,y3,z3,x2,y2,z2,x4,y4,z4,nx,ny,nz);

function addRectangle321_423(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz)

addTriangle(x3,y3,z3,x2,y2,z2,x1,y1,z1,nx,ny,nz);
% addTriangle(x2,y2,z2,x3,y3,z3,x4,y4,z4,nx,ny,nz);
addTriangle(x4,y4,z4,x2,y2,z2,x3,y3,z3,nx,ny,nz);




function addTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,nx,ny,nz)

global outFile;
global tricount;
global binary;
global mode;

if(mode ==2)
    % binary=1;
    if(binary ==1)
        fwrite(outFile,[nx ny nz],'float');
        fwrite(outFile,[x1 y1 z1],'float');
        fwrite(outFile,[x2 y2 z2],'float');
        fwrite(outFile,[x3 y3 z3],'float');
        fwrite(outFile,0,'ushort');
        %     fprintf(outFile,'  ');%uint16(0)

    else
        fprintf(outFile,'facet normal %f %f %f \n',nx,ny,nz);
        fprintf(outFile,'outer loop \n');
        fprintf(outFile,'vertex %f %f %f \n',x1,y1,z1);
        fprintf(outFile,'vertex %f %f %f \n',x2,y2,z2);
        fprintf(outFile,'vertex %f %f %f \n',x3,y3,z3);
        fprintf(outFile,'endloop \n');
        fprintf(outFile,'endfacet \n');
    end
end
tricount=tricount+1;


function LoopOverArray(completeStructure)

[ysize,xsize]=size(completeStructure);


total = ysize*xsize;
count =1;


for j = 1:ysize
    for i = 1:xsize
        
        value = completeStructure(j,i);
        
        if(mod(count,1000)==1)
            fprintf('cell %d  of %d, percent complete %0.01f with value %d\n',count,total,100*count/total,value);
        end
        
        count=count+1;
        
        topedge = 0;
        bottomedge = 0;
        leftedge =0;
        rightedge =0;
        
        if(value==1)
            % material here, add some top square, bottom square
            
            % -------------------------------------
            % check bottom
            % -------------------------------------
            
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
            if(j==ysize)
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
            if(i==xsize)
                rightedge=1;
            else
                % check cell right is empty
                value = completeStructure(j,i+1);
                if(value==0)
                    rightedge=1;
                end
            end
            
            
            % -------------------------------------
            % Add material
            % -------------------------------------
            addRectOneUnit_XYPlane(i,j,1,0,0,1); % add top, up,1 unit
            addRectOneUnit_XYPlane(i,j,0,0,0,-1); % add top, 0 plane, pointing down
            
            if(topedge==1)
                addRectOneUnit_xzPlane(i,j+1,1,0,1,0); % add the top,plus 1 y direction
            end
            
            if(bottomedge==1)
                addRectOneUnit_xzPlane(i,j,1,0,-1,0); % add the bottom, edge
            end
            
            if(leftedge==1)
                addRectOneUnit_yzPlane(i,j,1,1,0,0); % add left edge
            end
            
            if(rightedge==1)
                addRectOneUnit_yzPlane(i+1,j,1,-1,0,0); % add right edge , x plus 1
            end
            
            
            
        end % void regions
    end
end