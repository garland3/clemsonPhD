close all
coelhoMethod = 0;
DiffStore=[];
if(coelhoMethod==1)
    folderNum=0;
    for i = 1:20
        % get the density field
        outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,i);
        xnew = csvread(outname);
        
        if(i>1)
            DiffX = xold-xnew;
            sumOfDiffX = sum(sum(abs(DiffX)));
            
            %if(sumOfDiffX<=0.005)
            fprintf('%i,Diff = %f\n',i,sumOfDiffX)
            DiffStore=[DiffStore sumOfDiffX];
            %end
        end
        
        xold = xnew;
        
    end
    plot(DiffStore)
    title('Summed Difference in density values')
    xlabel('Iterations')
       ylabel('Summed Diff')
      %  nameGraph = sprintf('./gradTopOptimization%fwithmesh%i_load%i.png', config.w1,config.macro_meso_iteration,loadcaseIndex);
            print('SummedDiffPlotCoelhoMethod','-dpng');
end

matProp =MaterialProperties;
afterMode90runs = 0;
if(afterMode90runs==1)
    folderNum=0;
    for i = 1:5
        % get the density field
%         outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,i);
%         outname = sprintf('./out%i/ExxValues%i.csv',folderNum,i);
          outname = sprintf('./out%i/EyyValues%i.csv',folderNum,i);
        xnew = csvread(outname);
        
        if(i>1)
            DiffX = xold-xnew;
            sumOfDiffX = sum(sum(abs(DiffX)));
            sumOfDiffX=sumOfDiffX/matProp.E_material1;
            
            %if(sumOfDiffX<=0.005)
            fprintf('%i,Diff = %f\n',i,sumOfDiffX)
            DiffStore=[DiffStore sumOfDiffX];
            %end
        end
        
        xold = xnew;
        
        
        
    end
    plot(DiffStore)
    title('Summed Difference in density values')
    xlabel('Iterations')
   ylabel('Summed Diff')
   %  nameGraph = sprintf('./gradTopOptimization%fwithmesh%i_load%i.png', config.w1,config.macro_meso_iteration,loadcaseIndex);
    print('SummedDiffPlotCoelhoMethod','-dpng');
end

% ------------------------------
% Ground structure example plot
% ------------------------------
plotGroundstrExample=1;
if(plotGroundstrExample==1)
  
    
    [pointsx, pointsy] = meshgrid(0:3,0:3);
    pointsx=reshape(pointsx,1,[]);
      pointsy=reshape(pointsy,1,[]);
      
      
     c = ones(size(pointsy));
    sz = 100;
    scatter(pointsx,pointsy,sz,c,'filled')
    hold on
      
    
    xArrayToPlot = [];
     yArrayToPlot = [];
    
    nPoints = size(pointsx,2);
    for i =1:nPoints
        xhome = pointsx(i);
        yhome =pointsy(i);
        for j =1:nPoints
                newPointX = pointsx(j);
                 newPointY =pointsy(j);
%                  if(xhome~=newPointX &&yhome~=newPointY)
                     xArrayToPlot=[xArrayToPlot xhome];
                     yArrayToPlot=[yArrayToPlot yhome];
                     
                     xArrayToPlot=[xArrayToPlot newPointX];
                     yArrayToPlot=[yArrayToPlot newPointY];
                     
%                  end
        
        end
        
    end
%     
%     xlines=[ xhome 1 xhome 1 0 2 0 2]
%     ylines=[ 0 1 0 2 0 1 0 2]
    plot(xArrayToPlot, yArrayToPlot)
    
%      % node 1
%     xlines=[ 1 0 1 1 1 2]
%     ylines=[ 0 1 0 1 0 1]
%     plot(xlines, ylines)
%     
%       % node 2
%     xlines=[ 2 0 1 1 1 2]
%     ylines=[ 0 1 0 1 0 1]
%     plot(xlines, ylines)
    
    hold off
end