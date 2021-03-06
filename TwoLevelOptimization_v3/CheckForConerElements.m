function [xnew, numChanged] = CheckForConerElements(x, nelx, nely, cutoff)
% config = mesoConfig;


% Get the density field
x(x>cutoff)=1;
x(x<cutoff)=0.001;
xnew = x;
% xold=x;

 p = plotResults;
 masterDoPlot = 0;

beforeSumGlobal=sum(sum(x));
densityChangedLocal=1;
loopcount=1;
while densityChangedLocal~=0
    xold=xnew;
    localSumBefore=sum(sum(xnew));
    % preform the logic tests  to see if we need to add material
    % has_top=0;
    % has_top_right=0;
    % has_top_left=0;
    % has_bottom_right=0;
    % has_bottom_left=0;
    % has_bottom=0;
    % has_left=0;
    % has_right=0;
    
    for elx = 1:nelx
        for ely = 1:nely
            
            xCurrentSimp=xold(ely,elx);
            
            doPlot = 0;
            xCurrent=elx;
            xRight=xCurrent+1;
            xLeft=xCurrent-1;
            
            yCurrent=ely;
            yUp=yCurrent+1;
            yDown=yCurrent-1;
            
            has_top=0;
            has_top_right=0;
            has_top_left=0;
            has_bottom_right=0;
            has_bottom_left=0;
            has_bottom=0;
            has_left=0;
            has_right=0;
            type='';
            
            % -----------------------------
            % Check up
            % -------------------
            if(yUp<=nely)
                densityUp = xold(yUp, xCurrent);
                if(densityUp>cutoff)
                    has_top=1;
                end
            end
            
            
            % -----------------------------
            % Check Right
            % -------------------
            if(xRight<=nelx)
                
                densityRight = xold(yCurrent, xRight);
                if(densityRight>cutoff)
                    has_right=1;
                end
            end
            
            % -----------------------------
            % Check Right TOP
            % -------------------
            if(xRight<=nelx &&  yUp<=nely)
                
                densityTopRight = xold(yUp, xRight);
                if(densityTopRight>cutoff)
                    has_top_right=1;
                end
            end
            
            
            % -----------------------------
            % Check  Left TOP
            % -------------------
            if(xLeft>0 &&  yUp<=nely)
                
                densityTopLeft = xold(yUp, xLeft);
                if(densityTopLeft>cutoff)
                    has_top_left=1;
                end
            end
            
            
            % -----------------------------
            % Check Left
            % -------------------
            if(xLeft>0)                
                densityLeft = xold(yCurrent, xLeft);
                if(densityLeft>cutoff)
                    has_left=1;
                end
            end
            
            % -----------------------------
            % Check Left Bottom
            % -------------------
            if(xLeft>0 &&  yDown>0)                
                densityLeftBottom = xold(yDown, xLeft);
                if(densityLeftBottom>cutoff)
                    has_bottom_left=1;
                end
            end
            
            
            % -----------------------------
            % Check  Right Bottom
            % -------------------
            if(xRight<=nelx &&   yDown>0)
                
                densityRightBottom = xold(yDown, xRight);
                if(densityRightBottom>cutoff)
                    has_bottom_right=1;
                end
            end
            
            
            % -----------------------------
            % Check down
            % -------------------
            if(yDown>0)
                
                densityDown = xold(yDown, xCurrent);
                if(densityDown>cutoff)
                    has_bottom=1;
                end
            end
            
            
            valueForChangedLocations=0.5;
            if(xCurrentSimp>cutoff)
                % ---------------------
                %
                %         ADD MATERAIL
                %
                %
                % if has a corner, but no top or side, then add material.
                
                if(has_top_right==1)
                    if(has_top==0 && has_right==0)
                        xnew(yUp, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xRight)=valueForChangedLocations;
                        doPlot=1;
                        type='top right';
                    end
                end
                
                if(has_top_left==1)
                    if(has_top==0 && has_left==0)
                        xnew(yUp, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xLeft)=valueForChangedLocations;
                         doPlot=1;
                          type='top left';
                    end
                end
                
                if(has_bottom_left==1)
                    if(has_bottom==0 && has_left==0)
                        xnew(yDown, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xLeft)=valueForChangedLocations;
                         doPlot=1;
                             type='bottom left';
                    end
                end
                
                if(has_bottom_right==1)
                    if(has_bottom~=1 && has_right~=1)
                        xnew(yDown, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xRight)=valueForChangedLocations;
                         doPlot=1;
                             type='bottom right';
                    end
                end
                
                
                
            else
                % If an isolated hole, then fit it in.
                if(has_top==1 && has_right==1 && has_bottom==1 && has_left==1)
                    xnew(yCurrent,xCurrent)=valueForChangedLocations;
                     doPlot=1;
                        type='fill hole';
                end
                
            end
            
            if(masterDoPlot ==1 && doPlot==1)
                range = 15;
                xCoordinates = max(1,elx-range):min(nelx,elx+range);
                yCoordinates = max(1,ely-range):min(nely,ely+range);
                xSubSection = xnew(yCoordinates,xCoordinates);
                xOldSubSection = xold(yCoordinates,xCoordinates);
                
                subplot(1,2,1);
                p.PlotArrayGenericWithBlueWhiteColors(xOldSubSection,sprintf('before, element %i %i type: %s',elx,ely,type)); % plot the results.
                subplot(1,2,2);
                p.PlotArrayGenericWithBlueWhiteColors(xSubSection,'after'); % plot the results.
                tttt=1;
                
            end
        end
    end
    
    afterSum=sum(sum(xnew));
    densityChangedLocal=afterSum- localSumBefore;
    fprintf('Fixing Corners: Loop %i Density change %f\n',loopcount,densityChangedLocal);
    loopcount=loopcount+1;
    
    
end
doPlot =0;
if(doPlot ==1)
   
    figure(2)
    
    
    [idum,hostname]= system('hostname');
    hostname=strtrim(hostname);
    mycomputerName = 'LAPTOP-KQHSCJB1';
    
    if(strcmp(hostname,mycomputerName)~=1) % if NOT running on my laptop
        p.PlotArrayGeneric(xnew-x,'diff'); % plot the results.
        print('./CheckForCornersDiffPlot.png','-dpng', '-r1200')
    else
        subplot(1,2,1);
        p.PlotArrayGenericWithBlueWhiteColors(x,'before'); % plot the results.
        subplot(1,2,2);
        p.PlotArrayGenericWithBlueWhiteColors(xnew,'after'); % plot the results.
        %         subplot(1,3,3);
        %         p.PlotArrayGeneric(xnew-x,'diff'); % plot the results.
        %         print('./CheckForCornersDiffPlot.png','-dpng');
    end
    
    %       close all
end
fprintf('Finished Fixing Corners\n');
numChanged=afterSum-beforeSumGlobal;