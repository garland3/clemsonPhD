function [xnew, numChanged] = CheckForConerElements(x, nelx, nely, cutoff)
% config = mesoConfig;


% Get the density field
x(x>cutoff)=1;
x(x<cutoff)=0.001;
xnew = x;
% xold=x;

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
            if(xCurrentSimp>cutoff)
                
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
                
                % -----------------------------
                % Check up
                % -------------------
                if(yUp<=nely)
                    density = xold(yUp, xCurrent);
                    if(density>cutoff)
                        has_top=1;
                    end
                end
                
                
                % -----------------------------
                % Check Right
                % -------------------
                if(xRight<=nelx)
                    
                    density = xold(yCurrent, xRight);
                    if(density>cutoff)
                        has_right=1;
                    end
                end
                
                % -----------------------------
                % Check Right TOP
                % -------------------
                if(xRight<=nelx &&  yUp<=nely)
                    
                    density = xold(yUp, xRight);
                    if(density>cutoff)
                        has_top_right=1;
                    end
                end
                
                
                % -----------------------------
                % Check  Left TOP
                % -------------------
                if(xLeft>0 &&  yUp<=nely)
                    
                    density = xold(yUp, xLeft);
                    if(density>cutoff)
                        has_top_left=1;
                    end
                end
                
                
                % -----------------------------
                % Check Left
                % -------------------
                if(xLeft>0)
                    
                    density = xold(yCurrent, xLeft);
                    if(density>cutoff)
                        has_left=1;
                    end
                end
                
                % -----------------------------
                % Check Left Bottom
                % -------------------
                if(xLeft>0 &&  yDown>0)
                    
                    density = xold(yDown, xLeft);
                    if(density>cutoff)
                        has_bottom_left=1;
                    end
                end
                
                
                % -----------------------------
                % Check  Right Bottom
                % -------------------
                if(xRight<=nelx &&   yDown>0)
                    
                    density = xold(yDown, xRight);
                    if(density>cutoff)
                        has_bottom_right=1;
                    end
                end
                
                
                % -----------------------------
                % Check down
                % -------------------
                if(yDown>0)
                    
                    density = xold(yDown, xCurrent);
                    if(density>cutoff)
                        has_bottom=1;
                    end
                end
                
                
              
                
                % ---------------------
                %
                %         ADD MATERAIL
                %
                %
                % if has a corner, but no top or side, then add material.
                valueForChangedLocations=1;
                if(has_top_right==1)
                    if(has_top~=1 && has_right~=1)
                        xnew(yUp, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xRight)=valueForChangedLocations;
                    end
                end
                
                if(has_top_left==1)
                    if(has_top~=1 && has_left~=1)
                        xnew(yUp, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xLeft)=valueForChangedLocations;
                    end
                end
                
                if(has_bottom_left==1)
                    if(has_bottom~=1 && has_left~=1)
                        xnew(yDown, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xLeft)=valueForChangedLocations;
                    end
                end
                
                if(has_bottom_right==1)
                    if(has_bottom~=1 && has_right~=1)
                        xnew(yDown, xCurrent)=valueForChangedLocations;
                        xnew(yCurrent, xRight)=valueForChangedLocations;
                    end
                end
                
                
            end
        end
    end
    
    afterSum=sum(sum(xnew));
    densityChangedLocal=afterSum- localSumBefore;
    fprintf('Fixing Corners: Loop %i Density change %f\n',loopcount,densityChangedLocal);
    loopcount=loopcount+1;
    
    
end
doPlot =1;  
if(doPlot ==1)
    p = plotResults;
    figure(2)
    
    
    [idum,hostname]= system('hostname');
    hostname=strtrim(hostname);
    mycomputerName = 'LAPTOP-KQHSCJB1';
    
    if(strcmp(hostname,mycomputerName)~=1) % if NOT running on my laptop
        p.PlotArrayGeneric(xnew-x,'diff'); % plot the results.
        print('./CheckForCornersDiffPlot.png','-dpng', '-r1200')
    else
        subplot(1,3,1);
        p.PlotArrayGeneric(x,'before'); % plot the results.
        subplot(1,3,2);
        p.PlotArrayGeneric(xnew,'after'); % plot the results.
        subplot(1,3,3);
        p.PlotArrayGeneric(xnew-x,'diff'); % plot the results.
        print('./CheckForCornersDiffPlot.png','-dpng');
    end
    
    %       close all
end
fprintf('Finished Fixing Corners\n');
numChanged=afterSum-beforeSumGlobal;