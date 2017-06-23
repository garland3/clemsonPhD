function [xnew] = CheckRemoveOrphanedSegments(x, nelx, nely, cutoff)

% search the main array
% find the first solid element
% AAA: search up,down, left right for other solid elements, if Solid then
% 1. Add to future search array
% 2. Remove itself from future search (just finished searching this one)
% 3. Addself to SolidBody Array
% 3. Look at the future search array and go to step AAA

% if total area of solidbody array<10%, then find a new starting element
% that is NOT in solidbody
% else, remove all elements that are not part of solidbody

% Trying to avoid recursive search, might result in memory overflow
totalArea = nelx*nely;

xold=x;
% begin search for 1st solid element
 futureSerachArray = x*0;
AlreadySearchedArray=x*0;
AllSearchesAlreadyComplete = x*0;
for elx = 1:nelx
    for ely = 1:nely
        InPreviousBody=AllSearchesAlreadyComplete(ely,elx);
        
        xCurrentSimp=xold(ely,elx);
        if(xCurrentSimp>cutoff && InPreviousBody==0)
            futureSerachArray = x*0;
            AlreadySearchedArray=x*0;
            
            
            
            % add myself to future search
            % keep searching over the future search
            futureSerachArray(ely, elx)=1;
            
            MoreElementsToSearch=   sum(sum(futureSerachArray));
            while(MoreElementsToSearch>0)
                % search the futureSerachArray
                for elx2 = 1:nelx
                    for ely2 = 1:nely
                        if(futureSerachArray(ely2,elx2)==1)
                            
                            xCurrent=elx2;
                            xRight=xCurrent+1;
                            xLeft=xCurrent-1;
                            
                            yCurrent=ely2;
                            yUp=yCurrent+1;
                            yDown=yCurrent-1;
                            
                            % Add to the AlreadySearchedArray so we dont'
                            % check this one again. 
                            AlreadySearchedArray(ely2,elx2)=1;
                            futureSerachArray(ely2,elx2)=0;
                            
                            
                            
                            % -----------------------------
                            % Check up
                            % -------------------
                            if(yUp<=nely)
                                density = xold(yUp, xCurrent);
                                checkIfSeached= AlreadySearchedArray(yUp, xCurrent);
                                if(density>cutoff && checkIfSeached==0)
                                    futureSerachArray(yUp, xCurrent)=1;
                                    
%                                     has_top=1;
                                end
                            end
                            
                            
                            % -----------------------------
                            % Check Right
                            % -------------------
                            if(xRight<=nelx)
                                
                                density = xold(yCurrent, xRight);
                                checkIfSeached= AlreadySearchedArray(yCurrent, xRight);
                                if(density>cutoff && checkIfSeached==0)
                                    has_right=1;
                                    futureSerachArray(yCurrent, xRight)=1;
                                end
                            end
                            
                            
                            % -----------------------------
                            % Check Left
                            % -------------------
                            if(xLeft>0)
                                checkIfSeached= AlreadySearchedArray(yCurrent, xLeft);
                                density = xold(yCurrent, xLeft);
                                if(density>cutoff && checkIfSeached==0)
                                    has_left=1;
                                       futureSerachArray(yCurrent, xLeft)=1;
                                end
                            end
                            
                            
                            % -----------------------------
                            % Check down
                            % -------------------
                            if(yDown>0)
                                 checkIfSeached= AlreadySearchedArray(yDown, xCurrent);
                                density = xold(yDown, xCurrent);
                                if(density>cutoff && checkIfSeached==0)
                                    has_bottom=1;
                                     futureSerachArray(yDown, xCurrent)=1;
                                end
                            end
                        end
                    end
                end
                
                 totalAreaOfBody = sum(sum(AlreadySearchedArray));
                 MoreElementsToSearch=   sum(sum(futureSerachArray));
                 fprintf('searched %i %i, Number of more elements to search %i and body size %i of %i totalArea possible\n',xCurrent, yCurrent,MoreElementsToSearch,totalAreaOfBody,totalArea);
            end
            
            p = plotResults;
            p.PlotArrayGeneric(AlreadySearchedArray,'body');

            totalAreaOfBody = sum(sum(AlreadySearchedArray));
            if(totalAreaOfBody>totalArea*0.05)
                % we found the main body,so we can exit. 
                xnew=AlreadySearchedArray;
                
                fprintf('Body has more than 5 percent of areay. Finished\n');
                return;
            else
                fprintf('Body has less  than 5  percent of areay.  NOTFinished\n');
            end
            AllSearchesAlreadyComplete=AllSearchesAlreadyComplete+AlreadySearchedArray;
            
        end
    end
end



% search up,down, left right for other solid elements, if Solid and NOT in OK array
% then
% add each element to an OK array
% When finished, if total areay is more than 10%, then assum main body
% if main body, then
% -- if NOT in the OK array, then turn to empty



%-- if more than 10% of total area, then assume this is the main body
% if not, then STart over
% start
