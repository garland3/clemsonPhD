function xSensitivity = checkForSupports(xCutOff)

[m,n]=size(xCutOff);

% 0 sensitivy means material is likely to be added
% higher values mean it is less likely

xSensitivity = xCutOff;

for yy = 1:n
    for xx = 1:m
       if (yy == 1) % then bottom row
             xSensitivity(xx,yy) = 0;
       else % if not the bottom row
           
           % if the 3 to the left, midddle, and right on the row below it
           % are empty then discourage adding material. AND, find out how
           % much it costs to add a support
           
           % if on the left wall
           if(xx ==1)
               if xCutOff(xx,yy-1) == 1 || xCutOff(xx+1, yy-1) ==1
                    xSensitivity(xx,yy) = 0;
               else
                     xSensitivity(xx,yy) = 1;
               end
               
           % if on the right wall
           elseif (xx == m)
                if xCutOff(xx-1,yy-1) == 1 || xCutOff(xx, yy-1) ==1
                    xSensitivity(xx,yy) = 0;
               else
                     xSensitivity(xx,yy) = 1;
               end
               
               
           % if a normal case
           else
                if xCutOff(xx-1,yy-1) == 1 || xCutOff(xx, yy-1) ==1  || xCutOff(xx+1, yy-1) ==1
                    xSensitivity(xx,yy) = 0;
               else
                     xSensitivity(xx,yy) = 1;
               end
               
           end
       
       end
    
    end
end