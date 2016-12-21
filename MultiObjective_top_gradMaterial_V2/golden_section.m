function [ min, pointsArray ] = golden_section(funct, sym_var, lower_bracket,higher_bracket,epsilon )
% display('golden section');
% format('shortG');

debug = false;
verbosity = 0;
if(debug ==true)
    display('golden_section search');
end

gr =  (1+sqrt(5))/2 -1 ; % golden ratio  0.618033988749895

leng = higher_bracket - lower_bracket; % leng is length
grleng = leng*gr ; % golden ratio lenth

x0 = lower_bracket;
x1 = higher_bracket - grleng;
x2 = lower_bracket + grleng;
x3 = higher_bracket;




%fx0 = subs(funct, sym_var, x0);
fx1 = subs(funct, sym_var, x1);
fx2 = subs(funct, sym_var, x2);
%fx3 = subs(funct, sym_var, x3);

new_row =[0 x0 x1 x2 x3 fx1 fx2] ; % new row with values x & y
pointsArray =  new_row;

n = 1;
while(1 == 1)
    if(debug == true && verbosity ==1)
        str = sprintf('loop# = %d, x0 = %f, x1 = %f, x2 = %f, x3 = %f, fx1 = %f, fx2 = %f\n', n, x0, x1, x2, x3, fx1, fx2); display(str);     
    end
     
     
    if(fx1<=fx2) % less than or equal
        
        % x0 = x0; % x0 stays the same
        x3 = x2; % the old x2 is now x3        
        x2 = x1; % the old x1 is now x2
        fx2 = fx1;
       
        
        leng = x3 - x0; % find the length of the interval
        x1 = x3 - leng*gr; % find golden ratio of length, subtract it from the x3 value
        fx1 = subs(funct, sym_var, x1); % calculate the fx
        
    elseif(fx1>fx2) % greater than
        x0 = x1; % the old x1 is now x0
        x1 = x2; % the old x2 is now the new x1
        fx1 = fx2;
        % x3 = x3; % x3 stays the same.
        
        leng = (x3 - x0); % find the length of the interval
        x2 = x0 + leng*gr; % find golden ratio of length, subtract it from the x3 value
        fx2 = subs(funct, sym_var, x2); % calculate the fx
    end
    
    % check to see if we are as close as we want
   
    
    if(leng < epsilon) 
        break;
    end
   
    
   
    
    % Save the results. 
    new_row =[n x0 x1 x2 x3 fx1 fx2] ; % new row with values x & y
    pointsArray = [pointsArray ; new_row];
     n = n +1; % increment
     
     if(n>100)
         
         break;
     end
    
end

if(debug == true)
   display(pointsArray);
end

min = (x2 + x3)/2;
end

