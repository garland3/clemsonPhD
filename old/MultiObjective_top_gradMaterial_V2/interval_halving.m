function [ min, pointsArray ] = interval_halving(funct, sym_var, lower_bracket,higher_bracket,epsilon )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% 3 Point equal interval search
% Consider interval [a, b], f(a), f(b); b >a
display('interval_halving');
a = lower_bracket;
b = higher_bracket;
xm = (b - a)/2 + a; % xm is x middle
pointsArray = zeros(1, 8);
found_min = false;
n = 0;



while(found_min == false)
    leng = b - a; % leng is length
    
    x_1qtr = a+leng/4; % Find the x value that is a quarter of the way between the low and high
    x_3qtr = a+3*leng/4; % Find the x value that is 3 quarters of the way between the low and high
    
    fxm = subs(funct, sym_var, xm); % evaluate the function at the middle
    fx_1qtr = subs(funct, sym_var, x_1qtr); % evaluate the function at 1 quarter
    fx_3qtr = subs(funct, sym_var, x_3qtr); % evaluate the function at 3 quarters
    
    new_row =[a x_1qtr xm x_3qtr b fx_1qtr fxm  fx_3qtr] ; % new row with values x & y
    pointsArray = [pointsArray ; new_row];

    %str = sprintf('loop# = %d, a = %f, x_1qtr = %f, xm = %f, x_3qtr = %f, b = %f, fx_1qtr = %f , fxm = %f,  fx_3qtr = %f,\n', n, a, x_1qtr, xm, x_3qtr, b, fx_1qtr, fxm,  fx_3qtr); display(str);     
    
    if(fx_1qtr < fxm)
        b = xm;
        xm = x_1qtr;
    elseif(fx_3qtr < fxm)
        a = xm;
        xm = x_3qtr;
    elseif((fx_1qtr > fxm) && (fxm < fx_3qtr))
        a = x_1qtr;
        b = x_3qtr;
        % xm = xm;
    else
        error('Flat function, no minimum exists');
    end
    
    if(leng < epsilon)
        break;
    end
    n = n +1;
end

min = xm;


end

