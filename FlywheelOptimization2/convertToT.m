function u= convertToT(x,min,max)
% Coverts the x value between xMin and xMax, to a u value between 0 and 1
u = (x-min)/(max-min);