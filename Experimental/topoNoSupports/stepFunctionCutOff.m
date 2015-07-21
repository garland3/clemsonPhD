function xOut = stepFunctionCutOff(x)

[m,n]=size(x);
xOut = x;
cutOff = 0.7;

for yy = 1:n
    for xx = 1:m
        if x(xx,yy) >= cutOff
            xOut(xx,yy) = 1;
        else
            xOut(xx,yy) = 0;
        end
    end
end
            
            
