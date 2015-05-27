function testArgs(x,y)

disp(x)
disp(y)

% qsub -I 
%   mcc -v -N -R -nodisplay  -R  -nojvm -m testArgs.m
