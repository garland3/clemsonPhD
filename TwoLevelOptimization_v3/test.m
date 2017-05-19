function [outvalue] =test()
a=2;
b=2;



config = Configuration;
config.mode =65;

outvalue=a+b+config.mode;
end