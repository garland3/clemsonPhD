close all; clc; clear;
input = zeros(11,16);

input(1,:) = csvread('results0.csv');
input(2,:) = csvread('results1.csv');
input(3,:) = csvread('results2.csv');
input(4,:) = csvread('results3.csv');
input(5,:) = csvread('results4.csv');
input(6,:) = csvread('results5.csv');
input(7,:) = csvread('results6.csv');
input(8,:) = csvread('results7.csv');
input(9,:) = csvread('results8.csv');
input(10,:) = csvread('results9.csv');
input(11,:) = csvread('results10.csv');

deltaRSize = 1/40;
t = linspace(0,1,1/deltaRSize);

% % ----------------------
% % Plot the first 1-6 graphs
% %----------------------------

% % ----------------------
% % Sub Plot 1
% %----------------------------

firsgraphNums = 1:6;
subplot(2,2,1)
for i = firsgraphNums    
    [hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(input(i,:));
     h =  bezierInter(hpoints,t);
     generateCurveForOpenSCAD(h,i)
   plot(h(1,:),h(2,:));
   hold all
end

legendNum = 0:1/10:1
legend2 = legendNum(firsgraphNums)
C=num2str(transpose(legend2));
H=legend(C);
set(H,'Interpreter','none');
xlabel('radial distance, meters') % y-axis label
ylabel('Height') % x-axis label
tti= strcat('Profile  of flywheel for each weight');
title(tti);

% % ----------------------
% % Sub Plot 2
% %----------------------------

subplot(2,2,2)
for i = firsgraphNums    
    [hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(input(i,:));
     h =  bezierInter(hpoints,t);   
    volFrac =  bezierInter(VolFractpoints,t);
     plot(volFrac(1,:),volFrac(2,:));
       hold all
end

xlabel('radial distance, meters') % y-axis label
ylabel('Vol Fract SN, heavy') % x-axis label
tti= strcat('Vol Fract SN, heavy');
title(tti);

legend2 = legendNum(firsgraphNums)
C=num2str(transpose(legend2));
H=legend(C);
set(H,'Interpreter','none');

% % ----------------------
% % Plot the first 7-11 graphs
% %----------------------------

% % ----------------------
% % Sub Plot 3
% %----------------------------

secondgraphNums = 7:11
subplot(2,2,3)
for i = secondgraphNums    
    [hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(input(i,:));
     h =  bezierInter(hpoints,t);
   plot(h(1,:),h(2,:));
   hold all
end

legend2 = legendNum(secondgraphNums)
C=num2str(transpose(legend2));
H=legend(C);
set(H,'Interpreter','none');

xlabel('radial distance, meters') % y-axis label
ylabel('Height') % x-axis label
tti= strcat('Profile  of flywheel for each weight');
title(tti);

% % ----------------------
% % Sub Plot 4
% %----------------------------

subplot(2,2,4)
for i = secondgraphNums    
    [hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(input(i,:));
     h =  bezierInter(hpoints,t);   
    volFrac =  bezierInter(VolFractpoints,t);
     plot(volFrac(1,:),volFrac(2,:));
       hold all
end

xlabel('radial distance, meters') % y-axis label
ylabel('Vol Fract SN, heavy') % x-axis label
tti= strcat('Vol Fract SN, heavy');
title(tti);

legend2 = legendNum(secondgraphNums)
C=num2str(transpose(legend2));
H=legend(C);
set(H,'Interpreter','none');
