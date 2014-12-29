close all; clc; clear;
input = zeros(11,16);

innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2
global w;
w = 630; %'omega'

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


for i = 1:11
      [hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(input(i,:));
      stress = CalculateMaxStress(hpoints,VolFractpoints,0) ;
      energy = CalculateEnergy(hpoints, VolFractpoints,w,innerR, outerR);
      mass =  CalculateMass(hpoints, VolFractpoints,innerR, outerR);
 
      
    output(i,:) = [  stress energy mass]
end


csvwrite('Stress_energy_mass2.csv',output)
output