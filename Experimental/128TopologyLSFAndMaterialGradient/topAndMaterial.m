function [  ] = topAndMaterial(  )
%UNTITLED Trying to implement the gradient and material optimization in the
%paper Simultaneous optimization of the material properties and the topology of functionally graded structures
%   Detailed explanation goes here
% Following the step by step process that is used in the paper, starting on
% page 9. 
% Also, I plan to reuse some parts of the short level set top optimization
% presented in A discrete level-set topology optimization code.pdf

% Step 1
% Step 1: Initialize the embedding level set function ? ( x , 0 )
% at t = 0. A general treatment is to set ? ( x ) to be the
% signed distance to the given boundary of the initial design
% ?. Initialize the volume fraction ? . Set the parameters
% V1 , V2 , ?01, ?02, µ 1 , µ 2 , ? m , ? M , ?, ?, ?t?.

clear
clc
% Setup plotting stuff to help with visualization of what is happening. 
subplotX = 2;
subplotY = 1;
subplotCount = 1;
doPlot = 1


nelx = 10;
nely = 5;
time = 0;
v1 = 0.3; % amount of material 1 to allow where  1 = 100%
v2 = 0.3; % amount of mateiral 2 to allow
lambda1  = 0.01; % Set the lambda1 penalty/ lagrangian
lambda2 = 0.01; % set the lambda2 penality/lagrangian
mu1 = 0.002; % set the penality term close to zero
mu2 = 0.2; % set the second penality close to zero.     
omegaMin = 0.1; % set the minimum allowed vol fraction of material 1 (stronger)
omegaMax = 0.9; % set the max allowed vol fraction of material 2 (weaker)
alpha = 0.5; % set the term that influenes the smoothness of the vol fraction
beta = 5; % set the perimeter regularization term. 


% initialize the domain as ones
structure = ones(nely, nelx); 

% plotting stuff
if(doPlot ==1)
    figure(1)
    subplot(subplotX,subplotY,subplotCount); subplotCount=subplotCount+1;
    colormap(gray); imagesc(-structure); axis equal; axis tight; axis off;pause(1e-6);
end

% make the signed distance for the phi function
[lsf] = reinit(structure);

% plottting stuff
if(doPlot ==1)
    figure(1)
    subplot(subplotX,subplotY,subplotCount); subplotCount=subplotCount+1;
    colormap(gray); imagesc(-lsf); axis equal; axis tight; axis off;pause(1e-6);
end


% Define the signed distance transform. 
function [lsf] = reinit(struc)    
% Add a boundary around the actual domain and set it to zero. 
strucFull = zeros(size(struc)+2); strucFull(2:end-1,2:end-1) = struc;
% Use "bwdist" (Image Processing Toolbox)
% lsf = (~strucFull).*(bwdist(strucFull)-0.5) - strucFull.*(bwdist(strucFull-1)-0.5);
if(strucFull <=0)
    strucFull = 0;
end
lsf = bwdist(~strucFull);


end
end

