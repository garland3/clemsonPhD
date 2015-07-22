function [  ] = topAndMaterial(  )
%topAndMaterial Trying to implement the gradient and material optimization in the
%paper Simultaneous optimization of the material properties and the topology of functionally graded structures
%   Detailed explanation goes here
% Following the step by step process that is used in the paper, starting on
% page 9. 
% Also, I plan to reuse some parts of the short level set top optimization
% presented in A discrete level-set topology optimization code.pdf


clear
clc
close all

recvid = 1; % record view, 1 = yes
if recvid==1
    vidObj = VideoWriter('results.avi');    %Prepare the new file for video
    vidObj.FrameRate = 5;
    vidObj.Quality = 100;
    open(vidObj);
    vid=1;
end
% ------------------------------
% Step 1
% ------------------------------

% Step 1: Initialize the embedding level set function ? ( x , 0 )
% at t = 0. A general treatment is to set ? ( x ) to be the
% signed distance to the given boundary of the initial design
% ?. Initialize the volume fraction ? . Set the parameters
% V1 , V2 , ?01, ?02, µ 1 , µ 2 , ? m , ? M , ?, ?, ?t?.


% Setup plotting stuff to help with visualization of what is happening. 
subplotX = 2;
subplotY = 1;
subplotCount = 1;
doPlot = 0;


nelx = 12;
nely = 6;
time = 0;
v1 = 0.5; % amount of material 1 to allow where  1 = 100%
v2 = 0.5; % amount of mateiral 2 to allow
lambda1  = -5; % Set the lambda1 penalty/ lagrangian
lambda2 = 0.01; % set the lambda2 penality/lagrangian
mu1 = 0.1; % set the penality term close to zero
mu2 = 0.2; % set the second penality close to zero.     
omegaMin = 0.1; % set the minimum allowed vol fraction of material 1 (stronger)
omegaMax = 0.9; % set the max allowed vol fraction of material 2 (weaker)
alpha = 0.2; % set the term that influenes the smoothness of the vol fraction
beta = 5; % set the perimeter regularization term. 
timestep = 0.04; 


% initialize the domain as ones
structure = ones(nely, nelx); 
volFraction = structure*0; % initialize the volfraction composition
 volFraction(1:3:nely, 1:3:nelx) = 1;

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


% ------------------------------
% Step 2
% ------------------------------
%  Solve the state equation (9) via the finite element
% method to find the displacement u. Calculate the sensitivity G1.
% Update omega according to Eq. (26).

for count = 1:100

    [U, g1_local_square, volFracV1, volFracV2]  =  FEALevelSet_2D(structure,volFraction,  doPlot, alpha);
    G1 = g1_local_square - lambda1 +1/(mu1)*(v1-volFracV1);
    volFraction_proposedUpdate = volFraction+timestep*G1;

    volFraction = max(min(volFraction_proposedUpdate,omegaMax),omegaMin);

    lambda1 =  lambda1 -1/(mu1)*(v1-volFracV1);
    [volFracV1, count, lambda1]
   
    drawnow
    
    if recvid==1
        F(vid) = getframe(figure(1)); %Get frame of the topology in each iteration
        writeVideo(vidObj,F(vid)); %Save the topology in the video
        vid=vid+1;
    end
    
    

end


if recvid==1
    close(vidObj);  %close video
end










% ------------------------------
% Supporting functions. 
% ------------------------------

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

