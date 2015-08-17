function [  ] = topAndMaterial(  )
%topAndMaterial Trying to implement the gradient and material optimization in the
%paper Simultaneous optimization of the material properties and the topology of functionally graded structures
%   Detailed explanation goes here
% Following the step by step process that is used in the paper, starting on
% page 9.
% Also, I plan to reuse some parts of the short level set top optimization
% presented in A discrete level-set topology optimization code.pdf
% Copyright Anthony Garland 2015

clear
clc
close all

% ------------------------------
% Declare all configuration varriables
% -------------------------------

% ------------------------
% Setup plotting stuff to help with visualization of what is happening.
% ------------------------
subplotX = 2; % number of subplots in X direction
subplotY = 3; % number of subplots in the Y direction
subplotCount = 1; % current subplot count (do not modify)
doPlot = 0; % Controls plotting or not
plotStructure = doPlot;
recvid = 0; % Record a video of the figure 1, record view, 1 = yes

% ------------------------
% Algorithm configurations
% ------------------------

mode = 1; % 1 = optimize only material, 2 optimize only topology, 3 = both
nelx = 40; % number of elements in the x direction
nely = 10; % number of elements in the y direction
v1 = 0.5; % amount of material 1 to allow where  1 = 100%
v2 = 0.5; % amount of mateiral 2 to allow
lambda1  = 0; % Set the lambda1 penalty/ lagrangian
dampingLambda1 = 0.5;
timestep = 0.1; % step size for vol fraction update, influenced by lambda1

lambda2 = -0.1; % set the lambda2 penality/lagrangian
mu1 = 1; % set the penality term close to zero
mu2 = 20; % set the second penality close to zero.
omegaMin = 0.1; % set the minimum allowed vol fraction of material 1 (stronger)
omegaMax = 0.9; % set the max allowed vol fraction of material 2 (weaker)
alpha = 5; % set the term that influenes the smoothness of the vol fraction
beta = 2; % set the perimeter regularization term.


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
structure = ones(nely, nelx); % initialize the structure as the whole domain
volFraction = structure*0.5; % volFraction(1:3:nely, 1:3:nelx) = 1; % initialize the volfraction composition
[lsf] = reinit(structure); % make the lsf whic is the signed distance for the structure edge


if mode ==1 % optimize material distribution only, 50%
    v1 = 0.5;
    v2 = 0.5; 
elseif mode ==2 % optimize the topology only, 50% material 1
    v1 = 0;
    v2 = 0.5;
    volFraction = structure*0;
elseif mode ==3
    v1 = 0.3;
    v2 = 0.3;
end

% plotting the structure
if(plotStructure ==1)
    figure(1)
    subplot(subplotX,subplotY,subplotCount); subplotCount=subplotCount+1;
    imagesc(-structure); axis equal; axis tight; axis off;pause(1e-6);
end

% plottting the lsf
if(doPlot ==1)
    figure(1)
    subplot(subplotX,subplotY,subplotCount); subplotCount=subplotCount+1; %#ok<NASGU>
    colormap(gray); imagesc(-lsf); axis equal; axis tight; axis off;pause(1e-6);
end

% ------------------------------
% Step 2
% ------------------------------
%  Solve the state equation (9) via the finite element
% method to find the displacement u. Calculate the sensitivity G1.
% Update omega according to Eq. (26).

for count = 1:500
    % --------------------
    % Run the FEA
    % --------------------
    [U, g1_local_square,g2_local_square, volFracV1, volFracV2] =  FEALevelSet_2D(structure,lsf,volFraction,  doPlot, alpha,beta, count, mode); %#ok<ASGLU>
    
    % --------------------
    % Update the volume fraction composition
    % --------------------
    if (mode ==1 || mode ==3)
        if(mod(count, 5) ==0)
           G1 = g1_local_square - lambda1 +1/(mu1)*(v1-volFracV1)^2;
           volFraction_proposedUpdate = volFraction+timestep*G1;
           volFraction = max(min(volFraction_proposedUpdate,omegaMax),omegaMin);
           
           lambda1 =  lambda1 -1/(mu1)*(v1-volFracV1)*dampingLambda1;
        end
    end
   
    
    % --------------------
    % Update the boundary (level set)
    % Calculate the propogation of the levelset boundary
    % --------------------
    if (mode ==2 || mode ==3)   
        if (count < 1)
            G2 = g2_local_square+1/(mu2)*(v2-volFracV2);
         
            
        else
                 G2 = g2_local_square-lambda2+1/(mu2)*(v2-volFracV2);
                damping = 5;
                lambda2 =  lambda2 -1/(mu2)*(v2-volFracV2)*damping;
        end
%         
      
%         G2 = -g2_local_square - lambda2 + 1/mu2*(v2-volFracV2)^2;
%         lambda2 = lambda2 - 10/mu2 * (v2 - volFracV2)^2;
%         alpha = 0.95;
%         mu2 = mu2*alpha;
         
         
            g2Full = zeros(size(G2)+2); g2Full(2:end-1,2:end-1) = G2;
        
        % Choose time step for evolution based on CFL value
        dt = 1/max(abs(G2(:)));
        stepLength = 1;
        % Evolve for total time stepLength * CFL value:
        %for i = 1:(15*stepLength)
                % Calculate derivatives on the grid
                if 1==1
                diff_forward_x = circshift(lsf,[0,1])-lsf;
                diff_backward_x = lsf - circshift(lsf,[0,-1]);

                diff_forward_y = circshift(lsf,[1,0])-lsf;
                diff_backward_y = lsf - circshift(lsf,[-1,0]);

                positiveGradientUpwind = (max(diff_backward_x,0).^2 +...
                    min(diff_forward_x,0).^2 + ...
                    max(diff_backward_y,0).^2 +...
                    min(diff_forward_y,0).^2 ).^(1/2);

                negGradientUpwind = (min(diff_backward_x,0).^2 +...
                    max(diff_forward_x,0).^2 + ...
                    min(diff_backward_y,0).^2 +...
                    max(diff_forward_y,0).^2 ).^(1/2);

                moveSlope = (max(-g2Full,0).*positiveGradientUpwind + min(-g2Full,0).*negGradientUpwind);


                lsf = lsf - dt*moveSlope;
                else
                    
                    dt = 0.1/max(abs(g2Full(:)));
                    % Evolve for total time stepLength * CFL value:
                    for i = 1:(10*stepLength)
                     % Calculate derivatives on the grid
                     dpx = circshift(lsf,[0,1])-lsf;
                     dmx = lsf - circshift(lsf,[0,+1]);
                     dpy = circshift(lsf,[1,0]) - lsf;
                     dmy = lsf - circshift(lsf,[-1,0]);
                     % Update level set function using an upwind scheme 
                     lsf = lsf - dt * min(g2Full,0).* ...
                       sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2 ) ...
                       - dt * max(g2Full,0) .*...
                       sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2 );
                    end
                end

        
        strucFull2 = (lsf>0);
        structure = strucFull2(2:end-1,2:end-1);
        
     
        if(mod(count, 5) ==0)
            [lsf] = reinit(structure)   ;
        end
    end
    
    fprintf('Volfrac1, Volfra2, iter, lambda1, lambda2, mu2,  %0.2f , %0.2f, %d, %0.2f , %0.2f, %0.2f \n', volFracV1,volFracV2, count, lambda1, lambda2, mu2);
    %[volFracV1,volFracV2, count, lambda1, lambda2]
    drawnow
    if recvid==1
        F(vid) = getframe(figure(1)); %#ok<AGROW> %Get frame of the topology in each iteration
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
            strucFull = -0.1;
        end
        lsf = bwdist(~strucFull);
        lsfn = -bwdist(strucFull);
        
%         [y,x]=size(lsf);
%         for i = 1: x
%             for j = 1:y
%                 if(lsf(j,i) ==0)
%                     lsf(j,i) = lsfn(j,i);
%                 end
%             end
%         end
%       lsf  
        
    end
end

