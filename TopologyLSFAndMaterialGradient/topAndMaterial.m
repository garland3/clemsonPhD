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
recordResultsInFile = 1;

% ------------------------
% Algorithm configurations
% ------------------------

mode = 3; % 1 = optimize only material, 2 optimize only topology, 3 = both, mode 4 = heat transfer only
nelx = 20; % 40 % number of elements in the x direction
nely = 9; % 18 % number of elements in the y direction
v1 = 0.5; % amount of material 1 to allow where  1 = 100%
v2 = 0.5; % amount of mateiral 2 to allow
lambda1  = 0; % Set the lambda1 penalty/ lagrangian
dampingLambda1 = 0.5; % 0.1
g1Multiplier = 20;
timestep = 0.1; % step size for vol fraction update, influenced by lambda1

dampingtopology = 5; %0.5
g2Multiplier = 1;
lambda2 = -0.1; % set the lambda2 penality/lagrangian
mu1 = 5; % set the penality term close to zero
mu2 = 5; % set the second penality close to zero.
omegaMin = 0.1; % set the minimum allowed vol fraction of material 1 (stronger)
omegaMax = 0.9; % set the max allowed vol fraction of material 2 (weaker)
alpha = 0; % set the term that influenes the smoothness of the vol fraction
beta = 0; % set the perimeter regularization term.


stepsVolfraction = 15;
stepsTopology = 15;

maxFEAcalls = 30;


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
rM = 1;
structure = ones(nely*rM, nelx*rM); % initialize the structure as the whole domain
volFraction = structure*0.5; % volFraction(1:3:nely, 1:3:nelx) = 1; % initialize the volfraction composition

if mode ==1 % optimize material distribution only, 50%
    v1 = 0.5;
    v2 = 0.5; 
elseif mode ==2 % optimize the topology only, 50% material 1
    v1 = 0;
    v2 = 0.39;
    volFraction = structure*0.1/(0.1+0.29);
   
elseif mode ==3
    v1 = 0.10;
    v2 = 0.29;
    
%     if( 1== 0)
%         holesize = 5;
%         tempCount = 0;
%         for j = 1:5:nely-holesize
%            for i = 1:5:nelx-holesize
%                for jj = 1:holesize
%                    for ii = 1:holesize
%                        if(mod(tempCount ,2) ==0)
%                          structure(j+jj,i+ii)  = 0;
%                        end
%                    end
%                end
%                tempCount = tempCount+1;
%            end
%         end
%     end
    
elseif mode ==4 % optimize the topology only, using heat transfer
    v1 = 0;
    v2 = 0.39;
    volFraction = structure*0.1/(0.1+0.29);
    lambda2 = 0;
    mu2 = 0.5;
    
    dampingtopology = 10;
    g2Multiplier = 0.001;
    
    
    for j = 1:nely
       for i = 1:nelx
           if(randi(100)<v2*100)
             
                         structure(j,i)  = 1;
           else
                       structure(j,i)  = 0;
           end
               
       end
    end
  
end

totalVol = v1+v2; % the total volume of the whole structure
percentV1 = v1/totalVol*100;


[lsf] = reinit(structure); % make the lsf whic is the signed distance for the structure edge

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

FEAcount = 1;
count = 0;
while(FEAcount<maxFEAcalls)
    count = count +1;
    % --------------------
    % Run the FEA
    % --------------------
   % [U, g1_local_square,g2_local_square, volFracV1, volFracV2] =  FEALevelSet_2D(structure,lsf,volFraction,  doPlot, alpha,beta, count, mode, rM); %#ok<ASGLU>
    
    % --------------------
    % Update the volume fraction composition
    % --------------------
    if (mode ==1 || mode ==3)
        if(mod(count, 3) ==0 || mode ==1)
            for subcount1 = 1:stepsVolfraction
                FEAcount= FEAcount+1;
                [U, g1_local_square,g2_local_square, g3_local_square, volFracV1, volFracV2] =  FEALevelSet_2D(structure,lsf,volFraction,  doPlot, alpha,beta, FEAcount, mode, rM); %#ok<ASGLU>
                   
               totalStainE = sum(g2_local_square(:));
                totalVolLocal = volFracV1+ volFracV2;
               
                percentV1Local = volFracV1/totalVolLocal*100;
%                 G1 = g1_local_square*g1Multiplier - lambda1 +1/(mu1)*(v1-volFracV1)^2;
                  G1 = g1_local_square*g1Multiplier - lambda1 +1/(mu1)*(percentV1-percentV1Local);
                   volFraction_proposedUpdate = volFraction+timestep*G1;



                   volFraction = max(min(volFraction_proposedUpdate,omegaMax),omegaMin);           
%                    lambda1 =  lambda1 -1/(mu1)*(v1-volFracV1)*dampingLambda1;
                    lambda1 =  lambda1 -1/(mu1)*(percentV1-percentV1Local)*dampingLambda1;

                    for i = 1:nelx
                        for j = 1:nely   
                             structureLocal = structure(j,i);
                            if(structureLocal == 0) % if void region
                                % volFraction_proposedUpdate(j,i) = 0;  
                                 volFraction(j,i) =0;
                                % volFraction(j,i) =infill;
                            end
                       end
                    end
                    
                    if recvid==1
                        F(vid) = getframe(figure(1)); %#ok<AGROW> %Get frame of the topology in each iteration
                        writeVideo(vidObj,F(vid)); %Save the topology in the video
                        vid=vid+1;
                    end
                    
                     fprintf('Volfrac1, Volfra2, feacount, lambda1, lambda2, strainE,  %0.2f , %0.2f, %d, %0.2f , %0.2f, %0.2f \n', volFracV1,volFracV2, FEAcount, lambda1, lambda2, totalStainE);
                    
                    
            end
        end
    end
   
    
    % --------------------
    % Update the boundary (level set)
    % Calculate the propogation of the levelset boundary
    % --------------------
    if (mode ==2 || mode ==3 || mode == 4)   
     
         for subcount2 = 1:stepsTopology
             FEAcount= FEAcount+1;
                   [U, g1_local_square,g2_local_square, g3_local_square, volFracV1, volFracV2] =  FEALevelSet_2D(structure,lsf,volFraction,  doPlot, alpha,beta, FEAcount, mode, rM); %#ok<ASGLU>
                  
                    if (mode ==2 || mode ==3 )
                        totalStainE = sum(g2_local_square(:));
                          fprintf('Volfrac1, Volfra2, feacount, lambda1, lambda2, strainE,  %0.2f , %0.2f, %d, %0.2f , %0.2f, %0.2f \n', volFracV1,volFracV2, FEAcount, lambda1, lambda2, totalStainE);

                         totalVolLocal = volFracV1+volFracV2;
                          G2 = g2_local_square* g2Multiplier -lambda2+1/(mu2)*(totalVol-totalVolLocal);

                        % topologySens_square =  topologySens_square +pi*(lambda2 -1/(mu2)*(totalVol-totalVolLocal)*dampingtopology);

    %                      shapeSens = shapeSens - la + 1/La*(volCurr-volReq);
    %                      topSens = topSens + pi*(la - 1/La*(volCurr-volReq));

                        lambda2 =  lambda2 -1/(mu2)*(totalVol-totalVolLocal)*dampingtopology;
                        
                         % ----------------------
                         % mode 4 is heat transfer topology only
                         %-----------------------
                    elseif (mode ==4)
                        
                        totalStainE = sum(g3_local_square(:));
                        fprintf('Volfrac1, Volfra2, feacount, lambda1, lambda2, heat strainE,  %0.2f , %0.2f, %d, %0.2f , %0.2f, %0.2f \n', volFracV1,volFracV2, FEAcount, lambda1, lambda2, totalStainE);

                        totalVolLocal = volFracV1+volFracV2;
                        G2 = g3_local_square*g2Multiplier -lambda2+1/(mu2)*(totalVol-totalVolLocal);
                        lambda2 =  lambda2 -1/(mu2)*(totalVol-totalVolLocal)*dampingtopology;  

                    end
                    
                 topSensFull = zeros(size(topologySens_square)+2); topSensFull(2:end-1,2:end-1) = topologySens_square;
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


                   % lsf = lsf - dt*moveSlope- topologySensWeight*dt*topSensFull;
                      lsf = lsf - dt*moveSlope ;%- topologySensWeight*dt*topSensFull;
                    elseif (1==0)

                        dt = 0.1/max(abs(g2Full(:)));
                        % Evolve for total time stepLength * CFL value:
                        for i = 1:(10*stepLength)
                         % Calculate derivatives on the grid
                          dpx = circshift(lsf,[0,-1])-lsf;
                         dmx = lsf - circshift(lsf,[0,1]);
                         dpy = circshift(lsf,[-1,0]) - lsf;
                         dmy = lsf - circshift(lsf,[1,0]);
                         % Update level set function using an upwind scheme 
                         lsf = lsf - dt * min(g2Full,0).* ...
                           sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2 ) ...
                           - dt * max(g2Full,0) .*...
                           sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2 );
                        end
                    elseif(1 == 0)
                         diff_forward_x = circshift(lsf,[0,-1])-lsf;
                        diff_backward_x = lsf - circshift(lsf,[0,1]);

                        diff_forward_y = circshift(lsf,[-1,0])-lsf;
                        diff_backward_y = lsf - circshift(lsf,[+1,0]);

                        pNabla = (max(diff_backward_x,0).^2+min(diff_forward_x,0).^2 + max(diff_backward_y,0).^2 +min(diff_forward_y,0).^2).^(1/2);
                        nNabla = (max(diff_forward_x,0).^2+min(diff_backward_x,0).^2 + max(diff_forward_y,0).^2 + min(diff_backward_y,0).^2).^(1/2);

                        lsf = lsf - dt*(max(g2Full, 0).*pNabla+min(g2Full,0).*nNabla);


                    end


            strucFull2 = (lsf>0.001);
            structure = strucFull2(2:end-1,2:end-1);


            if(mod(FEAcount, 6) ==0)
                [lsf] = reinit(structure)   ;
            end
            
            if recvid==1
                F(vid) = getframe(figure(1)); %#ok<AGROW> %Get frame of the topology in each iteration
                writeVideo(vidObj,F(vid)); %Save the topology in the video
                vid=vid+1;
            end
         end
    end
    
    %fprintf('Volfrac1, Volfra2, feacount, lambda1, lambda2, mu2,  %0.2f , %0.2f, %d, %0.2f , %0.2f, %0.2f \n', volFracV1,volFracV2, FEAcount, lambda1, lambda2, mu2);
    %[volFracV1,volFracV2, count, lambda1, lambda2]
   
    
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
        lsf = bwdist(~strucFull)*0.5;
%         lsfn = -bwdist(strucFull);
        
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

