
clear
clc
close all


t = 5; % mm, thickness of the beam
FappliedLoad = 200; % N.
Enylon = 1700; % elastic mod of nylon, Mpa
KheatNylon = 2; % W/degreeC
Epla =  3368; % elastic mod of pla, MPa
KheatPLA = 1; % W/degreeC
E_empty = 1; % elastic mod for the region with no material to prevent sigularity of the FEA
K_empty = 0.001;
v = 0.3; % Assume 0.3 for Poisson's ratio for both materials

recvid = 1;

if recvid==1
    vidObj = VideoWriter('resultsPalmetto.avi');    %Prepare the new file for video
    vidObj.FrameRate = 5;
    vidObj.Quality = 100;
    open(vidObj);
    vid=1;
end



for i = 2:256
    name = sprintf('./StuctGradOutput/gradAndStuct%i.csv',i);
   
    structGradArray = csvread(name);
    figure(1)
    imagesc(structGradArray); axis equal; axis tight; axis off;
    scale = 1;
   % surf(structGradArray/scale); axis equal; axis tight; axis off;
        set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
        % caxis([-1 1 ]);
         titletext = sprintf('Structure and Elastic Mod Gradient Step %i',i)
        title(titletext)
        colormap winter
        %  cmap = colormap;
        rgbSteps = Epla- Enylon +1 ; % plus 1 for 1 below the Enylon for void
        % [cmin,cmax] = caxis;
        caxis([Enylon-100,Epla]./scale)
        map = colormap; % current colormap
        %map = [colormap(1,1):1/rgbSteps:colormap(1:-1)
        for zz =    1:rgbSteps
            map(1,:) = [1,1,1];
            map(zz,:) = [0,       zz*7/(8*rgbSteps)+1/8,          0.5];
        end    
        colormap(map)   
        colorbar
        
        
         if recvid==1
             drawnow
            F(vid) = getframe(figure(1)); %#ok<AGROW> %Get frame of the topology in each iteration
            writeVideo(vidObj,F(vid)); %Save the topology in the video
            vid=vid+1;
        end
        
end


if recvid==1
    close(vidObj);  %close video
end