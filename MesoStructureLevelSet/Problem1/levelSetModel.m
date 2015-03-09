% Copyright Anthony Garland 2015 

targetDisplacment = -0.015;
tolerance = 0.003;
maxNumIterations = 29;
moveLimitLarge = 0.2;
moveLimitsmall = 0.2;




% ----------------------------------------------------------------
% Start
% --------------------------------------------------------------

% original model
xpoints = 0:0.5:10; % X position of the control points

% ypoints = 1-xpoints.^2/80 % deceasing as you go right. Starts at 0.75

%ypoints = ones(1,21); %displacment is   -0.0104
 %ypoints = zeros(1,21);%    -0.0207
% ypoints = 0.5-0.5*xpoints.^2/25;
ypoints = 0+xpoints/10

[ydisplacment, strainEnergies,maxVonMisesInRegion,maxDisplacementInRegion]= FEALevelSet(xpoints,ypoints); % Call FEA for first time

% figure(2) 
% figure(2)
%subplot(2,2,4)

%strainEnergies
historyOfTotalStrainEnergy = zeros(maxNumIterations);

for i=1:maxNumIterations
    
    % update
    if(ydisplacment<(targetDisplacment-tolerance))
        % Deflecting too far
        % move all control points up to make more stiff
        ypoints = ypoints+moveLimitLarge;
    elseif(ydisplacment>(targetDisplacment+tolerance))
        % Deflection is not enough
       
          ypoints = ypoints-moveLimitLarge;
        
    else
        % Deflection target met, so tweak the design
        
         % regions with high strain energy move up (20)
         % Regions with low strain energy move down
         numRegions = size(strainEnergies,1);
        [C,ia,ic] = unique(strainEnergies,'sorted');
        for j =1:numRegions
            rank = ic(j);
            
            ideaNumber = 2;
            if(ideaNumber ==1)
                % if top 20% strain energy region
                if(rank>(numRegions-0.2*numRegions))
                    ypoints(j)  = ypoints(j) + moveLimitsmall;

                % if bottom 20% strain energy region
                elseif(rank<numRegions*0.2)
                    ypoints(j)  = ypoints(j) - moveLimitsmall;
                end
            elseif(ideaNumber ==2)
                % move based on the regions ranking. 
                move = -moveLimitsmall+rank*(2*moveLimitsmall/numRegions);
                
                 ypoints(j)  = ypoints(j) + move;
                
            elseif(ideaNumber ==3)
                % move based on the strain energy magnitude. Everything
                % with above average strain energy you move up based on its
                % distance from the average. 
            end
            
        end
        
    end
    
    % record the total strain energy
    
    historyOfTotalStrainEnergy(i) = sum(strainEnergies);
   %  if(i>3)
    
    % Call FEA
    [ydisplacment, strainEnergies,maxVonMisesInRegion,maxDisplacementInRegion]= FEALevelSet(xpoints,ypoints); 
    %close all
  %  figure(2)
  %  bar(strainEnergies)
end




