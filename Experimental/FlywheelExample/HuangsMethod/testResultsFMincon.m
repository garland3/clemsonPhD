% just plot the results 

% optimization setup 5.12, max energy, st. max stress <= 45
inputvars = [   0.7569    0.1432    0.4097    0.4089    0.4268    0.4668    0.3171    0.4557    0.4593    0.4733]
 

 
  %[objective] = callCaculateEnergy(x)

  numpoints = 5;

hpoints = inputvars(1,1:numpoints);
VolFractpoints =  inputvars(1,(numpoints+1):2*numpoints);

innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2

hpointsR = [innerR,0.04,0.1,0.12,outerR]; 
h = ones(2,5);
h(1,:) = hpointsR(1,:);
h(2,:) = hpoints(1,:);
h = transpose(h);

vf = ones(2,5);
vf(1,:) = hpointsR(1,:);
vf(2,:) = VolFractpoints(1,:);
vf = transpose(vf);


% ----------------------
% Just plot the curves. 
%----------------------------


doplot = 1;
hpoints = h;
 VolFractpoints = vf;
% ----------------------
% Model the height with the bezier curve
%----------------------------
deltaRSize = 1/20;



% hpoints = [innerR,0.5;
%             0.04,0.5;
%              0.08,0.5;
%               0.12,0.5;
%              outerR,0.5]
t = linspace(0,1,1/deltaRSize);
h =  bezierInter(hpoints,t);



% ----------------------
% Plot the height and control poiints
%----------------------------
if(doplot ==1)
   % subplot(2,2,1);
    plot(h(1,:),h(2,:),'--r')
    hold all
    x =hpoints(:,1);
    y = hpoints(:,2);
    scatter(x,y,150,'rd','fill');
end

% ----------------------
% Model the material volume composition
%----------------------------
%deltaRSize = 1/10; % Resuse, the other one

% Vol Fraction models the Sn in the fly wheel
%global VolFractpoints;
% VolFractpoints = [innerR,0.7;
%             0.04,0.7;
%              0.08,0.7;
%               0.12,0.7;
%              outerR,0.7]
%t = linspace(0,1,1/deltaRSize) % Resuse, the other one
volFrac =  bezierInter(VolFractpoints,t)

% ----------------------
% Plot the volume and control poiints
%----------------------------
%subplot(2,2,1);
if(doplot ==1)
    plot(volFrac(1,:),volFrac(2,:),'--g')
    % hold all
    x =VolFractpoints(:,1);
    y = VolFractpoints(:,2);
    scatter(x,y,150,'gd','fill');
    t= strcat('Green is Vol Fraction of Sn (heavy), Red is height');
    title(t);

    hold off
end