function combinedTopologyOptimization()

actualTopElasticCombined(60,30,0.4,3,2);
%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function [x] =  actualTopElasticCombined(nelx,nely,volfrac,penal,rmin)
% recvid=1;       %turn on or off the video recorder
% %% FEA and Elastic problem initialization
% if recvid==1
%     vidObj = VideoWriter('results_homog_level_set.avi');    %Prepare the new file for video
%     vidObj.FrameRate = 50;
%     vidObj.Quality = 100;
%     open(vidObj);
%     vid=1;
% end

w1 = 0.999; % weight elastic for multi-objective
w2 = 1- w1; % weight heat transfer


% INITIALIZE
x(1:nely,1:nelx) = volfrac; % density of the elements
rotM(1:nely,1:nelx) = zeros; % rotation matrix
loop = 0; 
change = 1.;
% START ITERATION
while change > 0.01  
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=FE_elastic(nelx,nely,x,penal);   
    [U_heatColumn]=temperatureFEA2(nelx,nely,x,penal);   
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = elK_elastic;
   [KEHeat] = lkHeat;
   
     % Normalize the U's for now
          uemax = max(max(U));
          uHeatmax = max(max(U_heatColumn));
          %U = U/uemax;
          %U_heatColumn = U_heatColumn/uHeatmax;
  
  
elementsInRow = nelx+1;
c = 0.; % c is the objective. Total strain energy
for ely = 1:nely
    rowMultiplier = ely-1;
    for elx = 1:nelx
%           n1 = (nely+1)*(elx-1)+ely; 
%           n2 = (nely+1)* elx   +ely;
%           Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
           nodes1=[rowMultiplier*elementsInRow+elx;
                    rowMultiplier*elementsInRow+elx+1;
                    (rowMultiplier +1)*elementsInRow+elx+1;
                    (rowMultiplier +1)*elementsInRow+elx];
        
         
              
            xNodes = nodes1*2-1;
            yNodes = nodes1*2;
    
         % NodeNumbers = union(xNodeNumbers,yNodeNumbers);
          NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    
         % NodeNumbers = union(xNodeNumbers,yNodeNumbers);
          Ue = U(NodeNumbers,:);
          U_heat = U_heatColumn(nodes1,:);
          
        
         
          c = c + w1*x(ely,elx)^penal*Ue'*KE*Ue;
          c = c + w2*x(ely,elx)^penal*U_heat'*KEHeat*U_heat;
          
          % for the x location
          % The first number is the row - "y value"
          % The second number is the column "x value"
          
        %  dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; % objective sensitivity, partial of c with respect to x
        
        
         temp1 = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; % objective sensitivity, partial of c with respect to x
          temp2 = -penal*x(ely,elx)^(penal-1)*U_heat'*KEHeat*U_heat;
           dc(ely,elx) = w1*temp1+w2*temp2;
    end
end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc); 
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES  
subplot(2,2,4);
colormap(winter); 
imagesc(x); 
set(gca,'YDir','normal'); % Flips the image so that row 1 is not at the top
% image(x); 
axis equal; axis tight; axis off;pause(1e-6);
%    F(vid) = getframe; %Get frame of the topology in each iteration
%      writeVideo(vidObj,F(vid)); %Save the topology in the video
%      vid=vid+1;
end 

% if recvid==1
%          close(vidObj);  %close video

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  
%   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
 
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lkHeat

KE= [0.6667   -0.1667   -0.3333   -0.1667;
   -0.1667    0.6667   -0.1667   -0.3333;
   -0.3333   -0.1667    0.6667   -0.1667;
   -0.1667   -0.3333   -0.1667    0.6667];

% function [U]=FE(nelx,nely,x,penal,F,fixeddofs)
% [KE] = lk; 
% K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
% %F = sparse(2*(nely+1)*(nelx+1),1);
% U = zeros(2*(nely+1)*(nelx+1),1);
% for elx = 1:nelx
%   for ely = 1:nely
%     n1 = (nely+1)*(elx-1)+ely; 
%     n2 = (nely+1)* elx   +ely;
%     edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
%     K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
%   end
% end
% % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% % F(2,1) = -1;
% % fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)])
% alldofs     = [1:2*(nely+1)*(nelx+1)];
% freedofs    = setdiff(alldofs,fixeddofs);
% % SOLVING
% U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
% U(fixeddofs,:)= 0;
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [KE]=lk
% 
% KE= [0.6667   -0.1667   -0.3333   -0.1667;
%    -0.1667    0.6667   -0.1667   -0.3333;
%    -0.3333   -0.1667    0.6667   -0.1667;
%    -0.1667   -0.3333   -0.1667    0.6667];