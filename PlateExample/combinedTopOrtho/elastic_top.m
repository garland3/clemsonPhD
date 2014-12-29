function elastic_top()

actualTopElastic(60,30,0.5,3,2)
%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function [x] =  actualTopElastic(nelx,nely,volfrac,penal,rmin)
% recvid=1;       %turn on or off the video recorder
% %% FEA and Elastic problem initialization
% if recvid==1
%     vidObj = VideoWriter('results_homog_level_set.avi');    %Prepare the new file for video
%     vidObj.FrameRate = 50;
%     vidObj.Quality = 100;
%     open(vidObj);
%     vid=1;
% end
E1 = 2.7e9; % Young's mod
E2 = E1*.8;  % 80 percent of E1
v1 = 0.3; % Piossons ratio
v2 = 0.25; % Piossons ratio

avgE = (E1+E2)/2;
avgV = (v1+v2)/2;
G = avgE/(2*(1+avgV));

nn = (nelx+1)*(nely+1); % number of nodes
ne = nelx*nely; % number of elements
count = 1;
numNodesInRow = nelx+1;
numNodesInColumn = nely+1;
IEN = zeros(nn,4); % Index of element nodes (IEN)
% Each row, so nely # of row
for i = 1:nely
     rowMultiplier = i-1;
    % Each column, so nelx # of row
    for j= 1:nelx        
        IEN(count,:)=[rowMultiplier*numNodesInRow+j, ...
                      rowMultiplier*numNodesInRow+j+1, ...
                      (rowMultiplier +1)*numNodesInRow+j+1, ...
                       (rowMultiplier +1)*numNodesInRow+j];
        count = count+1;
    end
end


% INITIALIZE
x(1:nely,1:nelx) = volfrac; % density of the elements
rotM = ones(nely,nelx)*pi/2; % rotation matrix
loop = 0; 
change = 1.;
% START ITERATION
while change > 0.01  
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U,rotationResultMatrix]=FEA_elastic_V2(nelx,nely,x,penal,rotM,E1,E2,v1,v2,G);
  
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
%  [KE] = elK_elastic;
  
  
elementsInRow = nelx+1;
c = 0.; % c is the objective. Total strain energy
% % loop over the elements
elx = 1;
ely = 1;
for e = 1:ne   

         nodes1 = IEN(e,:);
         xDOF = nodes1*2-1;
         yDOF = nodes1*2;

         % I cannot use the union, or else the order get messed up. The order
         % is important. Same in the actual topology code when you are
         % calculating the objectiv
          dofNumbers = [xDOF(1) yDOF(1) xDOF(2) yDOF(2) xDOF(3) yDOF(3) xDOF(4) yDOF(4)];
    
         % NodeNumbers = union(xNodeNumbers,yNodeNumbers);
        
         Ue = transpose(U(1,dofNumbers));
         
         theta = rotM(ely,elx);
        KE =  ke_elementV2(theta,E1,E2,v1,v2,G );
         
         c = c + x(ely,elx)^penal*Ue'*KE*Ue;
          
         % for the x location
         % The first number is the row - "y value"
         % The second number is the column "x value"
         
         %  dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; % objective sensitivity, partial of c with respect to x
         dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; % objective sensitivity, partial of c with respect to x
         elx = elx+1;
        if(elx>nelx)
          elx = 1;
          ely = ely+1;
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
  
 % rotM = convertRotMResultToMatrix(rotationResultMatrix,nelx,nely,ne);
  
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

function rotM = convertRotMResultToMatrix(rotationResultMatrix,nelx,nely,ne)
xLoc = 1;
yLoc = 1;
rotM = zeros(nely,nelx); % rotation matrix
for e = 1:ne
    rotM(yLoc,xLoc) = rotationResultMatrix(e);
      xLoc = xLoc+1;
      if(xLoc>nelx)
          xLoc = 1;
          yLoc = yLoc+1;
      end
 end

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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