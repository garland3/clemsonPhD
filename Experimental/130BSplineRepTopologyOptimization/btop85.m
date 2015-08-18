% AN 85 LINE CODE FOR TOPOLOGY OPTIMIZATION with B-SPLINES. Xiaoping Qian %%%
function btop85(nelx,nely,volfrac,penal,n_mu, n_mv, mp, mq)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% setting knot intervals and coefficient number
ui = nelx/n_mu; vi = nely/n_mv;
ncp_u = n_mu + mp; ncp_v = n_mv+mq;
loc = reshape(ones(mq+1,1)*(0:mp)*ncp_v+(1:mq+1)'*ones(1,mp+1),1,(mp+1)*(mq+1));
%% Preparing Bspline information
dv = zeros(ncp_u*ncp_v,1);
i = 0;
for elx = 1:nelx
    for ely = 1:nely
     i = i+1;
     ix = floor((elx-0.5)*n_mu/nelx); iy = floor((ely-0.5)*n_mv/nely);
     mElements(i,:) = loc + ix*ncp_v + iy;
     tNus(i,:)=getUniformBsplineBasis(mp, rem(elx-0.5, ui)/ui);    
     tNvs(i,:)=getUniformBsplineBasis(mq, rem(ely-0.5, vi)/vi);   
     Nuvs(i,:) = reshape((tNus(i,:)'*tNvs(i,:))', 1, (mp+1)*(mq+1));
     dv(mElements(i,:)) = dv(mElements(i,:))+ Nuvs(i,:)';
    end
end
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
x = ones(ncp_u*ncp_v,1)*volfrac;
xPhys = reshape(sum(Nuvs.*x(mElements(1:nelx*nely,:)),2),nely,nelx);
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  tdc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  Nuv_tdc=Nuvs.*(reshape(tdc, nelx*nely,1)*ones(1,(mp+1)*(mq+1))); 
  dc = zeros(length(x),1);
  for i = 1:(mp+1)*(mq+1)
      dc(mElements(:,i)) = dc(mElements(:,i))+Nuv_tdc(:,i);
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    xPhys = reshape(sum(Nuvs.*xnew(mElements(1:nelx*nely,:)),2),nely,nelx);
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f \n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
%%%%% B-spline basis from recursive definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [base] = getUniformBsplineBasis(p,t)    
base = 1; z = 0;
ts = t*ones(1,p+1)+base*[0:p];
for d=1:p
  base = (ts((d+1):-1:1).*[z, base]+(d+1-ts((d+1):-1:1)).*[base, z])/d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Xiaoping Qian at the Illinois Institute  %
% of Technology.                                                           %
% Please send your comments to the author: qian@iit.edu                    %
%                                                                          %
% The code is intended for educational purposes and  theoretical details   %
% are discussed in the paper "Topology Optimization in B-spline Space"     %                           %
% by Xiaoping Qian in Computer Methods in Applied Mechanics and            %
% Engineering in press, 2013.                                              %
% This code was modified based on the 88 line of code (www.topopt.dtu.dk)  %
% from Sigmund et al with the goal to illustrate the use of B-splines in   % 
% topology optimization.                                                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

