function [G_23,U3,Comp,dComp,dG] = trans_shear_analysis(coord, D, le, lm, ne, nn, nny, nnx, nelx, nely, id, nsdof, x,Dim,penal,K,He,Bm)

%%% FINITE ELEMENT ANALYSIS WITH TRANSVERSE SHEAR BOUNDARY CONDITIONS

strain_def = 0.1; % GLOBAL STRAIN

% for numerical integration
a    = 1/sqrt(3);
xint = [-a -a; a -a; a  a; -a  a];

% Locate Boundary Nodes
allnodes  = (1:nn); % all nodes
tnodes    = allnodes( coord(2,:) == max(coord(2,:)));
bnodes    = allnodes( coord(2,:) == min(coord(2,:)));
lnodes    = allnodes( coord(1,:) == min(coord(1,:)));
rnodes    = allnodes( coord(1,:) == max(coord(1,:)));

% Set Degrees of Freedom
tudofs = id(1,tnodes); tvdofs = id(2,tnodes);
budofs = id(1,bnodes); bvdofs = id(2,bnodes);
ludofs = id(1,lnodes); lvdofs = id(2,lnodes);
rudofs = id(1,rnodes); rvdofs = id(2,rnodes);

% All Degrees of Freedom
alldofs  = (1:1:2*nn);

%%%%%%%%%%%%%%%%%%%%%%%% Loading Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constraint_method = 'lagrange';

switch lower(constraint_method)
    case{'lagrange','penalty','transformation'}

        % Transverse Shear (Section 2.2 from Drago and Pindera)
        constraintdofs3 = [];
        %        constraintdofs3(:,2) = [];
        displdofs3(:,1) = [budofs bvdofs rudofs rvdofs tudofs tvdofs ludofs lvdofs];
        n = 1;
        %S_1
        for i = 1:length(budofs)
            displdofs3(n,2) = coord(2,bnodes(i))*strain_def;
            n = n+1;
        end
        for i = 1:length(budofs)
            displdofs3(n,2) = coord(1,bnodes(i))*strain_def;
            n = n+1;
        end
        %S_2
        for i = 1:length(rudofs)
            displdofs3(n,2) = coord(2,rnodes(i))*strain_def;
            n = n+1;
        end
        for i = 1:length(rvdofs)
            displdofs3(n,2) = coord(1,rnodes(i))*strain_def;
            n = n+1;
        end
        %S_3
        for i = 1:length(tudofs)
            displdofs3(n,2) = coord(2,tnodes(i))*strain_def;
            n = n+1;
        end
        for i = 1:length(tvdofs)
            displdofs3(n,2) = coord(1,tnodes(i))*strain_def;
            n = n+1;
        end
        %S_4
        for i = 1:length(ludofs)
            displdofs3(n,2) = coord(2,lnodes(i))*strain_def;
            n = n+1;
        end
        for i = 1:length(lvdofs)
            displdofs3(n,2) = coord(1,lnodes(i))*strain_def;
            n = n+1;
        end

        fixdofs3 = [];
        forcedofs3 = [];
end

freedofs3  = setdiff(alldofs,fixdofs3);

[U3, F3, stress, strain, Comp,dComp,Ce] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
    a, xint, freedofs3, fixdofs3, x, penal, ...
    constraint_method, constraintdofs3, forcedofs3, displdofs3, alldofs,K,He,Bm);

L = Dim(1); W = Dim(2); T = Dim(3);
Vol = L*W*T;
avg_stress = mean(stress,1);
avg_strain = mean(strain,1);
G_23 = avg_stress(1,3)/avg_strain(1,3);
% disp(' ')
% disp(['G_23 = ', num2str(G_23)])
dComp=penal*Ce./x*T;
dG=dComp/Vol/(avg_strain(1,3))^2/(ne^2);
% disp(' ')
% disp(['G_23 = ', num2str(G_23)])

%%%%%%%%%%%%%%%%%%%%%%%% Face Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function facenodes = getfacesnodes(le, element, face) %#ok<DEFNU>
%                3
%         ----------------
%         |              |
%       4 |              | 2
%         |              |
%         ----------------
%                1
if face == 1
    facenodes = le([1,2], element);
elseif face == 2
    facenodes = le([2,3], element);
elseif face == 3
    facenodes = le([3,4], element);
elseif face == 4
    facenodes = le([4,1], element);
else
    facenodes = [];
end
%%%%%%%%%%%%%%%%%%%%%%%% Stiffness Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ke, Re, Be, detJ] = make_elemstiff(xi, eta, local_coord, D)
%Stiffness Matrix evaluate it this way since derivates become complicated
%for rotated material properties
N1_xi  = -1/4*(1-eta);   N2_xi  =  1/4*(1-eta);  N3_xi  =  1/4*(1+eta);   N4_xi  = -1/4*(1+eta);
N1_eta = -1/4*(1-xi);    N2_eta = -1/4*(1+xi);   N3_eta =  1/4*(1+xi);    N4_eta =  1/4*(1-xi);

der_xi_eta = [N1_xi  N2_xi  N3_xi  N4_xi; N1_eta N2_eta N3_eta N4_eta;];
Jacobian   = der_xi_eta*local_coord';
Jinv       = inv(Jacobian);
detJ       = det(Jacobian);
dxy        = Jinv*der_xi_eta;

B_xy       = [dxy(1,1) 0        dxy(1,2) 0        dxy(1,3) 0        dxy(1,4) 0 ;
    0        dxy(2,1) 0        dxy(2,2) 0        dxy(2,3) 0        dxy(2,4);
    dxy(2,1) dxy(1,1) dxy(2,2) dxy(1,2) dxy(2,3) dxy(1,3) dxy(2,4) dxy(1,4)];
Ke         = B_xy'*D*B_xy*detJ;
Re         = B_xy'*D*detJ; % For Homogenization Purposes
Be         = B_xy*detJ;
%%%%%%%%%%%%%%%%%%%%%%%% FEA Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U, F, stress, strain, Comp, dU,Ce] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
    a, xint, freedofs, fixdofs, x, penal, ...
    constraint_method, consdofs, forcedofs, displdofs, alldofs,K,He,Bm)

R  = sparse(nsdof*nn, 1);
U = sparse(nsdof*nn, 1);

switch lower(constraint_method)
    case('lagrange')
        warning off % error is normal for Lagrange Multiplier see pg 432 Cook
        [mm n_fdofs] = size(fixdofs);
        [pp n_adofs] = size(alldofs);
        [n_consdofs,qq] = size(consdofs);
        displdofs = unique(displdofs, 'rows');
        [n_displdofs,qqs] = size(displdofs);
        %         for i = 1:n_displdofs
        %             for j = i+1:n_displdofs
        %                 if displdofs(i,1) == displdofs(j,1);
        %                     displdofs(j,:) = [];
        %                     n_displdofs = n_displdofs-1;
        %                 end
        %             end
        %         end

        %Build constraint matrix
        C = sparse(n_fdofs+n_consdofs+n_displdofs,n_adofs);
        Q = sparse(n_fdofs+n_consdofs+n_displdofs,1);
        n = 1;
        local_cons = [1 -1];
        for e = 1:n_consdofs
            C(n, consdofs(e,:)) = local_cons;
            n = n+1;
        end
        for e = 1:n_fdofs
            C(n,fixdofs(e)) = 1;
            n = n+1;
        end
        for e = 1:n_displdofs
            C(n,displdofs(e,1)) = 1;
            Q(n,1) = displdofs(e,2);
            n = n+1;
        end
        c_dofs = [displdofs(:,1)'];
        [pp n_cdofs] = size(c_dofs);

        lagrange_dofs   = (nsdof*nn+1: nsdof*nn+n_cdofs);
        zeromatrix      = sparse(n_cdofs, n_cdofs);
        zerovec         = sparse(n_cdofs, 1);

        KT              = [K C'; C zeromatrix];
        UT             = sparse(nsdof*nn + n_cdofs, 1);
        RT              = [R; Q];

        UT = KT\RT;

        for i = 1:nsdof*nn
            U(i,1) = UT(i,1);
        end

        %         cond_K = cond(K);
%                 cond_KT = cond(KT);
        F = K*U;

end
Comp=U'*F;
strain=zeros(ne,3);
stress=zeros(ne,3);
for e = 1:ne
    strain(e,:) = Bm*U(lm(:,e));
    stress(e,:) = D*x(e)^(penal)*strain(e,:)';
end

Ce=sparse(length(x),1);
for e=1:length(x)
    GlobalDof   = lm(:,e);    
    Ce(e)=U(GlobalDof)'*He(:,:,e)*U(GlobalDof);
end
dU=Ce;

% dd=c_dofs;
% ff=setdiff(alldofs,c_dofs);
% dUg=zeros(length(U),ne);
% Kffinv=inv(K(ff,ff));
% for e=1:ne
%     i = le(:,e);
%     local_coord = coord(:,i);
%     Km = 0;
%     for int = 1:4
%         [Ke, Re, Be, detJ] = make_elemstiff(xint(int,1), xint(int,2), local_coord, D);
%         Km = Km + Ke;
%         Bm = Bm + Be;
%     end
%     GlobalDof   = lm(:,e);
%     H = sparse(nsdof*nn,nsdof*nn);
%     H(GlobalDof,GlobalDof) = x(e)^penal*Km;
%     dUg(ff,e)=-Kffinv*(H(ff,ff)*U(ff)+H(ff,dd)*U(dd))*penal/x(e);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%MESH Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcnew]= meshfilter(dc, x, le, ne, coord, rmin)
%Mesh Filter Taken from Ole. Sigmund
%Programmmed a little differently
dcnew  = zeros(ne,1);
xcoord = coord(1,le([ 1 2],[1:ne])); %Tricky it places the first two x coordinates element
ycoord = coord(2,le([ 1 4],[1:ne]));
xbase  = (xcoord([1:2:2*ne]) + xcoord([2:2:2*ne]))/2;
ybase  = (ycoord([1:2:2*ne]) + ycoord([2:2:2*ne]))/2;

dist  = zeros(ne,1);
sum   = zeros(1);
% for f= 1:ne
%      dist   =  sqrt((xbase(f) -xbase).^2 + (ybase(f) -ybase).^2);
%      int    =  (rmin -dist);  int(int>rmin)= 0;
%      sum    =  sum + max(0, int )';
%      dcnew(:,1) = dcnew(:,1) + max(0, rmin -dist)'*x(f)*dc(f,1);
% end
for f= 1:ne
    dist   =  sqrt((xbase(f) -xbase).^2 + (ybase(f) -ybase).^2);
    int    =  (rmin -dist);  int(dist>rmin)= 0;
    sum    =  sum + int';
    dcnew(:,1) = dcnew(:,1) + max(0, rmin -dist)'*x(f)*dc(f,1);
end
dcnew(:,1) = dcnew(:,1)./(x(:).*sum);
%%%%%%%%%%%%%%%%%%%% Objective and Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,df] = fob(x,ne)

f = sum(x)/ne;
df = 1/ne*ones(ne,1);

%%%%%%%%%%%%%%%%%%%% Constraints and Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,h,dg,dh]=gh(x,coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
    a, xint, freedofs3, fixdofs3, penal, ...
    constraint_method, constraintdofs3, forcedofs3, displdofs3, alldofs,Ds,Dim,rmin,dx,dy)
global G_23 stress 
[U3, F3, stress, strain, Comp,dU,Bm,dUg] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
    a, xint, freedofs3, fixdofs3, x, penal, ...
    constraint_method, constraintdofs3, forcedofs3, displdofs3, alldofs);
% f=-Comp;
L=Dim(1);
W=Dim(2);
T=Dim(3);
Vol=L*W*T;
% Comp_target=0.5*Ds(3,3)**mean(strain(:,3))^2;
% keyboard
avg_stress = mean(stress(:,3));
avg_strain = mean(strain(:,3));
G_23 = avg_stress/avg_strain;
Ce=dU;
dC=penal*Ce./x*T;
dG=dC/Vol/avg_strain^2/(ne^2);

% filter = fspecial('gaussian', [3 3],1);
% dGf = filter2(filter,dG);
%dGf = meshfilter(dG, x, le, ne, coord, rmin); 

% figure
% subplot(211),plot(1:length(x),[dG dGf]);
% subplot(212),plot(1:length(x),[x x.^3]);
% pause(1)
% close 

inequal_tol = 10^-6;

f = (G_23-Ds(3,3))^2;
df = 2*(G_23-Ds(3,3))*dG;

df = meshfilter(df, x, le, ne, coord, rmin);

f1 = f-inequal_tol;
df1 = df;

% f = (G_23-Ds(3,3))/D(3,3);
% df = dG/D(3,3);
%von_mises = sqrt(stress(:,1).^2+stress(:,2).^2-stress(:,1).*stress(:,1)+3*stress(:,3).^2);
%sum_VM = sum(von_mises,1);

disp(G_23)

g=[f1];
h=[];
dg = [df1];
dh = [];






