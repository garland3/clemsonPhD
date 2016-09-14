function [E_xi,v_23,U,Comp,dComp,dE] = trans_normal_analysis(xi,coord, D, le, lm, ne, nn, nny, nnx, nelx, nely, id, nsdof, x,Dim,penal,K,He,Bm)

%%% FINITE ELEMENT ANALYSIS FOR TRANSVERSE NORMAL LOADING BOUNDARY
%%% CONDITIONS

strain_def = 0.1; % GLOBAL STRAIN
% DIMENSIONS OF UNIT CELL
L=Dim(1);
W=Dim(2);
T=Dim(3);

% for numerical integration
a    = 1/sqrt(3);
xint = [-a -a; a -a; a  a; -a  a];

% Locate Boundary Nodes
allnodes  = (1:nn); % all nodes
tnodes    = allnodes( coord(2,:) == max(coord(2,:)));
bnodes    = allnodes( coord(2,:) == min(coord(2,:)));
lnodes    = allnodes( coord(1,:) == min(coord(1,:)));
rnodes    = allnodes( coord(1,:) == max(coord(1,:)));

% Locate Boundary Elements
allelems = (1:ne); % all elements
telems = allelems(nelx*(nely-1)+1:ne);
belems = allelems(1:nelx);
lelems = allelems(1:nelx:ne);
relems = allelems(nelx:nelx:ne);

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
        
        if xi == 1
        % E_11
        displdofs(:,1) = [budofs rudofs tudofs ludofs];

        n = 1;
        %S_1 
        for i = 1:length(budofs) 
            displdofs(n,2) = coord(1,bnodes(i))*strain_def;
            n = n+1;
        end
        %S_2
        for i = 1:length(rudofs) 
            displdofs(n,2) = coord(1,rnodes(i))*strain_def;
            n = n+1;
        end
        %S_3
        for i = 1:length(tudofs) 
            displdofs(n,2) = coord(1,tnodes(i))*strain_def;
            n = n+1;
        end
        %S_4
        for i = 1:length(ludofs) 
            displdofs(n,2) = coord(1,lnodes(i))*strain_def;
            n = n+1;
        end
        
        conselems = [telems;belems];
        consdofs = [2 4 6 8];
        fixdofs = [];
        forcedofs = [];
        
        elseif xi == 2
        % E_22
        displdofs(:,1) = [bvdofs rvdofs tvdofs lvdofs];

        n = 1;
        %S_1 
        for i = 1:length(budofs) 
            displdofs(n,2) = coord(2,bnodes(i))*strain_def;
            n = n+1;
        end
        %S_2
        for i = 1:length(rudofs) 
            displdofs(n,2) = coord(2,rnodes(i))*strain_def;
            n = n+1;
        end
        %S_3
        for i = 1:length(tudofs) 
            displdofs(n,2) = coord(2,tnodes(i))*strain_def;
            n = n+1;
        end
        %S_4
        for i = 1:length(ludofs) 
            displdofs(n,2) = coord(2,lnodes(i))*strain_def;
            n = n+1;
        end
        
        conselems = [lelems;relems];
        consdofs = [1 3 5 7];
        fixdofs = [];
        forcedofs = [];

        end
    
end

freedofs  = setdiff(alldofs,fixdofs);

%E_1 = D(1,1); E_2 = D(2,2); v_12 = D(1,2); v_21 = D(2,1); G = D(3,3);

%C = system(E_1,E_2,G,v_12,v_21);

[U, F, stress, strain, Comp,Ce] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
    a, xint, freedofs, fixdofs, x, penal,constraint_method, consdofs, forcedofs, displdofs, alldofs, conselems,strain_def,xi,Dim,K,He,Bm);

L = Dim(1); W = Dim(2); T = Dim(3);
Vol = L*W*T;
avg_stress = mean(stress,1);
avg_strain = mean(strain,1);
E_xi = avg_stress(1,xi)/avg_strain(1,xi);
% disp(' ')
% disp(['E_',num2str(xi+1) num2str(xi+1),' = ', num2str(E_xi)])
% disp(' ')
if xi == 1
    v_23 = -avg_strain(1,2)/avg_strain(1,1);
%    disp(['v_',num2str(2) num2str(3),' = ', num2str(v_23)])
elseif xi == 2
    v_23 = -avg_strain(1,1)/avg_strain(1,2);
%    disp(['v_',num2str(3) num2str(2),' = ', num2str(v_23)])
end

dComp=penal*Ce./x*T;
dE=dComp/Vol/(avg_strain(1,xi))^2/(ne^2);

%%%% FINITE DIFFERENCES CODE FOR VERIFICATION %%%%
% [U30, F3, stress, strain, Comp, dU] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
%     a, xint, freedofs3, fixdofs3, x, penal, ...
%     constraint_method, constraintdofs3, forcedofs3, displdofs3, alldofs);
% E0=mean(stress(:,xi))/mean(strain(:,xi));
% Ci=dU;
% h=1e-6;
% dUdx=zeros(size(U30),length(x));
% dEdx=zeros(length(x),1);
% for i=1:length(x)
%     xp=x;
%     xp(i)=xp(i)+h;
%     [U3, F3, stress, strain, Comp] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
%     a, xint, freedofs3, fixdofs3, xp, penal, ...
%     constraint_method, constraintdofs3, forcedofs3, displdofs3, alldofs);
%     E1=mean(stress(:,xi))/mean(strain(:,xi));
%     dGdx(i)=(E1-E0)/h;
%     dUdx(:,i)=(U3-U30)/h;
% end
% i=1:length(x);
% %dGdx - FD derivatives, dGdx_2 - analytical derivatives
% dEdx_2=penal/(L*W)/mean(strain(:,xi))^2*Ci./x/(ne^2);
% figure
% plot(i,dEdx,i,dEdx_2),grid
%%%% END FD CODE %%%%

%keyboard

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
function [U, F, stress, strain, Comp,Ce] = fea_solve(coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
    a, xint, freedofs, fixdofs, x, penal,constraint_method, consdofs, forcedofs, displdofs, alldofs,conselems,strain_def,xi,Dim,K,He,Bm)

R  = sparse(nsdof*nn, 1);
U = sparse(nsdof*nn, 1);

switch lower(constraint_method)
    case('lagrange')
        warning off % error is normal for Lagrange Multiplier see pg 432 Cook
        [mm n_fdofs] = size(fixdofs);
        [pp n_adofs] = size(alldofs);
        [qq,n_consdofs] = size(consdofs);
        displdofs = unique(displdofs, 'rows');
        [n_displdofs,qqs] = size(displdofs);
        [n_cons,n_conselems] = size(conselems);
        
        %Build constraint matrix
        C = zeros(n_fdofs+n_displdofs+n_cons,n_adofs);
        Q = zeros(n_fdofs+n_displdofs+n_cons,1);
        n = 1;
        for e = 1:n_fdofs
            C(n,fixdofs(e)) = 1;
            n = n+1;
        end
        for e = 1:n_displdofs
            C(n,displdofs(e,1)) = 1;
            Q(n,1) = displdofs(e,2);
            n = n+1;
        end
        local_cons = ones(1,n_conselems);
        for i = 1:n_cons
            matrix = sparse(8,2*nn);
            m = 1;
            for e = 1:n_consdofs
                for ee = 1:n_conselems
                    vector(ee,1) = lm(consdofs(e),conselems(i,ee)); 
                end    
               matrix(consdofs(m),vector) = local_cons;
               m = m+1;
            end
            if xi == 1
                vector2 = D(2,2)*Bm(2,:)*matrix;
            elseif xi == 2
                vector2 = D(1,1)*Bm(1,:)*matrix;
            end
            C(n,:) = vector2;
            %Q(n,1) = -D(1,2)*strain_def/n_conselems;
            Q(n,1) = -D(1,2)*strain_def/(n_conselems/Dim(xi)^2);
            %Q(n,1) = -D(1,2)*strain_def;
            n = n+1;
        end    
        c_dofs = [displdofs(:,1)'];
        [pp n_cdofs] = size(c_dofs);

        lagrange_dofs   = (nsdof*nn+1: nsdof*nn+n-1);
        zeromatrix      = sparse(n-1,n-1);
        zerovec         = sparse(n-1, 1);

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
