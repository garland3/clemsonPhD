function trans_shear_multi_cell_non_pneu_2_1
%(refine,rmin_init,penal,W,ratio,target,x_cells,y_cells,mat_flag)

% EXAMPLE DRIVER: SET TO RUN STANDALONE IN ITS CURRENT FORM
% INPUT PARAMETERS:
% refine: unit cell mesh refinement - refine = 40 will create a minimum
% refinement of 40 elements per unit cell side
% rmin_init: mesh filter radius - rmin_init = 2 will smear the values at
% each step in the analysis by a radius of 2 element's side length
% penal: SIMP penalization factor (use 3 or 4)
% ratio: defines the ratio of the unit cell's geometry according to L = W/ratio;
% target: shear modulus target
% x_cells: number of unit cells in x-direction
% y_cells: number of unit cells in y-direction
% mat_flag: sets which base material properties to use


refine = 40; rmin_init = 2; penal = 4; height = 9; ratio = 2; target = 7.444; x_cells = 10; y_cells = 1; mat_flag = 1;

%design params
%x_cells = 2; y_cells = 1;
start = 2; %1- full, 2- sq. hole, 3-plus, 4- x, 5- rand (DON'T USE 1)
%refine = 10;
rmin_des = 2;
%rmin_init = 1.5;
%L = 1;
%W = 1;
T = 200; % plane stress thickness = 200mm or 0.2m
%target = 4; (units are MPa)
print_case = 'on';
if mat_flag == 1
    material = 'PC';
elseif mat_flag == 2
    material = 'steel';
elseif mat_flag == 3
    material = 'alum_7075-T6';
end
% Target shear moduli will range between 6 and 14 MPa
% Target W correspond to a range between 5 and 12 mm according to eq. G*W = 67.

%opt params
%rmin = 1.5;
% bounds on allowable finite element densities for topology optimization
% (decision variable bounds)
xmin = 0.001; 
xmax = 1;
volfrac = 1; %upper bound on objective function
%penal = 3;

%Material Properities
%Plane Stress (Units are of stress are MPa)
switch material
    case 'PC'
        E = 2700; v = 0.42;
    case 'steel'
        E = 210000; v = 0.29;
    case 'alum_7075-T6'
        E = 70000; v = 0.33;
end
% defines the base material tensor for the FEA code
D = E/(1-v^2)* [1 v 0; v 1 0; 0 0 0.5*(1-v)];
%D = [30 10 0; 10 30 0; 0 0 10];
% Desired Material Properties
Ds = [10 1.3716 0.00; 1.3716 10 0.00; 0.00 0.00 target];

tic %timer on

% set up parameters based on inputs
W = height/y_cells; % unit cell height
L = W/ratio; % unit cell length
Dim(1) = L; Dim(2) = W; Dim(3) = T; % place in vector array for FEA
nelx = refine; % number of elements in x-direction (per unit cell)
nely = refine*ratio; % number of elements in y-direction (per unit cell)
dx = L/nelx; % finite element dx
dy = W/nely; % finite element dy
rmin = rmin_init*dx; % scale the actual mesh filter radius to appropriate resolution
rmin_des = rmin_des*dx; 

% FINITE ELEMENT INTIALIZATION INFORMATION
nsdof = 2; % number of spatial dof
nen = 4; % number of element nodes
ne = nelx*nely; % total number of elements
xo = -L/2; 
yo = -W/2;
[coord, nnx, nny] = make_coord(xo,yo,L,W, nelx, nely); %Generate mesh
le = make_nodalconn(nelx, nely, nnx, nny, nen, ne);%Generate connectivity
[id,lm] = make_lm_and_id(le, nsdof, ne, nen, nnx, nny);
nn = nnx*nny;

% for numerical integration
a    = 1/sqrt(3);
xint = [-a -a; a -a; a  a; -a  a];

%Initialize density vector (initial point for optimization)
if start == 1
    initial_point = 'full';
elseif start == 2
    initial_point = 'kikuchi';
elseif start ==3
    initial_point = 'plus_init';
elseif start  == 4
    initial_point = 'x_init';
elseif start == 5
    initial_point = 'x_rand';
end
switch lower(initial_point)
    case('full')
        x = ones(ne,1)*xmax;  
    case('kikuchi') %test case for analysis (must be a square cell)
        n = 1;
        for i = 1:refine*ratio/2-3
            for j = 1:refine
                x(n,1) = xmax;
                n = n+1;
            end
        end
        for i = refine*ratio/2-2:refine*ratio/2+3
            for j = 1:refine/2-3
                x(n,1) = xmax;
                n = n+1;
            end
            for j = refine/2-2:refine/2+3
                x(n,1) = xmin;
                n = n+1;
            end
            for j = refine/2+4:refine
                x(n,1) = xmax;
                n = n+1;
            end
        end
        for i = refine*ratio/2+4:refine*ratio
            for j = 1:refine
                x(n,1) = xmax;
                n = n+1;
            end
        end
    case('x_init')
        x = x_xinit(xmin, xmax, nelx, nely, L, W, le, nen, ne, rmin_des,coord,volfrac);
    case('plus_init')
        x = plus_xinit(xmin, xmax, nelx, nely, L, W, le, nen, ne, rmin_des,coord,volfrac);    
    case('x_rand')
        x = make_xinit(nelx,nely,ne,displdofs',[],xmin,xmax,volfrac,lm,sizeforce,3);
end

% figure
% for e = 1:ne
%     i      = le(:,e);
%     xcoord = coord(1,i);
%     ycoord = coord(2,i);
%     patch(xcoord, ycoord,x(e))
% end
% colormap(flipud(gray)); axis tight; axis equal; axis off

xsym = zeros(nely,ne);
ysym = zeros(nelx,ne);
mm   = 0;
mmm  = 0;
jj   = 1;

%Enforce xy-symmetry
xsym = [];
ysym = [];
for e = 1:ne
    i      = le(:,e);
    xcoord = coord(1,i);
    ycoord = coord(2,i);
%    patch(xcoord, ycoord,x(e))
    for ee = e+1:ne
        ii      = le(:,ee);
        xcoorde = coord(1,ii);
        ycoorde = coord(2,ii);
        if (abs(xcoord(1)+xcoorde(2))<=10^-5) && (abs(ycoord(1)-ycoorde(1))<=10^-5)
            mm = mm + 1;
            xsym(mm,e)  =  1;
            xsym(mm,ee) = -1; 
        end  
        if (abs(ycoord(1)+ycoorde(4))<=10^-5) && (abs(xcoord(1)-xcoorde(4))<=10^-5)
            mmm = mmm + 1;
            ysym(mmm,e)  =  1;
            ysym(mmm,ee) = -1; 
        end  
    end
    hold on
%    plot(coord(1,[rnodes  tnodes]), coord(2,[rnodes  tnodes]),'.r')
end
disp('Symmetry Enforced')
n_x = x_cells;
n_y = y_cells;

E_22 = zeros(n_x,n_y);
E_33 = zeros(n_x,n_y);
Q_23 = zeros(n_x,n_y);
Q_32 = zeros(n_x,n_y);
G_23 = zeros(n_x,n_y);

y_cells = n_y;
x_cells = n_x;

uc_x = x_cells;
uc_y = y_cells;

%input parameters
n_uc = uc_x*uc_y;
g_nelx = uc_x*nelx;
g_nely = uc_y*nely;

g_ne = uc_x*nelx*uc_y*nely; %total number of global elements
g_xo = -L/2*uc_x;
g_yo = -W/2*uc_y;
[g_coord, g_nnx, g_nny] = make_coord(g_xo,g_yo,L*uc_x,W*uc_y, g_nelx, g_nely); %Generate Global mesh
g_le = make_nodalconn(g_nelx, g_nely, g_nnx, g_nny, nen, g_ne);%Generate Global connectivity
[g_id,g_lm] = make_lm_and_id(g_le, nsdof, g_ne, nen, g_nnx, g_nny);
g_nn = g_nnx*g_nny; %number of global nodes
lse = make_se_conn(g_nelx,g_nely,g_ne,nelx,nely,uc_x,uc_y,n_uc);

% Locate Boundary Nodes
allnodes  = (1:g_nn); % all nodes
tnodes    = allnodes(g_coord(2,:) == max(g_coord(2,:)));
bnodes    = allnodes(g_coord(2,:) == min(g_coord(2,:)));
lnodes    = allnodes(g_coord(1,:) == min(g_coord(1,:)));
rnodes    = allnodes(g_coord(1,:) == max(g_coord(1,:)));

% Set Degrees of Freedom
tudofs = g_id(1,tnodes); tvdofs = g_id(2,tnodes);
budofs = g_id(1,bnodes); bvdofs = g_id(2,bnodes);
ludofs = g_id(1,lnodes); lvdofs = g_id(2,lnodes);
rudofs = g_id(1,rnodes); rvdofs = g_id(2,rnodes);

% All Degrees of Freedom
alldofs  = (1:1:2*g_nn);

%%% Given a local solution, generate a global solution and analyze the
%%% solution in shear.
g_x = loc_to_glob(x,lse,n_uc,ne);

% Plot Undeformed Position
% h = figure;
% for e = 1: g_ne
%     i      = g_le(:,e);
%     xcoord = g_coord(1,i);
%     ycoord = g_coord(2,i);
%     patch(xcoord, ycoord,g_x(e))
% end
% colormap(flipud(gray)); axis tight; axis equal; axis off
% %name = sprintf('%7s%1g%1s%1g','undef_sol',uc_y,'x',uc_x);
% %print(h,'-djpeg',name);
% close(h)

%title(['Design Domain elements in X direction ', num2str(nelx), ',  elements in Y direction ', num2str(nely)]) 
g_Dim = [Dim(1)*uc_x Dim(2)*uc_y Dim(3)];
g_Vol = g_Dim(1)*g_Dim(2)*g_Dim(3);

%Optimization Algorithm
%Design Variable Bounds
lb=xmin*ones(ne,1);
ub=xmax*ones(ne,1);

%Inequality Constraints
A = [];
B = [];
% Volume Constraint
% A = [ones(1,ne)/ne];
% B = [volfrac];

%Equality Constraints
% Symmetry Constraints
Aeq = [];
Beq = [];
Aeq = [Aeq; xsym; ysym];
Beq = [Beq; zeros(mm+mmm,1)];
% MATLAB OPTIMIZATION PARAMETERS
options = optimset('Algorithm','active-set','Display','iter',...
    'TolX',1e-3,'TolFun',1e-4,'TolCon',1e-3','MaxIter',400,'DerivativeCheck','off','GradObj','On','GradConstr','on','MaxFunEvals',1000);
%options = optimset('Algorithm','active-set','Display','iter',...
    %'TolX',1e-8,'TolFun',1e-6,'MaxIter',200);
%  [x,f]=fmincon_mod(@(x) fob(x, coord, D, le, lm, ne, nn, nny, nnx, nsdof,...
%      a, xint, freedofs3, fixdofs3, penal, ...
%      constraint_method, constraintdofs3, forcedofs3, displdofs3, alldofs,Ds),x,ones(1,ne)/ne,volfrac,[],[],lb,ub,@(x) gh(x,Ds),options);

% RUN OPTIMIZER (THE CODE WILL SPEND THE MAJORITY OF ITS TIME ON THIS LINE
[x,f,exitflag,output,lambda]=fmincon(@(x)fob(x,ne),x,A,B,Aeq,Beq,lb,ub,@(x)gh(x,coord, D, le, ne, nsdof,penal, Ds,Dim,rmin,g_coord, g_le, g_lm, g_ne, g_nn, g_nny, g_nnx, g_nelx, g_nely, g_id, g_Dim,lse,n_uc),options);

% POST-PROCESSING ANALYSIS CALLS
% MAP OPT VALUES TO GLOBAL GEOM VECTOR
g_x = loc_to_glob(x,lse,n_uc,ne);

% GENERATE STIFFNESS MATRIX AND PERFORM FEA ONE FINAL TIME
[K,He,Bm] = make_stiff(g_coord, D, g_le, g_lm, g_ne, g_nn, nsdof,xint, g_x, penal);
[G_23,U3,Comp3,dComp3,dG_23] = trans_shear_analysis(g_coord, D, g_le, g_lm, g_ne, g_nn, g_nny, g_nnx, g_nelx, g_nely, g_id, nsdof, g_x,g_Dim,penal,K,He,Bm);

disp(' ')
disp([G_23 Ds(3,3)])

t = toc; % STOP TIMER

% RETRIEVE OPT DATA FROM FINAL OPT VALUES
iters = output.iterations;
func_evals = output.funcCount;

disp(['time = ',num2str(t)])

% figure
% for e = 1: g_ne
%     i      = g_le(:,e);
%     g_xcoord = g_coord(1,i);
%     g_ycoord = g_coord(2,i);
%     patch(g_xcoord, g_ycoord,g_x(e))
% end
% colormap(flipud(gray)); axis tight; axis equal; axis off

% PRINT OUTPUT 
switch print_case
    case 'on'
        name = (['tweel-refine=',num2str(refine),'_target=',num2str(target),'_penal=',num2str(penal),'_ratio=',num2str(ratio),'_rmin=',num2str(rmin_init),'_start=',num2str(start),'x_cells=',num2str(x_cells),'y_cells=',num2str(y_cells),'.txt']);
                 
        fout = fopen(name, 'w'); %print output to two different files

        fprintf(fout,'Material Tensor = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f\t%-5.3f\t%-5.3f',D(1,1),D(1,2),D(1,3));
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f\t%-5.3f\t%-5.3f',D(2,1),D(2,2),D(2,3));
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f\t%-5.3f\t%-5.3f',D(3,1),D(3,2),D(3,3));
        fprintf(fout,'\n');

        fprintf(fout,'Desired Modulus = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f',target);
        fprintf(fout,'\n');

        fprintf(fout,'Volume Filled = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-4.3f',sum(x)/ne);
        fprintf(fout,'\n');

        fprintf(fout,'Volume Upper Bound = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-4.3f',volfrac);
        fprintf(fout,'\n');

        fprintf(fout,'Domain Size: ');
        fprintf(fout,'\n');
        fprintf(fout,'nelx = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-2.0f',nelx);
        fprintf(fout,'\n');
        fprintf(fout,'nely = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-2.0f',nely);
        fprintf(fout,'\n');

        fprintf(fout,'Computation Time = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-8.7f',t);
        fprintf(fout,'\n');

        fprintf(fout,'G_23 = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-6.4f',G_23);
        fprintf(fout,'\n');
        
        fprintf(fout,'iterations = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-6.4f',iters);
        fprintf(fout,'\n');

        fprintf(fout,'functions_evals = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-6.4f',func_evals);
        fprintf(fout,'\n');

        fprintf(fout,'Termination = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f',exitflag);
        fprintf(fout,'\n');
        
        fprintf(fout,'initial point = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f',start);
        fprintf(fout,'\n');
        
        fprintf(fout,'rmin_initial point = ');
        fprintf(fout,'\n');
        fprintf(fout,'%-5.3f',rmin_des);
        fprintf(fout,'\n');
       
        fprintf(fout,'Solution = ');
        fprintf(fout,'\n');
        for e = 1:ne
            fprintf(fout,'%-5.4f',x(e));
            fprintf(fout,'\n');
        end
        fclose(fout);
end
%%%%%%%%%%%% Connectivity between local and global analysis %%%%%%%%%%%%%%%
% FUNCTIONS PROVIDE MAPPING FROM UNIT CELL TO FULL GEOMETRY AND BACK
function lse = make_se_conn(nelx,nely,ne,loc_nelx,loc_nely,nselx,nsely,nse)
lse = zeros(loc_nelx*loc_nely,nse);
e = 0;
for se_y = 1:nsely
    for j = 1:loc_nely 
        for se_x = 1:nselx
            for i = 1:loc_nelx
                e = e+1;
                lse(i+(j-1)*loc_nelx,se_x+(se_y-1)*nselx) = e;
            end
        end
    end
end    
%%%%%%%%%%%%%%%% Copy local densities to global density vector %%%%%%%%%%%%
% FUNCTIONS PROVIDE MAPPING FROM UNIT CELL TO FULL GEOMETRY AND BACK
function x = loc_to_glob(loc_x,lse,nse,loc_ne)
for se = 1:nse
    for i = 1:loc_ne
        el = lse(i,se);
        x(el,1) = loc_x(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%Generate Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coord, nnx, nny] = make_coord(xo, yo, L, W, nelx, nely)
dx   = L/nelx;    %x spacing of element
dy   = W/nely;    %y spacing of element
nel  = nelx*nely; %total nelx
nnx  = nelx+1;    %number of nodes in x
nny  = nely+1;    %number of nodes in y
n    =  0;        %sequential counter 
for j = 1:nny
    for i = 1:nnx
        n = n + 1;
        coord(1,n) = xo +dx*(i-1);
        coord(2,n) = yo+ dy*(j-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Le Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function le = make_nodalconn(nelx, nely, nnx, nny, nen, ne)

%first element nodal connectivity
le = zeros(nen, ne);
le(:,1)= [ 1 2 (nnx+1)+1 nnx+1]';

for j = 2:nelx
    le(:,j) = le(:,j-1) + 1;
end

n = nelx;
for j = 1:nely - 1
    for i = 1:nelx
        n = n + 1;
        le(:,n) = le(:,n-j*nelx ) + j*nnx;    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%% Id and LM Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id, lm] = make_lm_and_id(le, nsdof, ne, nen, nnx, nny)

%form id and lm arrays
for i = 1:nnx*nny; %for each node
    for a = 1:nsdof; %for each dof 
        id(a,i) = nsdof*(i-1) + a; %creates a nodal id array (global)
    end
end

for e = 1:ne %for each element
    for i = 1:nen %for each node in each element
        for a = 1:nsdof %for each dof (1 or 2)
            p = nsdof*(i-1) + a; %p = 1:8
            lm(p,e) = id(a, le(i,e)); % for each element, lm(1:8, array inside of array)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%% Face Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function facenodes = getfacesnodes(le, element, face) %#ok<DEFNU>
% FINDS FACE NODES OF THE UNIT CELL
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
function [K,He,Bm] = make_stiff(coord, D, le, lm, ne, nn, nsdof,xint, x, penal)
% MAKE FULL STIFFNESS MATRIX
K  = sparse(nsdof*nn, nsdof*nn);

He = zeros(8,8,length(x));
for e = 1:ne
    i = le(:,e);
    local_coord = coord(:,i);
    Km = 0;
    Bm = 0;
    for int = 1:4
        [Ke, Re, Be, detJ] = make_elemstiff(xint(int,1), xint(int,2), local_coord, D);
        Km = Km + Ke;
        Bm = Bm + Be;
    end
    GlobalDof   = lm(:,e);
    K(GlobalDof, GlobalDof) = K(GlobalDof, GlobalDof) + x(e)^penal*Km;
    He(:,:,e)=x(e)^penal*Km;
%     H(GlobalDof,GlobalDof,e)=x(e)^penal*Km;
end
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
% THIS IS THE VOLUME "CONSTRAINT" NECESSARY FOR T.O, EVEN THOUGH IT IS 
% TREATED AS THE OBJECTIVE FUNCTION HERE
f = sum(x)/ne;
df = 1/ne*ones(ne,1);

%%%%%%%%%%%%%%%%%%%% Constraints and Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,h,dg,dh]=gh(x,coord, D, le, ne, nsdof,penal, Ds,Dim,rmin,g_coord, g_le, g_lm, g_ne, g_nn, g_nny, g_nnx, g_nelx, g_nely, g_id, g_Dim,lse,n_uc)
% THE DESIRED MODULUS IS TARGETED BY TREATING THIS FUNCTION AS AN OPT CONSTRAINT
g_x = loc_to_glob(x,lse,n_uc,ne);

% for numerical integration
a    = 1/sqrt(3);
xint = [-a -a; a -a; a  a; -a  a];

% GENERATE STIFFNESS MATRIX
[K,He,Bm] = make_stiff(g_coord, D, g_le, g_lm, g_ne, g_nn, nsdof,xint, g_x, penal);
% CALL FINITE ELEMENT SOLVER
[G_23,U3,Comp3,dComp3,dG_23] = trans_shear_analysis(g_coord, D, g_le, g_lm, g_ne, g_nn, g_nny, g_nnx, g_nelx, g_nely, g_id, nsdof, g_x,g_Dim,penal,K,He,Bm);

% CALCULATE CONSTRAINT SENSITIVITIES
dG_23_g = zeros(ne,1);
for i = 1:ne    
    dG_23_g(i,1) = sum(dG_23(lse(i,:)))*n_uc^2*Dim(1)^2*Dim(2)^2;
end

%%% FINITE DIFFERENCES CODE FOR VERIFICATION %%%%
% G_23_0 = G_23;
% h=1e-7; dG_23dx = zeros(ne,1);
% for i=1:length(x)
%     xp=x;
%     xp(i)=xp(i)+h;
%     g_xp = loc_to_glob(xp,lse,n_uc,ne);
%     [K,He,Bm] = make_stiff(g_coord, D, g_le, g_lm, g_ne, g_nn, nsdof,xint, g_xp, penal);
%     [G_23_1,U3,Comp3,dComp3,dG_23] = trans_shear_analysis(g_coord, D, g_le, g_lm, g_ne, g_nn, g_nny, g_nnx, g_nelx, g_nely, g_id, nsdof, g_xp,g_Dim,penal,K,He,Bm);
%     dG_23dx(i,1)=(G_23_1-G_23_0)/h;
%     disp(['i = ',num2str(i)]) 
%     %dUdx(:,i)=(U3-U30)/h;
% end
% i=1:length(x);
% figure
% plot(i,dG_23dx,i,dG_23_g),grid
% keyboard
%%% END FD CODE %%%%

% PLACE VALUES INTO RETURN VECTORS
f = (G_23-Ds(3,3));
df = dG_23_g;

% SMEAR VALUES OUT USING MESHFILTER
df = meshfilter(df, x, le, ne, coord, rmin);

f1 = f;
df1 = df;

% FINAL RETURN VALUES
g=[];
h=[f1];
dg = [];
dh = [df1];