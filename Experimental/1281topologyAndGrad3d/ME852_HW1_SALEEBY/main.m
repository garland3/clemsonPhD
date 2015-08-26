
%--------------------------------------------------
% read the input files
%--------------------------------------------------
load element_types.dat;
load nodes.dat;        
load elements.dat;
load materials.dat;
load bcstraction.dat;
load bcsforce.dat;
load bcsdisp.dat;

%--------------------------------------------------
% bookkeeping 
%--------------------------------------------------
n_nodes=size(nodes,1);
n_elements=size(elements,1);
n_bcstraction=size(bcstraction,1);
if bcstraction(1,1)==0
    n_bcstraction=0;
end
n_bcsforce=size(bcsforce,1);
n_bcsdisp=size(bcsdisp,1);

fprintf('number of nodes= %d\n', n_nodes);
fprintf('number of elements= %d\n', n_elements);
fprintf('number of traction surfaces= %d\n', n_bcstraction);
fprintf('number of point forces= %d\n', n_bcsforce);
fprintf('number of disp BCS= %d\n', n_bcsdisp);


%--------------------------------------------------
% compute global K matrix
%--------------------------------------------------
K=CompK(nodes, elements, materials);
disp('K complete');


%--------------------------------------------------
% compute global F vector
%--------------------------------------------------
F=CompF(nodes, elements, materials, bcstraction, bcsforce);
disp('F complete');

%--------------------------------------------------
% apply displacement boundary condition
%--------------------------------------------------
coeff=abs(max(K(bcsdisp(1,1)*3-2,:)))*1e7;  %penalty factor

for i=1:n_bcsdisp
  node_id=bcsdisp(i,1);
  direction=bcsdisp(i,2);
  K(3*node_id-3 + direction, 3*node_id-3 + direction)=coeff;
  F(3*node_id-3 + direction, 1) = bcsdisp(i,3)*coeff;
end
disp('Setting BCs complete');



%--------------------------------------------------
% solve the global linear system
%--------------------------------------------------  
disp('Start solving linear system');
u=K\F;  
disp('Solution complete');

%--------------------------------------------------
% save the displacement results in file
%--------------------------------------------------
U=zeros(n_nodes,7);

for n=1:n_nodes
    U(n,1)=n;
    U(n,2)=nodes(n,2);
    U(n,3)=nodes(n,3);
    U(n,4)=nodes(n,4);
    U(n,5)=u(3*n-2,1);
    U(n,6)=u(3*n-1,1);
    U(n,7)=u(3*n,1);
end

save -ascii -double feU.dat U

disp('Displacement results stored in feU.dat');

figure(1)
plot_mesh;

figure(2)
plotdeformedmesh;






