
function F=CompF(nodes, elements, materials, bcstraction, bcsforce)
  
  n_nodes = size(nodes,1);
  n_elements = size(elements,1);
  n_traction_edges = size(bcstraction,1);
  n_force_nodes = size(bcsforce,1);
  n_nodes_per_element = size(elements,2)-1;
  n_edge_gauss_points=4;
  traction=zeros(3,1);
  
  F=zeros(n_nodes*3,1);
  fe=zeros(n_nodes_per_element*3,1);
  
  [N,Nx,Ny,Nz]=CompNDNQuad4EdgeGauss8x1();
  
  
  %------------------------------------------------------
  % Apply traction force
  %------------------------------------------------------
  
  %---if bcstraction(1,1)==0 then there is no traction force
  if bcstraction(1,1)==0
    n_traction_edges=0;
  end
  
  %---Loop over the number of edges affected by traction
  for t=1:n_traction_edges
    fe=zeros(n_nodes_per_element*3, 1);
    element_id=bcstraction(t,1);
    for i=1:n_nodes_per_element
      global_node_id=elements(element_id,i+1);
      element_nodes(i,1) = nodes(global_node_id,2);
      element_nodes(i,2) = nodes(global_node_id,3);
      element_nodes(i,3) = nodes(global_node_id,4);
      node_id_map(i,1) = global_node_id;
    end
    edge=bcstraction(t,2);
    traction(1,1)=bcstraction(t,5);
    traction(2,1)=bcstraction(t,6);
    traction(3,1)=bcstraction(t,7);

    for g=1:3
      gauss_pt_id=3*edge-3+g;
      J=CompEdgeJacobianQuad4Gauss8x1(element_nodes,gauss_pt_id,Nx,Ny,Nz);
      if  (edge==1);  
        lengthJ=sqrt(J(1,1)^2+J(1,2)^2)+J(1,3)^2;
      elseif (edge==2);
        lengthJ=sqrt(J(2,1)^2+J(2,2)^2)+J(2,3)^2;
      else
          lengthJ=sqrt(3,1)^2+J(3,2)^2+J(3,3)^2;
      end
      
      fxy=zeros(24,3);
      for k=1:n_nodes_per_element   
        fxy(3*k-2,1)=N(k, gauss_pt_id);
        fxy(3*k-1,2)=N(k, gauss_pt_id);
        fxy(3*k,3)=N(k, gauss_pt_id);
      end
      fe =fe+ fxy*traction*lengthJ;
    end

    %---Assemble F
    for j=1:n_nodes_per_element;
      row=node_id_map(j,1)*3-2;
      F(row:row+2,1) = F(row:row+2,1)+ fe(3*j-2:3*j, 1); 
    end
  end

  %---------------------------------------------------------
  % Apply point forces
  %---------------------------------------------------------
  
  %---if bcstraction(1,1)==0 then there is no point force
  if bcsforce(1,1)==0
    n_force_nodes=0;
  end
  
  for i=1:n_force_nodes
    row=3*bcsforce(i,1)-2;
    dir=bcsforce(i,2);
    F(row+dir-1,1) = F(row+dir-1,1)+bcsforce(i,3);%/materials(3);
%     F(row,1) = F(row,1)+ bcsforce(i,2)/materials(3);
%     F(row+1,1) = F(row+1,1)+ bcsforce(i,3)/materials(3);
    %---Note: materials(3) is thickness
  end
  
  