
function K=CompK(nodes, elements, materials)
  
  n_nodes = size(nodes,1);
  n_elements = size(elements,1);
  n_nodes_per_element = size(elements,2)-1;
  n_gauss_points=8;
  
  
  K=sparse(n_nodes*3,n_nodes*3);  % create empty sparse matrix K
  
  element_nodes=zeros(n_nodes_per_element,3); 
  node_id_map=zeros(n_nodes_per_element,1);
  C=CompCPlaneStress(materials);       
  H=CompH();
  [N,Nx,Ny,Nz]=CompNDNQuad4Gauss3x3();
  [gauss_points, gauss_weights]=GetQuadGauss3x3();
    
  %----------------------------------------------------------
  % Compute K matrix: loop over all the elements 
  %----------------------------------------------------------
  for e=1:n_elements
    ke=zeros(n_nodes_per_element*3, n_nodes_per_element*3);
    for i=1:n_nodes_per_element
      global_node_id=elements(e,i+1);
      element_nodes(i,1) = nodes(global_node_id,2);
      element_nodes(i,2) = nodes(global_node_id,3);
      element_nodes(i,3) = nodes(global_node_id,4);
      node_id_map(i,1) = global_node_id;
    end
    
    %---------
    % compute element stiffness matrix ke
    %---------     
    for g=1:n_gauss_points
      J=CompJacobianQuad4Gauss3x3(element_nodes,g, Nx, Ny, Nz);
      detJ=det(J);
      %if detJ<0.0001  return; end;
      Jinv=inv(J);
      Jb(1:3,1:3)=Jinv;
      Jb(4:6,4:6)=Jinv;
      Jb(7:9,7:9)=Jinv;
      B=CompB4x8Quad4Gauss2x2(g, Nx, Ny, Nz);
      HJB=(H*Jb*B);
      ke=ke+HJB'*C*HJB*detJ*gauss_weights(g);
    end
    
    %---------
    % assemble ke into global K
    %---------
    for j = 1:n_nodes_per_element    
      row_node= node_id_map(j,1);
      row=3*row_node-2;
      for k = 1:n_nodes_per_element  
        col_node= node_id_map(k,1);
        col=3*col_node-2;

        K(row, col)= K(row, col) + ke(3*j-2, 3*k-2);
        K(row+1, col) = K(row+1, col) + ke(3*j-1, 3*k-2);
        K(row+2, col) = K(row+2, col) + ke(3*j, 3*k-2);
        
        K(row, col+1)= K(row, col+1) + ke(3*j-2, 3*k-1);
        K(row+1, col+1) = K(row+1, col+1) + ke(3*j-1, 3*k-1);
        K(row+2, col+1) = K(row+2, col+1) + ke(3*j, 3*k-1);
        
        K(row, col+2)= K(row, col+2) + ke(3*j-2, 3*k);
        K(row+1, col+2) = K(row+1, col+2) + ke(3*j-1, 3*k);
        K(row+2, col+2) = K(row+2, col+2) + ke(3*j, 3*k);
      end
    end
    
  end
  