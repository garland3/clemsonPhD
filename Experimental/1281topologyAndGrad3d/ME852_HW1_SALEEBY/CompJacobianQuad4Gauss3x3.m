
%--------------------------------------------------------------------
% Compute Jacobian matrix of 4 node quadrilateral element
% with 2x2 Gauss points (i.e. 1,2,3,4 points shown below)
%
%       -------
%      | 2   1 |
%      |       |
%      | 3   4 |
%       -------
% 
% Input:
%   element_nodes: the nodal coordinates of the element
%   g: g-th Gauss point for which the Jacobian matrix is computed
% Output:
%   J: Jacobian matrix at the g-th Gauss point in the given element
%--------------------------------------------------------------------

function J= CompJacobianQuad4Gauss3x3(element_nodes, g, Nx, Ny, Nz)
  J=zeros(3,3);
  for j=1:3
    J(1,j) =  (Nx(:,g))' * element_nodes(:,j);
    J(2,j) =  (Ny(:,g))' * element_nodes(:,j);
    J(3,j) =  (Nz(:,g))' * element_nodes(:,j);
  end    
