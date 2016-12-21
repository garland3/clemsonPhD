function [ke]=elK_elastic_v2(Dgiven)
% Calculates the K stiffness matrix, when given the constiutive matrix
% if the strain matrix is not empty, then calculate the Force required to
% strain the element the strain amout.


% E = materialProp.E_material1;
% v = materialProp.v; % Piossons ratio
% G = materialProp.G;

% Plane stress problem.
%         D = [ 1 v 0;
%             v 1 0;
%             0 0 1/2*(1-v)]*E/(1-v^2);

D=Dgiven;

% 1 by 1 square of coordinates
coord = [0 0;
        1 0;
        1 1;
        0 1];

% etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
etaRow(1,:)= [0.5774    0.5774   -0.5774   -0.5774];
% zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
zetaRow(1,:) =  [  0.5774   -0.5774   -0.5774    0.5774];
weight = [ 1 1 1 1];

ke = zeros(8,8);
% F_meso = zeros(8,1);
% kExpansion_bar = 0;
% B_total = zeros(3,8);

% Loop over the guass points
for gu = 1:4
    eta = etaRow(gu);
    Zeta = zetaRow(gu);
    wght = weight(gu);    
    %           % Calculate the Shape functions.
    %           N(1) = 1/4*(1-eta)*(1-Zeta);
    %           N(2) = 1/4*(1+eta)*(1-Zeta);
    %           N(3) = 1/4*(1+eta)*(1+Zeta);
    %           N(4) = 1/4*(1-eta)*(1+Zeta);
    
    % B_hat (Derivative of N1 with respect to zeta and eta)
    B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
        -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];
    
    % Calculate the Jacobian
    J=B_hat*coord;
    
    % Calculate the determinate
    J_det = det(J);
    % J_inverse = inv(J);
    
    % Form the B matrix
    % B_2by4 = J_inverse*B_hat;
    B_2by4_v2 = J\B_hat;
    
    % Form B, which is an 3 by 8
    B = zeros(3,8);
    B(1,[1,3,5,7]) = B_2by4_v2(1,1:4);
    B(2,[2,4,6,8]) = B_2by4_v2(2,1:4);
    
    B(3,[1,3,5,7]) = B_2by4_v2(2,1:4);
    B(3,[2,4,6,8]) = B_2by4_v2(1,1:4);
%     B_total = B_total+ B*J_det*wght;    
    
    tempK = transpose(B)*D*B*J_det*wght;
    ke = ke + tempK;
    
    % Calcualte the thermal expansion matrix
%     f_bar_temp =  transpose(B)*D*transpose([1 1  0])*J_det*wght;
%     kExpansion_bar = kExpansion_bar+f_bar_temp;
end


% calculate the load force for a strain
% if(~isempty(strain))
%     F_meso = transpose(B_total)*D*transpose(strain);
% end







end
