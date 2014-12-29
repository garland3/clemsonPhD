function k = elementK_heat()

if(1==1)    
    % use the precalculated k
    % or if necessary we can actually do the calculations. 
    k = [0.6667   -0.1667   -0.3333   -0.1667;
           -0.1667    0.6667   -0.1667   -0.3333;
           -0.3333   -0.1667    0.6667   -0.1667;
           -0.1667   -0.3333   -0.1667    0.6667];
else
    kmaterial = 1; % W/degreeC
    % 1 by 1 square
    coord = [0 0;
            1 0;
            1 1;
            0 1];

    % Loop over the Guassian quadrature points. We will use 2 point
    % quadrature for a total of 4 points
    etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
    zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
    weight = [ 1 1 1 1];

    ke = zeros(4,4);
    btemp = zeros(2,4);

    for gu = 1:4
      eta = etaRow(gu);
      Zeta = zetaRow(gu);
      wght = weight(gu);

      % Calculate the Shape functions.
      N(1) = 1/4*(1-eta)*(1-Zeta);
      N(2) = 1/4*(1+eta)*(1-Zeta);
      N(3) = 1/4*(1+eta)*(1+Zeta);
      N(4) = 1/4*(1-eta)*(1+Zeta);

      % B_hat (Derivative of N1 with respect to zeta and eta)
      B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
                   -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];

      % Calculate the Jacobian
      J=B_hat*coord;

      % Calculate the determinate
      J_det = det(J);
    
       B = J\B_hat;
      btemp = btemp +B;

      tempK = transpose(B)*kmaterial*B*J_det*wght;          
       ke = ke + tempK;% +tempM;
    end  
    k = ke;
end