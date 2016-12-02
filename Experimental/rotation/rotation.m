theta = 45;
theta_r = deg2rad(theta);

m = cos(theta_r);
n = sin(theta_r);

T = [m^2 n^2 2*m*n; ...
    n^2 m^2 -2*m*n; ...
    -m*n m*n (m^2-n^2)];


v= 0.3;
E = 1000;
 D = [ 1 v 0;
    v 1 0;
    0 0 1/2*(1-v)]*E/(1-v^2)

T_inv = inv(T);

D_new = T_inv*D*inv(T_inv)

    
