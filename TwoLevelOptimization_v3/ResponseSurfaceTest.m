matProp = MaterialProperties;

% point 1, Exx extreme
point=1;
ResponseExx(1) = matProp.E_material1;
ResponseEyy(1) = 0;
ResponseDensity(1) = 1;


point = point+1;
ResponseExx(point) = matProp.E_material1;
ResponseEyy(point) = matProp.E_material1/2;
ResponseDensity(point) = 1;

% both extreme
point = point+1;
ResponseExx(point) = matProp.E_material1;
ResponseEyy(point) = matProp.E_material1;
ResponseDensity(point) = 1;

% Eyy exteme
point = point+1;
ResponseExx(point) = 0;
ResponseEyy(point) = matProp.E_material1;
ResponseDensity(point) = 1;


point = point+1;
ResponseExx(point) = matProp.E_material1/2;
ResponseEyy(point) = matProp.E_material1;
ResponseDensity(point) = 1;

point

% smallest possible. 
point = point+1;
ResponseExx(point) = matProp.E_material2/2;
ResponseEyy(point) =  matProp.E_material2/2;
ResponseDensity(point) = 0.2; % this is just a guess

ResponseExx = transpose(ResponseExx);
ResponseEyy = transpose(ResponseEyy);
ResponseDensity = transpose(ResponseDensity);

sf = fit([ResponseExx, ResponseEyy],ResponseDensity,'poly22')
plot(sf,[ResponseExx,ResponseEyy],ResponseDensity)
 xlabel('ResponseExx');
ylabel('ResponseEyy');
zlabel('ResponseDensity');