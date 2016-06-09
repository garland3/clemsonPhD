mat = MaterialProperties;

% ceramic Al2O3
%mat.E_material2 = 380e9; % The elastic mod of material 1, 
%mat.K_material2 = 8.4; % heat conduction of material 1

mat.E_material2 =0.25; % The elastic mod of material 1, 
mat.K_material2 =1; % heat conduction of material 1


% TiC
%mat.E_material1 = 450e9; % The elastic mod of material 2
% mat.K_material1 = 24.28; % heat conduction of material 2    

mat.E_material1 = 1; % The elastic mod of material 2
mat.K_material1 = 5; % heat conduction of material 2    
   
   
settings = Configuration;
settings.elasticMaterialInterpMethod = 2; % Hashin–Shtrikam law (average of upper and lower boundary)
settings.heatMaterialInterpMethod = 5; % Hashin–Shtrikam law (average of upper and lower boundary)




volfractionsM1 = 0:0.01:1;

count = 1;
for i = volfractionsM1
    eff(count) =  mat.effectiveElasticProperties( i, settings);
    effH(count) = mat.effectiveHeatProperties(i,settings);
    count = count+1;
end

subplot(2,1,1)
plot(volfractionsM1,eff)

subplot(2,1,2)
plot(volfractionsM1,effH)

