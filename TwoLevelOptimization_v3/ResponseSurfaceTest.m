function [] = ResponseSurfaceTest()

 matProp = MaterialProperties;

usetheta =0;
if(usetheta==0)
    
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
    
    ResponseExx=ResponseExx/ matProp.E_material1;
    ResponseEyy=ResponseEyy/ matProp.E_material1;
    
    sf = fit([ResponseExx, ResponseEyy],ResponseDensity,'poly22')
    plot(sf,[ResponseExx,ResponseEyy],ResponseDensity)
    xlabel('ResponseExx');
    ylabel('ResponseEyy');
    zlabel('ResponseDensity');
  coeffvalues(sf)
    
end
if(usetheta==1)
    % Exx and Eyy is normalizes so taht the max is 1
    
    % Symmetric 
%     Exx is always greater than Eyy
%     theta between 0 and pi/4
    
      % point 1, Exx extreme
    point=1;
    ResponseExx(1) = 1;
    ResponseEyy(1) = 0;
    ResponseTheta(1) = pi/4;
    ResponseDensity(1) = 0.8;
    
    % point 2
    point = point+1;
    ResponseExx(point) = 1;
    ResponseEyy(point) = 0;
    ResponseTheta(point) = 0;
    ResponseDensity(point) = 0.7;
    
    % point 3
    point = point+1;
    ResponseExx(point) = 1;
    ResponseEyy(point) = 1;
    ResponseTheta(point) = 0;
    ResponseDensity(point) = 1;
    
      % point 4
    point = point+1;
    ResponseExx(point) = 1;
    ResponseEyy(point) = 1;
    ResponseTheta(point) = pi/4;
    ResponseDensity(point) = 1;
    
       % point 5
    point = point+1;
    ResponseExx(point) = 0.5;
    ResponseEyy(point) = 0.5;
    ResponseTheta(point) = pi/4;
    ResponseDensity(point) = 0.5;
    
      % point 6
    point = point+1;
    ResponseExx(point) = 0.5;
    ResponseEyy(point) = 0;
    ResponseTheta(point) = pi/4;
    ResponseDensity(point) = 0.2;
    
        % point 6
    point = point+1;
    ResponseExx(point) = 0.5;
    ResponseEyy(point) = 0;
    ResponseTheta(point) = 0;
    ResponseDensity(point) = 0.15;
    
    
  
    
       x0=ones(1,10);
%        x0 = randi([-5,5],1,10);
        A = [];
        b = [];
               
       
        % theta the same
        % rho the same.
        % scale down Exx, Eyy
       
        X = ResponseExx;
        Y = ResponseEyy;
        
%         Z = thetaArray/(pi/4);
          Z = ResponseTheta;
         R = ResponseDensity;
        
        
        ub = ones(6,1)*10000;
        lb = -ub;
         o=Optimizer;
         config = Configuration;
         config= config.UpdateRunTimeComputedSettings( '1', '1', '2','60', 'na');
         
        [coefficients finalObjective]= fmincon(@(x) fitObjectiveV2(x,X,Y,Z,R,o,config,matProp),x0,A,b,[],[],lb,ub);
       
        finalObjective
        
        coefficients
        
        
        % -------------------------
        % plot the results. 
        % -------------------------
         tt  =0.1;
    
        [Xgrid, Ygrid, Zgrid]=meshgrid(0:tt:1,0:tt:1,0:pi/20:pi/4);
     
        % Reshape into columns
        E_xx=reshape(Xgrid,[],1);
        E_yy=reshape(Ygrid,[],1);
        theta=reshape(Zgrid,[],1);
        
            
        circleSize = ones(size(E_xx))*100; % circle size.
              [~, ~,rhoExperimental] = o.CalculateDensitySensitivityandRho(E_xx,E_yy,theta,coefficients,config,matProp);
        scatter3(E_xx,E_yy,theta,circleSize,rhoExperimental,'filled');
        title(sprintf('Response surface,Rho (the color) as a function of Exx, Eyy, theta: iter %i',config.macro_meso_iteration));
        colorbar
         xlabel('Exx');
        ylabel('Eyy');
        zlabel('Theta');
        
        hold on
        
      RhoColor=ResponseDensity;
     circleSize= ones(size(ResponseExx))*1000; % circle size.
          scatter3(ResponseExx,ResponseEyy,ResponseTheta,circleSize,RhoColor,'filled');
          hold off
   
     
        
    
    
    
    
    
    
    
    
    
    
    
    
    
end