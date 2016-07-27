function combinedTopologyOptimization(useInputArgs, w1text, iterationNum)
% input args, useInputArgs = 1 if to use the input args
% w1 weight1 for weighted objective.
% iterationNum, used to know where to output files.

% --------------------------------------
% %% Settings
% --------------------------------------------
settings = Configuration;
settings.elasticMaterialInterpMethod = 2; % Hashin–Shtrikam law (average of upper and lower boundary)
settings.heatMaterialInterpMethod = 5; % Hashin–Shtrikam law (average of upper and lower boundary)

% target volumes of material 1 and 2
settings.v1 = 0.2; 
settings.v2 = 0.2;

settings.mode = 4; % 1 = topology only, 2 = material optimization only. 3 = both, 4 = meso-structure testing     

% if using input args, then override some configurations.
% if using input args, then running on the cluster, so use high resolution,
% otherwise use low resolution
if(str2num(useInputArgs) ==1)
    settings.w1 = str2num(w1text);
    
    settings.iterationNum = str2num(iterationNum);
    settings.nelx = 80;
    settings.nely = 40;
    
    settings.doPlotVolFractionDesignVar = 0;
    settings.doPlotTopologyDesignVar = 0;
    settings.doPlotHeat = 0;
    settings.doPlotHeatSensitivityTopology = 0;
    settings.doPlotStress = 0;
    settings.doPlotFinal = 0;
    settings.doSaveDesignVarsToCSVFile = 1; % set to 1 to write to csv file instead
    settings.maxFEACalls = 350;
    settings.maxMasterLoops = 500; % make it so, the fea maxes out first. 
else
    
    settings.nelx = 40;
    settings.nely = 20;
    settings.w1 = 0; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
    settings.iterationNum = 0;
    settings.doSaveDesignVarsToCSVFile = 1;
    settings.terminationCriteria =0.1; % 10%
    
    % if meso structure designing, then make a smaller initial mesh
    if(settings.mode == 4)
       settings.nelx = 20;
       settings.nely = 20;
        
    end
    
end


settings= settings.UpdateVolTargetsAndObjectiveWeights();
settings
% material properties Object
matProp = MaterialProperties;

% ---------------------------------
% Initialization of varriables
% ---------------------------------
designVars = DesignVars(settings);
designVars.x(1:settings.nely,1:settings.nelx) = settings.totalVolume; % artificial density of the elements
designVars.w(1:settings.nely,1:settings.nelx)  = 1; % actual volume fraction composition of each element
fractionCurrent_V1Local =1;

designVars.temp1(1:settings.nely,1:settings.nelx) = 0;
designVars.temp2(1:settings.nely,1:settings.nelx) = 0;
designVars.complianceSensitivity(1:settings.nely,1:settings.nelx) = 0;
designVars.totalStress(1:settings.nely,1:settings.nelx) = 0;

designVars.g1elastic(1:settings.nely,1:settings.nelx) = 0;
designVars.g1heat(1:settings.nely,1:settings.nelx) = 0;

designVars = designVars.CalcIENmatrix(settings);
designVars =  designVars.CalcElementLocation(settings);
designVars = designVars.PreCalculateXYmapToNodeNumber(settings);

% if doing meso optimization, setup optimization configurations
if ( settings.mode == 4) 
%     homgSettings = Configuration;
%     homgSettings.nelx = 40;
%     homgSettings.nely = 40;
%     homgSettings.w1 = 0; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
%     homgSettings.iterationNum = 0;
%     homgSettings.doSaveDesignVarsToCSVFile = 1;
%     homgSettings.terminationCriteria =0.1; % 10%
    
    designVars = designVars.CalcElementNodeMapmatrixWithPeriodicXandY(settings);
    designVars =  designVars.CalcNodeLocationMeso(settings);
end


masterloop = 0;
FEACalls = 0;
status=0;


% START ITERATION
while status == 0  && masterloop<=settings.maxMasterLoops && FEACalls<=settings.maxFEACalls
    masterloop = masterloop + 1;
    
    % --------------------------------
    % Topology Optimization
    % --------------------------------
    if ( settings.mode == 1 || settings.mode == 3)
        for loopTop = 1:3
            designVars = designVars.CalculateSensitivies(settings, matProp, masterloop);
            [vol1Fraction, vol2Fraction] =  designVars.CalculateVolumeFractions(settings);
            
            FEACalls = FEACalls+1;
            % normalize the sensitivies  by dividing by their max values.
            temp1Max =-1* min(min(designVars.temp1));
            designVars.temp1 = designVars.temp1/temp1Max;
            temp2Max = -1* min(min(designVars.temp2));
            designVars.temp2 = designVars.temp2/temp2Max;
            
            designVars.dc = settings.w1*designVars.temp1+settings.w2*designVars.temp2; % add the two sensitivies together using their weights
            
            % FILTERING OF SENSITIVITIES
            [designVars.dc]   = check(settings.nelx,settings.nely,settings.rmin,designVars.x,designVars.dc);
            % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
            [designVars.x]    = OC(settings.nelx,settings.nely,designVars.x,settings.totalVolume,designVars.dc, designVars, settings);
            % PRINT RESULTS
            %change = max(max(abs(designVars.x-designVars.xold)));
            disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',designVars.c) ...
                ' Vol. 1: ' sprintf('%6.3f', vol1Fraction) ...
                ' Vol. 2: ' sprintf('%6.3f', vol2Fraction) ...
                ' Lambda.: ' sprintf('%6.3f',designVars.lambda1  )])
            
            densitySum = sum(sum(designVars.x));
            designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            
            p = plotResults;
            p.plotTopAndFraction(designVars,  settings, matProp, FEACalls); % plot the results.
            status = TestForTermaination(designVars, settings);
            if(status ==1)
                m = 'break in topology'
                  break;
             end
        end
    end
    
    % exit the master loop if we termination criteria are true. 
    if(status ==1)
            m = 'exiting master (break in topology)'
        break;
    end

    % --------------------------------
    % Volume fraction optimization
    % --------------------------------
    if ( settings.mode ==2 || settings.mode ==3)
        for loopVolFrac = 1:3
            designVars = designVars.CalculateSensitivies( settings, matProp, masterloop);
            FEACalls = FEACalls+1;
            
            % for j = 1:5
            [vol1Fraction, vol2Fraction] =  designVars.CalculateVolumeFractions(settings);
            
            totalVolLocal = vol1Fraction+ vol2Fraction;
            fractionCurrent_V1Local = vol1Fraction/totalVolLocal;
            targetFraction_v1 = settings.v1/(settings.v1+settings.v2);
            
            % Normalize the sensitives.
            temp1Max = max(max(abs(designVars.g1elastic)));
            designVars.g1elastic = designVars.g1elastic/temp1Max;
            temp2Max = max(max(abs(designVars.g1heat)));
            designVars.g1heat = designVars.g1heat/temp2Max;
            
            g1 = settings.w1*designVars.g1elastic+settings.w2*designVars.g1heat; % Calculate the weighted volume fraction change sensitivity.
            G1 = g1 - designVars.lambda1 +1/(designVars.mu1)*( targetFraction_v1-fractionCurrent_V1Local); % add in the lagrangian
            designVars.w = designVars.w+settings.timestep*G1; % update the volume fraction.
            
            designVars.w = max(min( designVars.w,1),0);    % Don't allow the    vol fraction to go above 1 or below 0
            designVars.lambda1 =  designVars.lambda1 -1/(designVars.mu1)*(targetFraction_v1-fractionCurrent_V1Local)*settings.volFractionDamping;
            
            %  designVars.w = OC_gradient(g1,designVars, settings)  ;
            
            
            
            
            % PRINT RESULTS
            %change = max(max(abs(designVars.x-designVars.xold)));
            designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            p = plotResults;
            p.plotTopAndFraction(designVars, settings, matProp,FEACalls ); % plot the results.
         
            
            disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',designVars.c) ...
                ' Vol. 1: ' sprintf('%6.3f', vol1Fraction) ...
                ' Vol. 2: ' sprintf(    '%6.3f', vol2Fraction) ...
                ' Lambda.: ' sprintf('%6.3f',designVars.lambda1  )])            
            % obj.storeOptimizationVar = [obj.storeOptimizationVar;obj.c,  obj.cCompliance, obj.cHeat ];    
              status = TestForTermaination(designVars, settings);
              if(status ==1)
                    m = 'break in vol fraction'
                  break;
              end
        end
    end
    
    if(settings.mode ==4) % meso-structure design
         settings.maxMasterLoops = 3; % make it small
         
         % this vector is just an element diplacement vector pulled from
         % the main displacement vector U
         
       
          % need to convert the strain into a displacement B^-1
          % epsilon = B*d
        material1Fraction = 1;
        macroElemProps = macroElementProp;          
        macroElemProps.disp = [-0.1 0    0.1 0    0.1 0   -0.1 0] ; % make these up for now
        [macroElemProps.K ,~,macroElemProps.B] = matProp.effectiveElasticKEmatrix(material1Fraction, settings);
        macroElemProps.strain = macroElemProps.B* transpose(macroElemProps.disp); % transpose the disp to be vertical
        
        
        p = plotResults;
          for mesoLoop = 1:4
              designVars = designVars.CalculateSensitiviesMesoStructure( settings, matProp, masterloop,macroElemProps);
              designVars.dc=designVars.temp1;              
 
             % FILTERING OF SENSITIVITIES
             [designVars.dc]   = check(settings.nelx,settings.nely,settings.rmin,designVars.x,designVars.dc);
             % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD          
             [designVars.x]    = OC(settings.nelx,settings.nely,designVars.x,settings.totalVolume,designVars.dc, designVars, settings);
             
             p.plotTopAndFraction(designVars, settings, matProp,FEACalls ); % plot the results.
          end
            macroElemProps = designVars.GetHomogenizedProperties(settings,settings, matProp, masterloop,macroElemProps);
            macroElemProps.D_homog
            
         
           %  p.plotTopAndFraction(mesoDesignVars, homgSettings, matProp,FEACalls ); % plot the results.
        
      
        
    end
end

folderNum = settings.iterationNum;
outname = sprintf('./out%i/storeOptimizationVar.csv',folderNum);
csvwrite(outname,designVars.storeOptimizationVar);
status


% test for termaination of the function. 
% normalize the objectives, compare to a moving average. If moving average
% below target, then return 1
% if not terminated, then return 0
function status = TestForTermaination(designVars, settings)
    status = 0;
    y2 = designVars.storeOptimizationVar(:,2); % Elastic Compliance   
    
    t = settings.terminationAverageCount;
    if(size(y2)<(t+3))
        return;
    end
  
  
    y3 = designVars.storeOptimizationVar(:,3); % Heat Compliance
    y4 = designVars.storeOptimizationVar(:,4); % material 1
    y5 = designVars.storeOptimizationVar(:,5); % material 2
    
    avg_y2 = FindAvergeChangeOfLastValue(y2,settings); % elastic objective
    avg_y3 = FindAvergeChangeOfLastValue(y3,settings); % heat objective
    avg_y4 = FindAvergeChangeOfLastValue(y4,settings); % material 1
    avg_y5 = FindAvergeChangeOfLastValue(y5,settings); % material 2
    
    
   
    
    tt = settings.terminationCriteria;
    if(        avg_y2<tt ...
            && avg_y3<tt ...
            && avg_y4<tt ...
            && avg_y5<tt)
        status=1;
        % print the vars to screen
     [avg_y2 avg_y3 avg_y4 avg_y5 ] 
    end
    
    
function [averageDiffOfLastValues] = FindAvergeChangeOfLastValue(arrayOfValues, settings)
y2 = arrayOfValues;
t = settings.terminationAverageCount;
diff2 = y2-circshift(y2,1); % subtract the current one from the pervious on 
diff2 = diff2(2:end-1); 

 % normalize and take ABS
 % we want to normalize by max(abs) over all differences, so 
diff2 = abs(diff2)/(max(abs(diff2)));

diff2 = diff2(end-t:end);
averagey2diff = sum(diff2)/t;
averageDiffOfLastValues = averagey2diff;

% if recvid==1
%          close(vidObj);  %close video

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc ,designVar, settings)

% if the sensitivity is really small, then make it larger to help the
% optimal criteria method work better.
absL =  max(max(abs(dc)));
if absL <1000
    dc = dc*100000.0/absL;
end

l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    xnew = max(0.01,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
    
    %   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
    
    %[volume1, volume2] = designVar.CalculateVolumeFractions(settings);
    %currentvolume=volume1+volume2;
    
    %if currentvolume - volfrac > 0;
    
    if sum(sum(xnew)) - volfrac*nelx*nely > 0;
        l1 = lmid;
    else
        l2 = lmid;
    end
end


%%%%%%%%%% OPTIMALITY CRITERIA UPDATE Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wnew]=OC_gradient(g1,designVar, settings)
l1 = 0; l2 = 100000; move = 0.05;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    designVar.w = max(0.001,max(designVar.w-move,min(1.,min(designVar.w+move,designVar.w.*sqrt(-g1./lmid)))));
    
    %   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
    
    [volume1, volume2] = designVar.CalculateVolumeFractions(settings);
    currentTotal = volume1+volume2;
    %currentvolume=volume1+volume2;
    
    %if currentvolume - volfrac > 0;
    if currentTotal - settings.totalVolume > 0;
        l1 = lmid;
    else
        l2 = lmid;
    end
    
    wnew =  designVar.w ;
end

%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
    for j = 1:nely
        sum=0.0;
        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
            for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                fac = rmin-sqrt((i-k)^2+(j-l)^2);
                sum = sum+max(0,fac);
                dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
            end
        end
        dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
    end
end
