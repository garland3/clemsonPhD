function  macroEleProps = GetMacroElementPropertiesFromCSV(config,e)
% for the new meso design method using the Exx and Eyy and the consistency
% constraints I need
% 1. Topology var. Tells me if I need to make a new design or not. DONE
% 2. XY position map, DONE
% 3. D_system matrix for each element in the macro design.
% 4. Displacement field. So that I can calculate the strain on each. DONE
% element.


if(config.validationModeOn==0 && config.strainAndTargetTest==0)
    % ----------------------------------------------
    %
    %      Normal Multiscale optimization
    %
    % ----------------------------------------------
    
    mm_iteration = config.macro_meso_iteration;
    macroEleProps = macroElementProp;
    macroEleProps.elementNumber = e;
    
    folderNum = config.iterationNum;
    
    % Get element->node mapping
    outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,mm_iteration);
    IEN = csvread(outname);
    
    % % Get displacement fields
    outname = sprintf('./out%i/displacement%i.csv',folderNum,mm_iteration);
    U =  csvread(outname);
    
    % -----------------------------------
    %
    % 1 element per design var case on macro level
    %
    % -----------------------------------
    % GET the saved element to XY position map (needed for x and w vars retrival)
    outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,mm_iteration);
    elementXYposition=csvread(outname);
    results = elementXYposition(macroEleProps.elementNumber,:);
    macroEleProps.yPos = results(1);
    macroEleProps.xPos = results(2);
    
    % Get the density field
    outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
    x = csvread(outname);
    macroEleProps.densitySIMP = x(macroEleProps.yPos,macroEleProps.xPos );
    
    % Get the Exx field
    outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
    ExxMacro = csvread(outname);
    macroEleProps.Exx = ExxMacro(macroEleProps.yPos,macroEleProps.xPos );
    
    % Get the Eyy field
    outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
    EyyMacro =csvread(outname);
    macroEleProps.Eyy = EyyMacro(macroEleProps.yPos,macroEleProps.xPos );
    
    % Get the Theta field
    outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
    ThetaMacro = csvread(outname);
    macroEleProps.theta = ThetaMacro(macroEleProps.yPos,macroEleProps.xPos );
    
    w=1;
    matProp=MaterialProperties;
    D= matProp.getDmatMatrixTopExxYyyRotVars(config,macroEleProps.densitySIMP ,macroEleProps.Exx, macroEleProps.Eyy,macroEleProps.theta ,w);
    
    % outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    % D = csvread(outname);
    macroEleProps.D_sys =D;
    
    if(macroEleProps.densitySIMP>config.noNewMesoDesignDensityCutOff ||config.multiscaleMethodCompare==1)
        nodes1=  IEN(e,:);
        macroEleProps.elementNodes=nodes1;
        xNodes = nodes1*2-1;
        yNodes = nodes1*2;
        dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
        
        % plan for multi-loading cases.
        [~, t2] = size(config.loadingCase);
        for loadcaseIndex = 1:t2
            utemp = U(loadcaseIndex,:);
            u_local =   utemp(dofNumbers);
            
            %         offsetX = mean(u_local([1 3 5 7]));
            %         offsetY = mean(u_local([2 4 6 8]));
            offsetX = u_local(7);
            offsetY = u_local(8);
            u_local([1 3 5 7]) = u_local([1 3 5 7])-offsetX;
            u_local([2 4 6 8]) = u_local([2 4 6 8])-offsetY;
            macroEleProps.disp(loadcaseIndex,:)  = u_local;
        end
        
        % % Get the volume fraction field (but only if using the old method)
        if(config.useExxEyy~=1)
            outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,mm_iteration);
            w = csvread(outname);
            macroEleProps.material1Fraction = w(macroEleProps.yPos,macroEleProps.xPos );
        else
            macroEleProps.material1Fraction = 1; % set to 1 as a place holder
        end
    else
        fprintf('SIMP density is below threshold\n');
    end
elseif(config.validationModeOn==1)
    % -----------------------------------
    %
    %       Validation problem for meso designs.
    %
    % -----------------------------------
    folderNum=0;
    mm_iteration=1;
    % read the Exx field
    outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
    ExxSaved= csvread(outname);
    
    % read the Eyy field
    outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
    EyySaved= csvread(outname);
    
    % read the Theta field
    outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
    ThetaValues=csvread(outname);
    
    Exx=ExxSaved(e);
    Eyy=EyySaved(e);
    theta=ThetaValues(e);
    fprintf('Meso design Validation mode. Targets Exx %f Eyy %f Theta %f\n',Exx,Eyy,theta);
    
    macroEleProps = macroElementProp;
    macroEleProps.elementNumber = e;
    macroEleProps.yPos =1;
    macroEleProps.xPos = e;
    macroEleProps.material1Fraction=1;
    macroEleProps.disp=ones(1,8);
    
    
    macroEleProps.densitySIMP = 1;
    macroEleProps.Exx =Exx;
    macroEleProps.Eyy =Eyy;
    macroEleProps.theta=theta;
    
    
    w=1;
    matProp=MaterialProperties;
    D= matProp.getDmatMatrixTopExxYyyRotVars(config,macroEleProps.densitySIMP ,macroEleProps.Exx, macroEleProps.Eyy,macroEleProps.theta ,w);
    
    % outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    % D = csvread(outname);
    macroEleProps.D_sys =D;
    
    
    
elseif(config.strainAndTargetTest==1)
    % -----------------------------------
    %
    %      Psuedo Strain and target density test for ANN training data
    %
    % -----------------------------------
    folderNum=0;
    mm_iteration=1;
    % pseudoStrain
    outname = sprintf('./out%i/strain1%i.csv',folderNum,mm_iteration);
    strain1= csvread(outname);
    
    % pseudoStrain
    outname = sprintf('./out%i/strain2%i.csv',folderNum,mm_iteration);
    strain2= csvread(outname);
    
    % pseudoStrain
    outname = sprintf('./out%i/strain3%i.csv',folderNum,mm_iteration);
    strain3=csvread(outname);
    
    % targetDensity
    outname = sprintf('./out%i/densityTargets%i.csv',folderNum,mm_iteration);
    densityTarget=csvread(outname);
    
    %     Exx=strain1(e);
    %     Eyy=strain2(e);
    %     theta=ThetaValues(e);
    
    %  e=e+60000;
    
    macroEleProps = macroElementProp;
    macroEleProps.elementNumber = e;
    macroEleProps.yPos =1;
    macroEleProps.xPos = e;
    macroEleProps.material1Fraction=1;
    macroEleProps.disp=ones(1,8);
    
    macroEleProps.psuedoStrain=[1; 1; 1];
    macroEleProps.psuedoStrain(1) = strain1(e);
    macroEleProps.psuedoStrain(2) = strain2(e);
    macroEleProps.psuedoStrain(3) = strain3(e);
    temp1 =densityTarget(e);
    macroEleProps.targetDensity=temp1;
    
    
    macroEleProps.densitySIMP = 1;
    macroEleProps.Exx =1;
    macroEleProps.Eyy =1;
    macroEleProps.theta=1;
    
    if(    macroEleProps.psuedoStrain(2)> macroEleProps.psuedoStrain(1) )
        macroEleProps.targetDensity=-1;
        fprintf('Do not use in test, ps2>ps1 symmetric\n');
        return;
    end
    
    
    w=1;
    matProp=MaterialProperties;
    D= matProp.getDmatMatrixTopExxYyyRotVars(config,macroEleProps.densitySIMP ,macroEleProps.Exx, macroEleProps.Eyy,macroEleProps.theta ,w);
    
    % outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    % D = csvread(outname);
    macroEleProps.D_sys =D;
    
    %printf('Finished getting macro data from csv files');
end

if(config.UseLookUpTableForPsuedoStrain==1 && config.strainAndTargetTest~=1)
    
    % 1. REad in the data look up table.
    % 2. (Done) find the D_h for the coorospnding Exx, Eyy, Theta values
    % 3. Search the lookup table for the value that minimizes
    % Save this min value as the psuedo strain to use.
    macro_meso_iteration=1;
    if( config.mesoDesignInitalConditions==1)
        folder = '_random';
    elseif(config.mesoDesignInitalConditions==3)
        folder = '3_circle';
        folder='';
    elseif(config.mesoDesignInitalConditions==7)
        %              folder ='_circleSolid';
        folder='';
    end
    outname = sprintf('./data%s/D11%i_datat.data',folder,macro_meso_iteration);
    D11 =csvread(outname);
    
    outname = sprintf('./data%s/D12%i_datat.data',folder,macro_meso_iteration);
    D12 =csvread(outname);
    
    outname = sprintf('./data%s/D22%i_datat.data',folder,macro_meso_iteration);
    D22 =csvread(outname);
    
    outname = sprintf('./data%s/D33%i_datat.data',folder,macro_meso_iteration);
    D33 =csvread(outname);
    
    % read the targets
    outname = sprintf('./data%s/pstrain1%i_datat.data',folder,macro_meso_iteration);
    pstrain1 =csvread(outname);
    
    outname = sprintf('./data%s/pstrain2%i_datat.data',folder,macro_meso_iteration);
    pstrain2 =csvread(outname);
    
    outname = sprintf('./data%s/pstrain3%i_datat.data',folder,macro_meso_iteration);
    pstrain3 =csvread(outname);
    
    outname = sprintf('./data%s/etaTarget%i_datat.data',folder,macro_meso_iteration);
    etaTarget =csvread(outname);
    
    
    version=2;
    if(version==1)
%         %%         % ------------------------------------
%         %         %           Version 1
%         %         % ------------------------------------
%                 [~, t2]=size(etaTarget);
%                 minValue=1e11;
%                 indexOfMinValue = 1;
%                 flip = 1;
%                 D11sys = max(macroEleProps.D_sys(1,1),0.001 );
%                 D12sys = max(macroEleProps.D_sys(1,2),0.001 );
%                 D22sys =max( macroEleProps.D_sys(2,2) ,0.001 );
%                 D33sys = max(macroEleProps.D_sys(3,3),0.001 );
%                 for i = 1:t2
%                     D11_table = D11(i);
%                     D12_table = D12(i);
%                     D22_table = D22(i);
%                     D33_table = D33(i);
%                     etaLocal = etaTarget(i);
%                    if(etaLocal<-0.01)
%                         continue
%                     end
%                     for j = 1:2
%                         % try flipping the D11 and D22
%                         if(j==2)
%                             temp = D22_table;
%                             D22_table=D11_table;
%                             D11_table=temp;
%                         end
%                         diffValue = (D11_table-D11sys)^2+(D12_table-D12sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
%                         ps(1) = pstrain1(i);
%                         ps(2)= pstrain2(i);
%                         ps(3) = pstrain3(i);
%                         ps2 = sum(abs( ps));
%                         if(ps2<0.5 &&etaLocal<0.4 &&config.mesoDesignInitalConditions==3 && config.mesoDesignInitalConditions==3)
%                             continue
%                         end
%                         if(diffValue<minValue && etaLocal>config.MesoMinimumDensity )
%                             minValue=diffValue;
%                             indexOfMinValue=i;
%                             flip = j;
%                         end
%                     end
%                 end
%         
%                 D11_table = D11(indexOfMinValue);
%                 D12_table = D12(indexOfMinValue);
%                 D22_table = D22(indexOfMinValue);
%                 D33_table = D33(indexOfMinValue);
%                 originalEta = etaTarget(indexOfMinValue);
%                 if(flip==2)
%                     temp = D22_table;
%                     D22_table=D11_table;
%                     D11_table=temp;
%                 end
%                 matProp=MaterialProperties;
%                 smallestDiff = minValue;
%                 bestScale=0;
%                 for scale = -1:0.001:1
%                     D11_table2=D11_table+(scale)*D11_table;
%                     D12_table2 = D12_table+(scale)*D12_table ;
%                     D22_table2=D22_table+(scale)*D22_table;
%                     D33_table2=D33_table+(scale)*D33_table ;
%                     diffValue = (D11_table2-D11sys)^2+(D12_table2-D12sys)^2+(D22_table2-D22sys)^2+(D33_table2-D33sys)^2;
%                     if(diffValue<smallestDiff)
%                         bestScale=scale;
%                         smallestDiff=diffValue;
%                     end
%                 end
%                 D11_table2=D11_table+(bestScale)*D11_table;
%                 D12_table2 = D12_table+(bestScale)*D12_table ;
%                 D22_table2=D22_table+(bestScale)*D22_table;
%                 D33_table2=D33_table+(bestScale)*D33_table ;
%                 ps(1) = pstrain1(indexOfMinValue);
%                 ps(2)= pstrain2(indexOfMinValue);
%                 ps(3) = pstrain3(indexOfMinValue);
%                 if(flip==2)
%                     temp = ps(2);
%                     ps(2)= ps(1);
%                     ps(1)=temp;
%                 end
%                 % Make sure the shear sign is correct based on the target orientation.
%                 if(macroEleProps.Exx>macroEleProps.Eyy)
%                     ps(3)=abs( ps(3));
%                 else
%                     ps(3)=- abs(ps(3));
%                 end
%                 etaTargetLocal = originalEta+bestScale*originalEta;
%                 etaTargetLocal=max(config.MesoMinimumDensity,min(etaTargetLocal,1));
%                 fprintf('Targets  D11 %f D22 %f D33 %f\n', D11sys,D22sys,D33sys);
%                 fprintf('Expected D11 %f D22 %f D33 %f\nIndex is at %i\nbestscale %f\n%f\n', D11_table2,D22_table2,D33_table2,indexOfMinValue,bestScale,etaTargetLocal);
 elseif(version== 2)
                % ------------------------------------
                %           Version 2
                % ------------------------------------
                
                temp=config.targetTestVectorLen-1;
                term1 =0:1/temp:0.5; % -1 to 1 is the domain
                term2 =0:1/temp:0.5;  % -1 to 1 is the domain
                term3 =0:1/temp:0.3;  % to 1 is the domain
                densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
                
                t1=size(term1,2);
                t2=size(term2,2);
                t3=size(term3,2);
                t4=size(densityTargetsVector,2);
                
%                 D11=D11(1:end-2);
%                 D12=D12(1:end-2);
%                 D22=D22(1:end-2);
%                 D33=D33(1:end-2);
%                 etaTarget=etaTarget(1:end-2);
%                 pstrain1=pstrain1(1:end-2);
%                 pstrain2=pstrain2(1:end-2);
%                 pstrain3=pstrain3(1:end-2);
        
                D11 = reshape(D11,[t1 t2 t3 t4]);
                D12 = reshape(D12,[t1 t2 t3 t4]);
                D22 = reshape(D22,[t1 t2 t3 t4]);
                D33 = reshape(D33,[t1 t2 t3 t4]);
        
                etaTarget = reshape(etaTarget,[t1 t2 t3 t4]);
                pstrain1 = reshape(pstrain1,[t1 t2 t3 t4]);
                pstrain2 = reshape(pstrain2,[t1 t2 t3 t4]);
                pstrain3 = reshape(pstrain3,[t1 t2 t3 t4]);
        
                % logicArray = abs(strain1_before-strain1_after)<0.00001
        
        
        
                minValue=1e11;
                bestEtaForMinValue = 1;
                errorOnMinAllowed = 0.03;
                
                indexOfMinValue = 1;
                flip = 1;
                D11sys = max(macroEleProps.D_sys(1,1),0.001 );
                D12sys = max(macroEleProps.D_sys(1,2),0.001 );
                D22sys =max( macroEleProps.D_sys(2,2) ,0.001 );
                D33sys = max(macroEleProps.D_sys(3,3),0.001 );
                
                if(D22sys>D11sys)
                    flip = 2;
                else
                    flip=1;
                end
                
                for i= fliplr(1:t1) % scanning with D1 first
                    for j= fliplr(1:t2)
                        for k = 1:t3
                            for l = 1:t4
                                etaLocal = etaTarget(i,j,k,l);
                                if(etaLocal<-0.01)
                                    continue
                                end
                                if( etaLocal<config.MesoMinimumDensity )
                                    continue
                                end
                                
                                if(j>i)
                                    continue; % this data is not valid. since we didn't actually run it with knowing it is symmetric.
                                end
                                
                                D11_table = D11(i,j,k,l);
                                D12_table = D12(i,j,k,l);
                                D22_table = D22(i,j,k,l);
                                D33_table = D33(i,j,k,l);
                                %                                   for m = 1:2
                                % try flipping the D11 and D22
                                if(flip==2)
                                    temp = D22_table;
                                    D22_table=D11_table;
                                    D11_table=temp;
                                end
                                
                                diffValue = (D11_table-D11sys)^2+(D12_table-D12sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
                                
                                
                                if(diffValue<(minValue+minValue*errorOnMinAllowed)) % allow up to errorOnMinAllowed higher values.
                                    if(etaLocal<bestEtaForMinValue || diffValue<(minValue-minValue*errorOnMinAllowed)) % then compare the density
                                        minValue=diffValue;
                                        bestEtaForMinValue=etaLocal;
                                        indexOfMinValue=[i j k l];
%                                         flip = m;
                                    end
                                end
                                
                                %                                   end
                            end
                        end
                    end
                end
        
                  D11_table = D11(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
                D12_table = D12(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
                D22_table = D22(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
                D33_table = D33(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
                
%                  D11_table_densityHigher = D11(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),min(indexOfMinValue(4)+1,t4));
%                 D12_table_densityHigher = D12(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),min(indexOfMinValue(4)+1,t4));
%                 D22_table_densityHigher = D22(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),min(indexOfMinValue(4)+1,t4));
%                 D33_table_densityHigher = D33(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),min(indexOfMinValue(4)+1,t4));
%                 
%                   D11_table_densityLower = D11(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),max(indexOfMinValue(4)-1,1));
%                 D12_table_densityLower = D12(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),max(indexOfMinValue(4)-1,1));
%                 D22_table_densityLower = D22(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),max(indexOfMinValue(4)-1,1));
%                 D33_table_densityLower = D33(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),max(indexOfMinValue(4)-1,1));
                
                originalEta = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
        
                ps(1)=pstrain1(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
                ps(2)=pstrain2(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
                ps(3)=pstrain3(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
        
                if(flip==2)
                    temp = D22_table;
                    D22_table=D11_table;
                    D11_table=temp;
%                     
%                     temp = D22_table_densityHigher;
%                     D22_table_densityHigher=D11_table_densityHigher;
%                     D11_table_densityHigher=temp;
%                     
%                     temp = D22_table_densityLower;
%                     D22_table_densityLower=D11_table_densityLower;
%                     D11_table_densityLower=temp;
        
                       temp = ps(2);
                    ps(2)= ps(1);
                    ps(1)=temp;
                end
        
        
                if(macroEleProps.Exx>macroEleProps.Eyy)
                    ps(3)=abs( ps(3));
                else
                    ps(3)=- abs(ps(3));
                end
                
                
                scaleVersion=1;
                if(scaleVersion==1)
               
                    smallestDiff = minValue;
                    bestScale=0;
                    for scale = -1:0.001:1
                        D11_table2=D11_table+(scale)*D11_table;
                        D12_table2 = D12_table+(scale)*D12_table ;
                        D22_table2=D22_table+(scale)*D22_table;
                        D33_table2=D33_table+(scale)*D33_table ;
                        diffValue = (D11_table2-D11sys)^2+(D12_table2-D12sys)^2+(D22_table2-D22sys)^2+(D33_table2-D33sys)^2;
                        if(diffValue<smallestDiff)
                            bestScale=scale;
                            smallestDiff=diffValue;
                        end
                    end
                    D11_table2=D11_table+(bestScale)*D11_table;
                    %                 D12_table2 = D12_table+(bestScale)*D12_table ;
                    D22_table2=D22_table+(bestScale)*D22_table;
                    D33_table2=D33_table+(bestScale)*D33_table ;
                    
                    etaTargetLocal = originalEta+bestScale*originalEta;
                    etaTargetLocal=max(config.MesoMinimumDensity,min(etaTargetLocal,1));
                    
                else
%                     % --------------------
%                     % Version 2
%                     % --------------------
%                     
%                     % first scale from the lower to middle
%                     % then scale from middle to upper
%                     smallestDiffLower = minValue;
%                     bestScaleLower=0;
%                     for scale = 0:0.001:1
%                         weightLower=1-scale;
%                         D11_table2=D11_table_densityLower*weightLower+D11_table*scale;
%                         D12_table2=D12_table_densityLower*weightLower+D12_table*scale;
%                         D22_table2=D22_table_densityLower*weightLower+D22_table*scale;
%                         D33_table2=D33_table_densityLower*weightLower+D33_table*scale;
%                         
%                         
%                         diffValue = (D11_table2-D11sys)^2+(D12_table2-D12sys)^2+(D22_table2-D22sys)^2+(D33_table2-D33sys)^2;
%                         if(diffValue<smallestDiffLower)
%                             bestScaleLower=scale;
%                             smallestDiffLower=diffValue;
%                         end
%                     end
%                     
%                     % scale the upper
%                       smallestDiffUpper = minValue;
%                      bestScaleUpper=0;
%                     for scale = 0:0.001:1
%                         weightUpper=1-scale;
%                         D11_table2=D11_table_densityHigher*weightUpper+D11_table*scale;
%                         D12_table2=D12_table_densityHigher*weightUpper+D12_table*scale;
%                         D22_table2=D22_table_densityHigher*weightUpper+D22_table*scale;
%                         D33_table2=D33_table_densityHigher*weightUpper+D33_table*scale;
%                         
%                         
%                         diffValue = (D11_table2-D11sys)^2+(D12_table2-D12sys)^2+(D22_table2-D22sys)^2+(D33_table2-D33sys)^2;
%                         if(diffValue<smallestDiffUpper)
%                             bestScaleUpper=scale;
%                             smallestDiffUpper=diffValue;
%                         end
%                     end
%                     
%                     % check to see which scaleing was better. 
%                     if(smallestDiffLower<=smallestDiffUpper)
%                         etaLower = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),max(indexOfMinValue(4)-1,1));
%                         originalEta = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                         bestScale=-bestScaleLower+1;
%                         weightLower=1-bestScaleLower;
%                         etaTargetLocal=etaLower*weightLower+originalEta*bestScaleLower;
%                         
%                         scale=bestScaleLower;
%                         D11_table2=D11_table_densityLower*weightLower+D11_table*scale;
%                         D12_table2=D12_table_densityLower*weightLower+D12_table*scale;
%                         D22_table2=D22_table_densityLower*weightLower+D22_table*scale;
%                         D33_table2=D33_table_densityLower*weightLower+D33_table*scale;
%                         
%                     else
%                         etaUpper = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),min(indexOfMinValue(4)+11,t4));
%                         originalEta = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                         
%                         weightUpper=1-bestScaleUpper;
%                         etaTargetLocal=etaUpper*weightUpper+originalEta*bestScaleUpper;
%                         
%                         bestScale=1-bestScaleUpper;
%                         scale=bestScaleUpper;
%                         
%                          D11_table2=D11_table_densityHigher*weightUpper+D11_table*scale;
%                         D12_table2=D12_table_densityHigher*weightUpper+D12_table*scale;
%                         D22_table2=D22_table_densityHigher*weightUpper+D22_table*scale;
%                         D33_table2=D33_table_densityHigher*weightUpper+D33_table*scale;
%                         %                         etaTargetLocal=etaUpper*(1-bestScaleUpper)+bestScaleUpper*originalEta;
%                         
%                     end
%                     t=1;
                    
                end  
            
                fprintf('Targets  D11 %f D22 %f D33 %f\n', D11sys,D22sys,D33sys);
                i=indexOfMinValue;
                fprintf('Expected D11 %f D22 %f D33 %f\nIndex is at %i,%i,%i,%i\nbestscale %f\nEta %f\n', D11_table2,D22_table2,D33_table2,i(1),i(2),i(3),i(4),bestScale,etaTargetLocal);
            elseif(version== 3)
%                 % ------------------------------------
%                 %           Version 3
%                 % ------------------------------------
%         
%                 temp=config.targetTestVectorLen-1;
%                 term1 =0:1/temp:0.5; % -1 to 1 is the domain
%                 term2 =0:1/temp:0.5;  % -1 to 1 is the domain
%                 term3 =0:1/temp:0.3;  % to 1 is the domain
%                 densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
%         
%                 t1=size(term1,2);
%                 t2=size(term2,2);
%                 t3=size(term3,2);
%                 t4=size(densityTargetsVector,2);
%         
%                 D11=D11(1:end-2);
%                 D12=D12(1:end-2);
%                 D22=D22(1:end-2);
%                 D33=D33(1:end-2);
%                 etaTarget=etaTarget(1:end-2);
%                 pstrain1=pstrain1(1:end-2);
%                 pstrain2=pstrain2(1:end-2);
%                 pstrain3=pstrain3(1:end-2);
%         
%                 D11 = reshape(D11,[t1 t2 t3 t4]);
%                 D12 = reshape(D12,[t1 t2 t3 t4]);
%                 D22 = reshape(D22,[t1 t2 t3 t4]);
%                 D33 = reshape(D33,[t1 t2 t3 t4]);
%         
%                 etaTarget = reshape(etaTarget,[t1 t2 t3 t4]);
%                 pstrain1 = reshape(pstrain1,[t1 t2 t3 t4]);
%                 pstrain2 = reshape(pstrain2,[t1 t2 t3 t4]);
%                 pstrain3 = reshape(pstrain3,[t1 t2 t3 t4]);
%         
%                 % logicArray = abs(strain1_before-strain1_after)<0.00001
%         
%         
%         
%                 minValue=1e11;
%                 indexOfMinValue = 1;
%                 flip = 1;
%                 D11sys = max(macroEleProps.D_sys(1,1),0.001 );
%                 D12sys = max(macroEleProps.D_sys(1,2),0.001 );
%                 D22sys =max( macroEleProps.D_sys(2,2) ,0.001 );
%                 D33sys = max(macroEleProps.D_sys(3,3),0.001 );
%         
%                 for i= 1:t1 % scanning with D1 first
%                     for j= 1:t2
%                         for k = 1:t3
%                             for l = 1:t4
%                                 etaLocal = etaTarget(i,j,k,l);
%                                 if(etaLocal<-0.01)
%                                     continue
%                                 end
%         
%                                 D11_table = D11(i,j,k,l);
%                                 D12_table = D12(i,j,k,l);
%                                 D22_table = D22(i,j,k,l);
%                                 D33_table = D33(i,j,k,l);
%                                 for m = 1:2
%                                     % try flipping the D11 and D22
%                                     if(m==2)
%                                         temp = D22_table;
%                                         D22_table=D11_table;
%                                         D11_table=temp;
%                                     end
%         
%                                     diffValue = (D11_table-D11sys)^2+(D12_table-D12sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
%                                     if(diffValue<minValue && etaLocal>config.MesoMinimumDensity )
%                                         minValue=diffValue;
%                                         indexOfMinValue=[i j k l];
%                                         flip = m;
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%         
%                   D11_table = D11(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                 D12_table = D12(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                 D22_table = D22(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                 D33_table = D33(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                 originalEta = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         
%                 ps(1)=pstrain1(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                 ps(2)=pstrain2(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%                 ps(3)=pstrain3(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         
%                 if(flip==2)
%                     temp = D22_table;
%                     D22_table=D11_table;
%                     D11_table=temp;
%         
%                        temp = ps(2);
%                     ps(2)= ps(1);
%                     ps(1)=temp;
%                 end
%         
%         
%                 if(macroEleProps.Exx>macroEleProps.Eyy)
%                     ps(3)=abs( ps(3));
%                 else
%                     ps(3)=- abs(ps(3));
%                 end
%                 ps
%                 t=1;
%         
%                 % Generate an interpolation and then find the min
%                 start1 = max(1,indexOfMinValue(1)-1);
%                 last1 = min(t1,indexOfMinValue(1)+1);
%         
%                 start2 = max(1,indexOfMinValue(2)-1);
%                 last2 = min(t2,indexOfMinValue(2)+1);
%         
%                 start3 = max(1,indexOfMinValue(3)-1);
%                 last3 = min(t3,indexOfMinValue(3)+1);
%         
%                 start4 = max(1,indexOfMinValue(4)-1);
%                 last4 = min(t4,indexOfMinValue(4)+1);
%         
%                 vector1 = start1:last1 ;
%                  vector2 = start2:last2;
%                   vector3 =start3:last3;
%                    vector4 = start4:last4;
%         
%         
%                 [x,y,z,t] = ndgrid(vector1,vector2,vector3,vector4);
%                 diffArray = x*0;
%                    etaArray = x*0;
%                    p1_array= x*0;
%                    p2_array= x*0;
%                    p3_array= x*0;
%         
%         %          for i= vector1% scanning with D1 first
%         %             for j= vector2
%         %                 for k =vector3
%         %                     for l = vector4
%         
%                  for i=1: size(x,1)% scanning with D1 first
%                     for j= 1:size(x,2)
%                         for k =1:size(x,3)
%                             for l =1: size(x,4)
%         
%                                 xpoint = x(i,j,k,l);
%                                 ypoint = y(i,j,k,l);
%                                 zpoint = z(i,j,k,l);
%                                 tpoint = t(i,j,k,l);
%         
%                                   etaLocal = etaTarget(xpoint,ypoint,zpoint,tpoint);
%                                 if(etaLocal<-0.01)
%                                     diffArray(i,j,k,l)=minValue+1e3;
%                                         etaArray(i,j,k,l)=0;
%         
%                                     continue
%                                 end
%         
%         %                         D11_table = D11(i,j,k,l);
%         %                         D12_table = D12(i,j,k,l);
%         %                         D22_table = D22(i,j,k,l);
%         %                         D33_table = D33(i,j,k,l);
%         %
%         
%                                 D11_table = D11(xpoint,ypoint,zpoint,tpoint);
%                                 D12_table = D12(xpoint,ypoint,zpoint,tpoint);
%                                 D22_table = D22(xpoint,ypoint,zpoint,tpoint);
%                                 D33_table = D33(xpoint,ypoint,zpoint,tpoint);
%         
%                                 %                         for m = 1:2
%                                 % try flipping the D11 and D22
%                                 if(flip==2)
%                                     temp = D22_table;
%                                     D22_table=D11_table;
%                                     D11_table=temp;
%                                 end
%         
%                                 diffValue = (D11_table-D11sys)^2+(D12_table-D12sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
%                                 diffArray(i,j,k,l)=diffValue;
%         
%                                 etaArray(i,j,k,l)=etaLocal;
%                                 p1_array(i,j,k,l)=pstrain1(xpoint,ypoint,zpoint,tpoint);
%                                 p2_array(i,j,k,l)=pstrain2(xpoint,ypoint,zpoint,tpoint);
%                                 p3_array(i,j,k,l)=pstrain3(xpoint,ypoint,zpoint,tpoint);
%                                 %                         if(diffValue<minValue && etaLocal>config.MesoMinimumDensity )
%                                 %                             minValue=diffValue;
%                                 %                             indexOfMinValue=[i j k l];
%                                 %                             flip = m;
%                                 %                         end
%         %                         end
%                             end
%                         end
%                     end
%                  end
%                  maxOfAll= max(max(max(max(diffArray))));
%                  step=0.1;
%                  vector1_v2 = start1:step:last1 ;
%                  vector2_v2 = start2:step:last2;
%                  vector3_v2 =start3:step:last3;
%                  vector4_v2 = start4:step:last4;
%         
%         
%                  [xq,yq,zq,tq] = ndgrid(vector1_v2,vector2_v2,vector3_v2,vector4_v2);
%                  diffArray_Q = interpn(x,y,z,t,diffArray,xq,yq,zq,tq);
%         
%                  Eta_q = interpn(x,y,z,t,etaArray,xq,yq,zq,tq);
%         
%                  p1_q = interpn(x,y,z,t,p1_array,xq,yq,zq,tq);
%                  p2_q = interpn(x,y,z,t,p2_array,xq,yq,zq,tq);
%                  p3_q = interpn(x,y,z,t,p3_array,xq,yq,zq,tq);
%         
%         
%                  %            figure('renderer','zbuffer');
%                  %            nframes = size(tq, 4);
%                  %            for j = 1:nframes
%                  %                slice(yq(:,:,:,j),xq(:,:,:,j),zq(:,:,:,j),...
%                  %                    Vq(:,:,:,j),0,0,0);
%                  %                caxis([0 maxOfAll]);
%                  %                M(j) = getframe;
%                  %            end
%                  %            movie(M);
%                  makevide=1;
%                  if(makevide ==1)
%                      config.recvid=1;
%                      video = VideoManager;
%                      [vidObj, framedNumber] = video.InitializeVideo( config);
%                      F=getframe();
%                  end
%         
%                  maxOfAll= max(max(max(max(diffArray))));
%                  figure('renderer','zbuffer');
%                  nframes = size(diffArray_Q, 4);
%                  for j = 1:nframes
%                      %              slice(x(:,:,:,j),y(:,:,:,j),z(:,:,:,j),...
%                      %                  diffArray(:,:,:,j),0,0,0);
%                      %              caxis([0 maxOfAll]);
%                      %              RhoColumn=  diffArray(:,:,:,j); % color
%         
%                      %
%                      %              xtemp = reshape(xq(:,:,:,j),[],1);
%                      %              ytemp = reshape(yq(:,:,:,j),[],1);
%                      %              ztemp = reshape(zq(:,:,:,j),[],1);
%         
%                      xtemp = reshape(p1_q(:,:,:,j),[],1);
%                      ytemp = reshape(p2_q(:,:,:,j),[],1);
%                      ztemp = reshape(p3_q(:,:,:,j),[],1);
%                      RhoColumn = reshape(diffArray_Q(:,:,:,j),[],1);
%                      circleSize = ones(size(RhoColumn))*100; % circle size.
%                      %              scatter3(x(:,:,:,j),y(:,:,:,j),z(:,:,:,j),circleSize, RhoColumn,'filled','MarkerEdgeColor','k')
%                      scatter3(xtemp,ytemp,ztemp,circleSize, RhoColumn,'filled','MarkerEdgeColor','k')
%                      caxis([0 maxOfAll]);
%                      %                scatter3(xtemp,ytemp,ztemp,circleSize, RhoColumn,'filled','MarkerEdgeColor','k')
%                      %                  scatter3(xq(:,:,:,j),vq(:,:,:,j),zq(:,:,:,j),   Vq(:,:,:,j),circleSize, RhoColumn,'filled','MarkerEdgeColor','k')
%         
%                      % Zoomed out view.
%                      %              xlim([0 1]);
%                      %                ylim([0 1]);
%                      %                  zlim([0 1]);
%                      axis equal
%                      targetDensity=Eta_q(1,1,1,j);
%                      title(sprintf('diffArray as function of psuedo strains at density = %f',targetDensity));
%                      colorbar
%                      xlabel('pstrain1');
%                      ylabel('pstrain2');
%                      zlabel('pstrain3');
%                      %              colormap winter
%                      if(makevide ==1)
%                          [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
%                      end
%                  end
%                  %          movie(M);
%                  if(makevide ==1)
%                      video.CloseVideo( config, F,vidObj)
%                  end
        
        
    elseif(version== 4)
        % ------------------------------------
        %           Version 4
        % ------------------------------------
        
        
        temp=config.targetTestVectorLen-1;
        term1 =0:1/temp:0.5; % -1 to 1 is the domain
        term2 =0:1/temp:0.5;  % -1 to 1 is the domain
        term3 =0:1/temp:0.3;  % to 1 is the domain
        densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
        
        
        
        t1=size(term1,2);
        t2=size(term2,2);
        t3=size(term3,2);
        t4=size(densityTargetsVector,2);
        
        %         D11=D11(1:end-2); % remove the 2 extra valeues that I added.
        %         D12=D12(1:end-2);
        %         D22=D22(1:end-2);
        %         D33=D33(1:end-2);
        %         etaTarget=etaTarget(1:end-2);
        %         pstrain1=pstrain1(1:end-2);
        %         pstrain2=pstrain2(1:end-2);
        %         pstrain3=pstrain3(1:end-2);
        
        % flip everything. This
        
        % flip over the symmetric axis of ps1=ps2
        %         logictest = pstrain2>pstrain1;
        %         D11(logictest)=
        
        
        D11_table4D = reshape(D11,[t1 t2 t3 t4]);
        D12_table4D = reshape(D12,[t1 t2 t3 t4]);
        D22_table4D = reshape(D22,[t1 t2 t3 t4]);
        D33_table4D = reshape(D33,[t1 t2 t3 t4]);
        
        eta_table4D = reshape(etaTarget,[t1 t2 t3 t4]);
        ps1_table4D = reshape(pstrain1,[t1 t2 t3 t4]);
        ps2_table4D = reshape(pstrain2,[t1 t2 t3 t4]);
        ps3_table4D = reshape(pstrain3,[t1 t2 t3 t4]);
        
        
        
        
        %           temp=config.targetTestVectorLen-1;
        %     term1 =0:1/temp:0.5; % -1 to 1 is the domain
        %     term2 =0:1/temp:0.5;  % -1 to 1 is the domain
        %     term3 =0:1/temp:0.3;  % to 1 is the domain
        %     densityTargetsVector = 0:1/temp:1;  % 0 to 1 is the domain
        %     [strain1_temp, strain2_temp, strain3_temp, densityTargets_temp] = ndgrid(term1,term2,term3,densityTargetsVector);
        %
        %      strain1_temp2 = reshape(strain1_temp,[],1);
        %
        %
        %         logicArray = abs(ps1_table4D-strain1_temp)<0.00001
        
        
        
        minValue=1e11;
        indexOfMinValue = 1;
        flip = 1;
        D11sys = max(macroEleProps.D_sys(1,1),0.001 );
        D12sys = max(macroEleProps.D_sys(1,2),0.001 );
        D22sys =max( macroEleProps.D_sys(2,2) ,0.001 );
        D33sys = max(macroEleProps.D_sys(3,3),0.001 );
        
        
        % 1 reshape everything (done)        %
        % -- use ngrid to generate 4d arrays of the x,y,z,t indixes, use
        % these indixes when getting the data ans saving it.
        % 3. call fmincon
        % -- find the psuedo strains and densiteis that minimize the diffvalue
        % -- inputs to evalutor are pseudo strain values and density as 4D arrays, and the psuedo strain and density to evalute and the system values.
        % --------- interpn, the psuedo strains and density, then uses the interped values to find the diffValue
        % -- outputs of the evalutor are the diffValue (which is interpolated based on the data we have), the pseeudo strains and density target
        %
        
        % -----------
        % Find a start value
        % ------------------
%         for i= 1:t1 % scanning with D1 first
%             for j= 1:t2
%                 for k = 1:t3
%                     for l = 1:t4
%                         %                         etaLocal = eta_table4D(i,j,k,l);
%                         ps1= ps1_table4D(i,j,k,l);
%                         ps2 = ps2_table4D(i,j,k,l);
%                         if(ps2>ps1)
%                             continue
%                         end
%                         
%                         D11_table = D11_table4D(i,j,k,l);
%                         D12_table = D12_table4D(i,j,k,l);
%                         D22_table = D22_table4D(i,j,k,l);
%                          D33_table = D33_table4D(i,j,k,l);
%                         %                         for m = 1:2
%                         %                             % try flipping the D11 and D22
%                         %                             if(m==2)
%                         %                                 temp = D22_table;
%                         %                                 D22_table=D11_table;
%                         %                                 D11_table=temp;
%                         %                             end
%                         
%                         diffValue = (D11_table-D11sys)^2+(D12_table-D12sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
%                         if(diffValue<minValue  )
%                             minValue=diffValue;
%                             indexOfMinValue=[i j k l];
%                             %                             flip = m;
%                         end
%                         %                         end
%                     end
%                 end
%             end
%         end
%         %
%         %           D11_table = D11(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         %         D12_table = D12(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         %         D22_table = D22(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         %         D33_table = D33(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         %         originalEta = etaTarget(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         %
%         ps_start(1)=ps1_table4D(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         ps_start(2)=ps2_table4D(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         ps_start(3)=ps3_table4D(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
%         eta_start=eta_table4D(indexOfMinValue(1),indexOfMinValue(2),indexOfMinValue(3),indexOfMinValue(4));
        %
        %         if(flip==2)
        %             temp = D22_table;
        %             D22_table=D11_table;
        %             D11_table=temp;
        %
        %                temp = ps(2);
        %             ps(2)= ps(1);
        %             ps(1)=temp;
        %         end
        
        
        %
        % x has 4 values, 3 psuedo strains, and 1 etaTarget
        %         x0=[0.5 0.5 0.5 0.5];
        
        A = [];
        b = [];
        
        
        
        ub = [0.5 0.5 0.3 1];
        lb = [0 0 0 0];
                 x0=(ub+lb)/2;
%         x0=[  ps_start(1)  ps_start(2)  ps_start(3) eta_start];
        
        [x finalDiffValue]= fmincon(@(x) EvalutePseudoStrainAndDensityForFit(x,D11sys,D12sys,D22sys,D33sys, ...
            D11_table4D,D12_table4D,D22_table4D,D33_table4D,...
            ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D)...
            ,x0,A,b,[],[],lb,ub);
        
        ps(1)=x(1);
        ps(2)=x(2);
        ps(3)=x(3);
        etaTargetLocal=x(4);
        
        %         [D11_interp,D12_interp ,D22_interp,D33_interp]= InterpolatePseudoStrainsAndDensity(ps_new,eta_new, ...
        %     D11_table4D,D12_table4D,D22_table4D,D33_table4D,...
        %     ps1_table4D,ps2_table4D,ps3_table4D,eta_table4D)
        
        
        if(D22sys>D11sys)
            %             temp = D22_interp;
            %             D22_interp=D11_interp;
            %             D11_interp=temp;
            
            temp = ps(2);
            ps(2)= ps(1);
            ps(1)=temp;
        end
        
        if(macroEleProps.Exx>macroEleProps.Eyy)
            ps(3)=abs( ps(3));
        else
            ps(3)=- abs(ps(3));
        end
        
%         ps_start
        ps
        
%         eta_start
        etaTargetLocal
        
        
    end
    
    
    
    %     outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',0,1,indexOfMinValue);
    %     d = csvread(outname)
    
    
    %      macroEleProps.psuedoStrain(1) = ps(1) ;
    %     macroEleProps.psuedoStrain(2) = ps(2);
    %     macroEleProps.psuedoStrain(3) =  ps(3);
    macroEleProps.psuedoStrain=[ps(1);ps(2);ps(3)];
    macroEleProps.targetDensity=etaTargetLocal;
    
    if 1==0
        start = max(indexOfMinValue-500,1);
        limit = min(indexOfMinValue+500,t2);
        figure(1)
        subplot(7,1,1);
        plot(start:limit,D11(start:limit))
        ylabel('d11');
        
        subplot(7,1,2);
        plot(start:limit,D22(start:limit))
        ylabel('d22');
        
        subplot(7,1,3);
        plot(start:limit,D33(start:limit))
        ylabel('d33');
        
        subplot(7,1,4);
        plot(start:limit,pstrain1(start:limit))
        ylabel('p1');
        
        subplot(7,1,5);
        plot(start:limit,pstrain2(start:limit))
        ylabel('p2');
        
        subplot(7,1,6);
        plot(start:limit,pstrain3(start:limit))
        ylabel('p3');
        
        subplot(7,1,7);
        plot(start:limit,etaTarget(start:limit))
        ylabel('eta');
        hold on
        stairs([start indexOfMinValue limit],[0 1 1] );
        hold off
        
        
    end
    
    
    
    
end
% else
%
%     % -----------------------------------
%     %
%     % More than 1 element per design var case on macro level
%     % e = design var and not element for this case.
%     %
%     % -----------------------------------
%
%
%     %
%     %1. make an array of meso local X positions , and make another one for Y positions for each node that is a part of the elements
%     %2. make an zeros array that will hold the X displacements at each
%     %3. make an zeros array that will hold the Y displacements at each
%     %4. make a list of macro element node numbers
%     %5. Loop over the macro elements
%     %	1. check if this node has already been added to the lsit of node numbers
%     %	2. if not then add the X and Y displacements to the local displacement array
%
%     %     num Xnodes = config.numXElmPerDV
%     %      [X,Y] = ndgrid(0:config.numXElmPerDV,0:config.numYElmPerDV);
%     [Y,X] = ndgrid(0:config.numYElmPerDV,0:config.numXElmPerDV);
%     %     [t1, t2] = size(X);
%     %     displacementsX = zeros(t1,t2);
%     %     displacementsY = displacementsX;
%
%     macroEleProps.mesoXnodelocations=X;
%     macroEleProps.mesoYnodelocations=Y;
%
%     numVarsinRow =config.numVarsX;
%     %     numVarsinColumn =config.numVarsY;
%     ydesignVar =   floor( (e-1)/numVarsinRow)+1;
%     xdesignVar = mod(e-1,numVarsinRow)+1;
%
%     % trying to find the nodes and then dofs of the macro problem
%     shiftElementsX = (xdesignVar-1)*config.numXElmPerDV;
%     elementsInRow = config.nelx;
%     shiftNodesY = (ydesignVar-1)*config.numYElmPerDV*elementsInRow;
%
%     %     count1 = 1;
%     %     count2 = 1;
%     %     isfirstrow = 1;
%
%     nodelist=[];
%     for i  =1: config.numYElmPerDV
%         startElement = shiftElementsX+shiftNodesY+(i-1)*elementsInRow;
%         for j =1: config.numXElmPerDV
%             ee = startElement+j;
%
%             nodes1=  IEN(ee,:);
%             nodelist = [nodelist nodes1];
%         end
%     end
%
%     nodelist = unique(sort(nodelist));
%
%     % nodes1=  IEN(e,:);
%     macroEleProps.elementNodes=nodelist;
%
%     % multiloading cases.
%     [~, t2] = size(config.loadingCase);
%     for loadcaseIndex = 1:t2
%         utemp = U(loadcaseIndex,:);
%
%         displacementsX = utemp(nodelist*2-1);
%         displacementsY = utemp(nodelist*2);
%
%         macroEleProps.xDisplacements(loadcaseIndex,:)=displacementsX;
%         macroEleProps.yDisplacements(loadcaseIndex,:)=displacementsY;
%     end
%
%     macroEleProps.elementNumber = e;
%     % Save element to XY position map (needed for x and w vars retrival)
%     %         outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
%     %         elementXYposition=csvread(outname);
%     %         results = elementXYposition(macroEleProps.elementNumber,:);
%     macroEleProps.yPos = ydesignVar;
%     macroEleProps.xPos = xdesignVar;
%
%     % Get the density field
%     outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
%     x = csvread(outname);
%     macroEleProps.density = x(macroEleProps.yPos,macroEleProps.xPos );
%
%     % % Get the volume fraction field
%     outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
%     w = csvread(outname);
%     macroEleProps.material1Fraction = w(macroEleProps.yPos,macroEleProps.xPos );
%
%     % if(macro_meso_iteration>1)
%     outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
%     D = csvread(outname);
%     macroEleProps.D_given =D;
%     % end
%
%
%
% end