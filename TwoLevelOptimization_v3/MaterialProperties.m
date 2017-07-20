classdef MaterialProperties
    properties
        % steel AISI 2080 is material 1
        % copper is material 2
        %         E_material1 = 200000; % N/mm^2 The elastic mod of material 1
        %         E_material2 = 110000; % The elastic mod of material 2
        %
        %         K_material1 = 47/1000; %  W/ (mm*K)heat conduction of material 1
        %         K_material2 = 390/1000; % heat conduction of material 2
        %
        %         alpha1 = 1.5e-5; %thermal expansion coefficient for material 1
        %         alpha2 = 2.4e-5; % thermal expansion coefficient for material 2
        
        E_material1 =100000;% 100000; %4 N/mm^2 The elastic mod of material 1
        E_material2 = 25000;% 25000 E_material1/2; %4 The elastic mod of material 2
%          E_material1 =1;% 100000; %4 N/mm^2 The elastic mod of material 1
%         E_material2 = 0.25;% 25000 E_material1/2; %4 The elastic mod of material 2
        
        K_material1 = 0.02; %0.02  W/ (mm*K)heat conduction of material 1
        K_material2 = 0.04; %0.04 heat conduction of material 2
        
        alpha1 =0.001; %thermal expansion coefficient for material 1
        alpha2 =0.001; % thermal expansion coefficient for material 2
        
        v = 0.3; % Piossons ratio
        G ; % = E_material1/(2*(1+v));
        dKelastic; % derivative of K elastic matrix with respect to a material vol fraction change.
        dKheat;   % derivative of K heat matrix with respect to a material vol fraction change.
        
        
        
        SavedDmatrix ;%= zeros(ne,9);
    end
    
    methods
        
        % Constructor
        function obj = MaterialProperties
            obj.G = obj.E_material1/(2*(1+obj.v));
            E = 1;
            strain = []; % empty matrix,
            Dgiven = []; % empty
            obj.dKelastic = elK_elastic(E,obj.v, obj.G,strain,Dgiven)*(obj.E_material1-obj.E_material2);
            
            heatCoefficient = 1;
            obj.dKheat =  elementK_heat(heatCoefficient)*(obj.K_material1-obj.K_material2);
        end
        
        %% -------------------------------------------------------
        % w is vol fraction of material 1, x is density, at this meso
        % element or design var.
        % -------------------------------------------------------
        function [density]=CalculateDensityTargetforMeso(obj,w,x,Exx,Eyy,config)
            
            % if using the new method, then use a different equation.
            if(config.useExxEyy==1)
%                 density=0.5*(Exx+Eyy)/obj.E_material1;
                 density=max(Exx,Eyy)/obj.E_material1;
            else
                offset = obj.E_material2/obj.E_material1;
                density = offset+(w)*(1-offset);
                
            end
        end
        
        %% ------------------------------------------------
        %
        %             READS THE D MATRIXES AND SAVES THEM TO AN ARRAY FOR
        %             FUTURE USE.
        %
        % ------------------------------------------------
        function obj =  ReadConstitutiveMatrixesFromFiles(obj,  config)
            if(config.macro_meso_iteration>1 || config.mode==90)
                folderNum = config.iterationNum;
                if(config.mode==90)
                      oldIteration = config.macro_meso_iteration;
                else
                    oldIteration = config.macro_meso_iteration-1; % minus 1, because we want to get the previous iterationdesign.
                end
                
                ne = config.nelx*config.nely;
                obj.SavedDmatrix = zeros(ne,9);
                
                % Get the density field
                outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,oldIteration);
                xxx = csvread(outname);
                
                outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,oldIteration);
                elementXYposition=csvread(outname);
                
                
                for e = 1:ne
                    
                    results = elementXYposition(e,:);
                    yPos = results(1);
                    xPos = results(2);
                    densitySIMP =xxx(yPos,xPos);
                    
                    if(densitySIMP>config.voidMaterialDensityCutOff)
                        % Read the D_h
                        %                     outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,elementNumber);
                        outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,oldIteration,e);
                        if exist(outname, 'file') == 2
                            try
                            D_h= csvread(outname);
                            catch
                               fprintf('could not read, %s\n',outname);
                               continue
                            end
                        else
                            material1Fraction=1;
                            D_h = obj.calculateEffectiveConstitutiveEquation(material1Fraction, config,[]);
                        end
                     else
                          material1Fraction=1;
                          topDensity=densitySIMP;
                          Exx=1;
                          Eyy=1;
                          rotation=1;
                          E12=1;
                          E33=1;
                            D_h = obj.getDmatMatrixTopExxYyyRotVars(config,topDensity,Exx, Eyy,rotation,material1Fraction, E12, E33);
                     end
                    
                    % Read the D_old
%                     outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,oldIteration,e);
%                     if exist(outname, 'file') == 2
%                         D_old= csvread(outname);
%                     else
%                         material1Fraction=1;
%                         D_old = obj.calculateEffectiveConstitutiveEquation(material1Fraction, config,[]);
%                     end
                    
                    % --------------------
                    %----------------- OLID [ Average the 2 Ds (D_h and                    %D_old]
                    % --------------------
                    D = D_h ; %(D_h+D_old)/2;
                    
                    
                    D_flat = reshape(D,1,9);
                    obj.SavedDmatrix(e,:)=D_flat;
                end
            end
        end
        
        % --------------------------------------------
        % GETS A SAVED D MATRIX FROM THE SAVED ARRAY
        % --------------------------------------------
        function [D] = GetSavedDMatrix(obj,e)
            D_flat =  obj.SavedDmatrix(e,:);
            D = reshape(D_flat,3,3);
        end
        
        
        %% ---------------------------------------
        %
        %
        % 	ELASTIC CALCULATIONS
        %
        % ---------------------------------------
        
        
        
        % ---------------------------
        %
        % CALCULATE K MATRIX
        %    USE TOPOLOGY VAR
        %    USE Exx and Eyy
        %  HAVE SPECIAL CASES FOR OLD METHODS AS LEGACY CODE
        %  The material gradient is for legacy.
        % ---------------------------
        function K = getKMatrixTopExxYyyRotVars(obj,config,topDensity,Exx, Eyy,rotation,material1Fraction,E12, E33,e)
            
            % For validation
            if(config.validationModeOn==1 || config.mode==90)
                if(config.macro_meso_iteration>1 || config.mode==90)
                     D_out= obj.GetSavedDMatrix(e);

                     if(config.multiscaleMethodCompare==1)
                         D_out=D_out*topDensity;
                     end

                 else
                      D_out = [ 1 obj.v 0;
                                obj.v 1 0;
                                0 0 1/2*(1-obj.v)]*obj.E_material1/(1-obj.v^2);
                end
                [K]=elK_elastic_v2(D_out);
                return;
            end
            
            
            % -------------
            % New method
            % -------------
            if(config.useExxEyy==1 || config.useRotation==1 || config.anisotropicMat)
                D_out= obj.getDmatMatrixTopExxYyyRotVars(config,topDensity,Exx, Eyy,rotation,material1Fraction,E12, E33);
                [K]=elK_elastic_v2(D_out);
                return;
            else % END NEW METHOD
                % -------------
                % Old Method method
                % -------------
                D_out= obj.getDmatMatrixTopExxYyyRotVars(config,topDensity,Exx, Eyy,rotation,material1Fraction);
                [K, ~, ~, ~]=elK_elastic(1,obj.v, obj.G,[],D_out);
            end% END OLD METHOD
        end
        
        function D_out =  getDmatMatrixTopExxYyyRotVars(obj,config,topDensity,Exx, Eyy,rotation,material1Fraction, E12, E33)
            % -------------
            % New method
            % -------------
            
            if((config.useExxEyy==1 || config.useRotation==1 ) && config.anisotropicMat ~=1) 
                D_out = obj.getDmatrixForExxEyy(Exx, Eyy,config); % Exx and Eyy DISTRIBUTION
                D_out=D_out*topDensity^(config.penal); % TOPOLOGY;
                if(config.useRotation==1)  % ROTATION
                    D_out=obj.rotateDmatrix(config,rotation, D_out)  ;
                end
                
                % - anisotropic Case
            elseif(config.anisotropicMat ==1)
                 D_out = obj.getDmatrixForAnisotropic(Exx, Eyy, E12, E33); 
                D_out=D_out*topDensity^(config.penal); % TOPOLOGY;

            else % END NEW METHOD
                
                % -------------
                % Old Method method
                % -------------
                E0 = effectiveElasticProperties(obj, material1Fraction, config); % GRADIENT MATERIAL
                E0 = E0*topDensity^(config.penal); % TOPOLOGY
                if(config.macro_meso_iteration>1)
                    D_out =matProp.GetSavedDMatrix(count);
                else                   
                    D_out = [ 1 obj.v 0;
                                obj.v 1 0;
                                0 0 1/2*(1-obj.v)]*E0/(1-obj.v^2);
                end
            end% END OLD METHOD
        end
        
        %         function D_out = getDmatrixforElement(obj,config,topDensity,material1Fraction,orthD,rotation)
        %             %  -------------
        %             % New method
        %             % -------------
        %              E0 = effectiveElasticProperties(obj, material1Fraction, config); % GRADIENT MATERIAL
        %             E0 = E0*topDensity^(config.penal);
        %             if(config.useOrthDistribution==1 || config.useRotation==1 )
        %                 D_out = obj.getDmatrixForOrthDistrVar(orthD);
        %                 D_out=D_out*E0;
        %                 % Add rotation here!!!!!
        %                 if(config.useRotation==1)
        %                     D_out=obj.rotateDmatrix(config,rotation, D_out)  ;
        %                 end
        %             else % END NEW METHOD
        %                 % -------------
        %                 % Old Method method
        %                 % -------------
        %                 if(config.macro_meso_iteration>1)
        %                     D_out =matProp.GetSavedDMatrix(count);
        %                 else
        %
        %                     D_out =   [ 1 obj.v 0;
        %                                 obj.v 1 0;
        %                                 0 0 1/2*(1-obj.v)]*E0/(1-obj.v^2);
        %                 end
        %             end% END OLD METHOD D = reshape(D_flat,3,3);
        %         end
        
        % ---------------------------
        %
        % CALCULATE K MATRIX FOR Exx,Eyy,theta, when using the R ratio to
        % scale the shear components.
        %
        % Mode = 1, Exx sensitivity
        % Mode = 2, Eyy sensitivity
        %
        % ---------------------------
        function K = getKMatrixSensitivityTopExxYyyRotVarsWithRshear(obj,config,topDensity,Exx, Eyy,rotation,mode)
            term1 = obj.v/(1-obj.v^2);
            avgE = 0.5*(Exx+Eyy);
            if(Exx>Eyy) % Case 1, where the Exx value is greater
                if(mode==1) % Exx derivitive mode
                    dE12_dExx = term1*(0.5*Eyy/Exx+avgE*Eyy*(-1)*Exx^-2);
                    dE33_dExx=0.5*(1-obj.v)/(1-obj.v^2)*(-1*Eyy/Exx^2*avgE+0.5*Eyy/Exx);
                else                    
                    dE12_dEyy = term1*(0.5*Eyy/Exx+avgE*1/Exx);
                    dE33_dEyy=0.5*(1-obj.v)/(1-obj.v^2)*(1/Exx*avgE+0.5*Eyy/Exx);
                end
                
            else % Case 2, where the Eyy value is greater
                if(mode==1) % Exx derivitive mode
                    dE12_dExx = term1*(0.5*Eyy/Exx+avgE*1/Exx);
                    dE33_dExx=0.5*(1-obj.v)/(1-obj.v^2)*(1/Eyy*avgE+0.5*Exx/Eyy);
                else
                     dE12_dEyy = term1*(0.5*Exx/Eyy+avgE*Exx*(-1)*Eyy^-2);
                      dE33_dEyy=0.5*(1-obj.v)/(1-obj.v^2)*(-1*Exx/Eyy^2*avgE+0.5*Exx/Eyy);
                end                
            end
            
              if(mode==1) % Exx derivitive mode
                  D_out = [ 1/(1-obj.v^2)   dE12_dExx     0;
                            dE12_dExx       0             0;
                            0               0             dE33_dExx];
              else  % Eyy derivitive mode
                   D_out = [ 0             dE12_dEyy      0;
                            dE12_dEyy      1/(1-obj.v^2)  0;
                            0               0             dE33_dEyy];
              end
              
               if(config.useRotation==1)  % ROTATION
                    D_out=obj.rotateDmatrix(config,rotation, D_out)  ;
               end
                D_out=D_out*topDensity^config.penal;
            
            
            [K]=elK_elastic_v2(D_out);
            
        end
            
        
        
        % ---------------------------
        %
        % CALCULATE K MATRIX FOR GRADIENT MATERIAL SENSITIVITY
        %
        % ---------------------------
        function K = getKMatrixGradientMaterialSensitivity(obj,config,topDensity,Exx, Eyy,rotation)
            
            % -------------
            % New method
            % -------------
            %             if(config.useExxEyy==1 || config.useRotation==1 )
            %                 D_out = obj.getDmatrixForExxEyy(Exx, Eyy); % Exx and Eyy DISTRIBUTION
            %                 D_out=D_out*topDensity^(config.penal); % TOPOLOGY;
            %                 if(config.useRotation==1)  % ROTATION
            %                     D_out=obj.rotateDmatrix(config,rotation, D_out)  ;
            %                 end
            %                 [K]=elK_elastic_v2(D_out);
            %             else % END NEW METHOD
            
            
            % -------------
            % Old Method method
            % -------------
            %E0 = effectiveElasticProperties(obj, material1Fraction, config); % GRADIENT MATERIAL
            E0 = obj.E_material1-obj.E_material2;
            E0 = E0*topDensity^(config.penal); % TOPOLOGY
            if(config.macro_meso_iteration>1)
                Dgiven =matProp.GetSavedDMatrix(count);
            else
                Dgiven = [];
            end
            [K, ~, ~, ~]=elK_elastic(E0,obj.v, obj.G,[],Dgiven);
            %             end% END OLD METHOD
        end
        
        %---------------
        % ROTATES THE D MATRIX
        %
        %  USE RADIANS,
        %--------------
        function Dout = rotateDmatrix(obj,config,theta, D_in)
            c = cos(theta);
            s = sin(theta);
            
            T = [c^2 s^2 2*s*c;
                s^2 c^2 -2*s*c;
                -s*c s*c c^2-s^2];
            R = [1 0 0;
                0  1 0;
                0  0 2];
            
            
            Dout = T^-1*D_in*R*T*R^-1;
        end
        
        
        %---------------
        % Get D matrix for a particular Exx and Eyy
        % -----------------------------------------------
        function D_out = getDmatrixForExxEyy(obj,Exx, Eyy,config)
            
            %             E0 = 1;%effectiveElasticProperties(obj, material1Fraction, config);
            vv = obj.v;
            %             d = orthD;
            E12 = (Exx+Eyy)/2;
            
            Q12 = vv*E12/(1-vv^2);
            Q66 =  0.5*(1-vv)*E12/(1-vv^2);
            
            if(config.useRinOrthMaterialModel==0)
                  r=1;
            else
              %  r=1-(Exx/obj.E_material1-Eyy/obj.E_material1)^2;
              r = min(Exx/Eyy, Eyy/Exx);
            end
       
            
            D_out = [ Exx/(1-vv^2)        Q12*r                    0                   ;
                      Q12*r                 Eyy/(1-vv^2)           0                   ;
                      0                   0                      r*Q66  ] ;
            
        end
        
        
        %---------------
        % Get D matrix for anisotropic material
        % -----------------------------------------------
        function D_out = getDmatrixForAnisotropic(obj,Exx, Eyy,E12, E33)
            
            %             E0 = 1;%effectiveElasticProperties(obj, material1Fraction, config);
            vv = obj.v;
%             %             d = orthD;
             E12_local = E12/(1-vv^2);   
%             
%             Q12 = vv*E12/(1-vv^2);
%             Q66 =  0.5*(1-vv)*E12/(1-vv^2);
            D_out = [ Exx/(1-vv^2)        E12_local            0                        ;
                      E12_local           Eyy/(1-vv^2)         0                        ;
                      0                   0                    0.5*(1-vv)*E33/(1-vv^2)] ;
            
        end
        
        %---------------
        % Get D matrix sensitivity (Effectie constitutive equation) for a particualar
        % w1 and orthotropic Distribution value, d
        % -----------------------------------------------
        %         function D_out = calculateEffConsMatrixWithGradAndOrthDistrbution_sensitivity(obj, material1Fraction, config,orthD)
        %
        %             E0 = effectiveElasticProperties(obj, material1Fraction, config);
        %             vv = obj.v;
        %             d = orthD;
        %
        %             D_out = [ E0/(1-vv^2)       vv*E0*0.5/(1-vv^2)      0                   ;
        %                 vv*E0*0.5/(1-vv^2)     E0*(-1)/(1-vv^2)     0                   ;
        %                 0                   0                     vv*E0*0.5/(1-vv^2)  ] ;
        %
        %         end
        
        %------------------------------------
        % Calculate Elastic mod
        %------------------------------------
        function e =  effectiveElasticProperties(obj, material1Fraction, config)
            
            % linear case
            if(config. elasticMaterialInterpMethod == 1)
                e = material1Fraction*obj.E_material1+(1-material1Fraction)*obj.E_material2;
                % HashinShtrikamAverage case
            elseif(config. elasticMaterialInterpMethod == 2)
                e = HashinShtrikamAverage(obj,obj.E_material1,obj.E_material2,material1Fraction);
            end
        end
        
        %------------------------------------
        % Calculate element stiffness matrix.
        %------------------------------------
        function [ke, kexpansionBar,B] = effectiveElasticKEmatrix(obj, material1Fraction, config,Dgiven)
            % Calculate E, then calculate the K element matrix
            E = effectiveElasticProperties(obj, material1Fraction, config);
            [ke, kexpansionBar, B, ~]=elK_elastic(E,obj.v, obj.G,[],Dgiven);
        end
        
        
        
        %------------------------------------
        % MESO
        %------------------------------------
        function [ke, kexpansionBar,B_total, F_meso] = effectiveElasticKEmatrix_meso(obj, material1Fraction, config,strain)
            % Calculate E, then calculate the K element matrix
            E = effectiveElasticProperties(obj, material1Fraction, config);
          
            [ke, kexpansionBar,B_total,F_meso]=elK_elastic(E,obj.v, obj.G,strain,[]);
        end
        
        % Calculate thermal expansion coefficient
        function e =  effectiveThermalExpansionCoefficient(obj,material1Fraction,config)
            e = material1Fraction*obj.alpha1+(1-material1Fraction)*obj.alpha2;
        end
        
        % ---------------------------------------
        %
        % Thermal properties caculations
        %
        % ---------------------------------------
        
        % ---------------------------------------
        % Calculate heat transfer coefficient
        % ---------------------------------------
        function e =  effectiveHeatProperties(obj,material1Fraction,config)
            % linear case
            if(config. heatMaterialInterpMethod == 1)
                e = material1Fraction*obj.K_material1+(1-material1Fraction)*obj.K_material2;
                
                % HashinShtrikamAverage case
            elseif(config. heatMaterialInterpMethod == 5)
                e = HashinShtrikamAverage(obj,obj.K_material1,obj.K_material2,material1Fraction);
            end
        end
        
        % ---------------------------------------
        % Calculate the heat conduction matrix
        % ---------------------------------------
        function  kheat = effectiveHeatKEmatrix(obj, material1Fraction, config)
            heatCoefficient = obj.effectiveHeatProperties(material1Fraction,config);
            [kheat]=elementK_heat(heatCoefficient);
        end
        
        % ---------------------------------------
        % Calculate Hashin-Shtrikam upper and lower bound and then average
        % them
        % ---------------------------------------
        function effective =  HashinShtrikamAverage(obj,mat1property, mat2property, volmat1)
            w = volmat1;
            E1 = mat1property;
            E2 = mat2property;
            
            lower = ((2+w)*E1+(1-w)*E2)*E2/(2*(1-w)*E1+(1+2*w)*E2);
            upper = (w*E1+(3-w)*E2)*E1/((3-2*w)*E1 +2*w*E2);
            effective = (lower+upper)/2.0;
            
        end
    end
end