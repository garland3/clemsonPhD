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
        
        E_material1 = 100000; %4 N/mm^2 The elastic mod of material 1
        E_material2 =  50000;%E_material1/2; %4 The elastic mod of material 2
        
        K_material1 = 0.02; %  W/ (mm*K)heat conduction of material 1
        K_material2 = 0.04; % heat conduction of material 2
        
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
        
        % -------------------------------------------------------
        % w is vol fraction of material 1, x is density, at this meso
        % element or design var.
        % -------------------------------------------------------
        function [density]=CalculateDensityTargetforMeso(obj,w,x,settings)
            %             density=0.5+(w*x^settings.penal)/2;
            
            %             density = home much greater than nothing material 2 is (scaled to 0 to 1) + (vol fraction * density^p)/how much is between v2 and 1
            offset = obj.E_material2/obj.E_material1;
            density = offset+(w)*(1-offset);
            %   density=density-0.2;
        end
        
        % ------------------------------------------------
        %
        %             READS THE D MATRIXES AND SAVES THEM TO AN ARRAY FOR
        %             FUTURE USE.
        %
        % ------------------------------------------------
        function obj =  ReadConstitutiveMatrixesFromFiles(obj,  settings)
            
            folderNum = settings.iterationNum;
            oldIteration = settings.macro_meso_iteration-1; % minus 1, because we want to get the previous iterationdesign.
            
            ne = settings.nelx*settings.nely;
            obj.SavedDmatrix = zeros(ne,9);
            
            % if not a single meso design
%             if(settings.singleMesoDesign ~=1)
                for e = 1:ne
                    % Read the D_h
                    outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,oldIteration,e);
                    if exist(outname, 'file') == 2
                        D_h= csvread(outname);
                    else
                        material1Fraction=1;
                        D_h = obj.calculateEffectiveConstitutiveEquation(material1Fraction, settings,[]);
                    end
                    
                    
                    % Read the D_old
                    outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,oldIteration,e);
                    if exist(outname, 'file') == 2
                        D_old= csvread(outname);
                    else
                        material1Fraction=1;
                        D_old = obj.calculateEffectiveConstitutiveEquation(material1Fraction, settings,[]);
                    end
                    
                    % --------------------
                    % Average the 2 Ds (D_h and D_old
                    % --------------------
                    D = (D_h+D_old)/2;
                    
                    
                    D_flat = reshape(D,1,9);
                    obj.SavedDmatrix(e,:)=D_flat;
                end
%             else
%                 % ---------------------------------------
%                 % if  a single meso design
%                 % Read the D_h, get the one and only
%                 % ---------------------------------------
%                 e = 1
%                 outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,oldIteration,e);
%                 if exist(outname, 'file') == 2
%                     D_h= csvread(outname);
%                 else
%                     material1Fraction=1;
%                     D_h = obj.calculateEffectiveConstitutiveEquation(material1Fraction, settings,[]);
%                 end
%                 
%                 
%                 for e = 1:ne
%                     D_flat = reshape(D_h,1,9);
%                     obj.SavedDmatrix(e,:)=D_flat;
%                 end
%                 
%             end
        end
        
        function [D] = GetSavedDMatrix(obj,e)
            D_flat =  obj.SavedDmatrix(e,:);
            D = reshape(D_flat,3,3);
        end
        
        
        % ---------------------------------------
        %
        % Elastic caculations
        %
        % ---------------------------------------
        
        %---------------
        % Get D matrix (Effectie constitutive equation) for a particualar
        % w1
        % -----------------------------------------------
        function D = calculateEffectiveConstitutiveEquation(obj, material1Fraction, settings,Dgiven)
            if(isempty(Dgiven))
                E = effectiveElasticProperties(obj, material1Fraction, settings);
                D = [ 1.1 obj.v 0;
                    obj.v 1 0;
                    0 0 1/2*(1-obj.v)]*E/(1-obj.v^2);
            else
                D = Dgiven;
            end
        end
        
        %------------------------------------
        % Calculate Elastic mod
        %------------------------------------
        function e =  effectiveElasticProperties(obj, material1Fraction, settings)
            
            % linear case
            if(settings. elasticMaterialInterpMethod == 1)
                e = material1Fraction*obj.E_material1+(1-material1Fraction)*obj.E_material2;
                
                % HashinShtrikamAverage case
            elseif(settings. elasticMaterialInterpMethod == 2)
                e = HashinShtrikamAverage(obj,obj.E_material1,obj.E_material2,material1Fraction);
            end
        end
        
        %------------------------------------
        % Calculate element stiffness matrix.
        %------------------------------------
        function [ke, kexpansionBar,B] = effectiveElasticKEmatrix(obj, material1Fraction, settings,Dgiven)
            % Calculate E, then calculate the K element matrix
            E = effectiveElasticProperties(obj, material1Fraction, settings);
            [ke, kexpansionBar, B, ~]=elK_elastic(E,obj.v, obj.G,[],Dgiven);
        end
        
        %------------------------------------
        % MESO
        %------------------------------------
        function [ke, kexpansionBar,B_total, F_meso] = effectiveElasticKEmatrix_meso(obj, material1Fraction, settings,strain)
            % Calculate E, then calculate the K element matrix
            E = effectiveElasticProperties(obj, material1Fraction, settings);
            [ke, kexpansionBar,B_total,F_meso]=elK_elastic(E,obj.v, obj.G,strain,[]);
        end
        
        % Calculate thermal expansion coefficient
        function e =  effectiveThermalExpansionCoefficient(obj,material1Fraction,settings)
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
        function e =  effectiveHeatProperties(obj,material1Fraction,settings)
            % linear case
            if(settings. heatMaterialInterpMethod == 1)
                e = material1Fraction*obj.K_material1+(1-material1Fraction)*obj.K_material2;
                
                % HashinShtrikamAverage case
            elseif(settings. heatMaterialInterpMethod == 5)
                e = HashinShtrikamAverage(obj,obj.K_material1,obj.K_material2,material1Fraction);
            end
        end
        
        % ---------------------------------------
        % Calculate the heat conduction matrix
        % ---------------------------------------
        function  kheat = effectiveHeatKEmatrix(obj, material1Fraction, settings)
            heatCoefficient = obj.effectiveHeatProperties(material1Fraction,settings);
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