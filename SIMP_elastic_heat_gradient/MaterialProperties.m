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
        
        E_material1 = 4; % N/mm^2 The elastic mod of material 1
        E_material2 = 2; % The elastic mod of material 2
        
        K_material1 = 0.02; %  W/ (mm*K)heat conduction of material 1
        K_material2 = 0.04; % heat conduction of material 2
        
        alpha1 =0.001; %thermal expansion coefficient for material 1
        alpha2 =0.001; % thermal expansion coefficient for material 2 
        
        v = 0.3; % Piossons ratio
        G ; % = E_material1/(2*(1+v));
        dKelastic; % derivative of K elastic matrix with respect to a material vol fraction change.
        dKheat;   % derivative of K heat matrix with respect to a material vol fraction change.
    end
    
    methods
        
        % Constructor
        function obj = MaterialProperties
            obj.G = obj.E_material1/(2*(1+obj.v));
            E = 1;
            obj.dKelastic = elK_elastic(E,obj.v, obj.G)*(obj.E_material1-obj.E_material2);
            
            heatCoefficient = 1;
            obj.dKheat =  elementK_heat(heatCoefficient)*(obj.K_material1-obj.K_material2);
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
        function D = calculateEffectiveConstitutiveEquation(obj, material1Fraction, settings)
            E = effectiveElasticProperties(obj, material1Fraction, settings);
            D = [ 1 obj.v 0;
                obj.v 1 0;
                0 0 1/2*(1-obj.v)]*E/(1-obj.v^2);
        end
        
        %------------------------------------
        % Calculate Elastic mod
        %------------------------------------
        function e =  effectiveElasticProperties(obj, material1Fraction, settings)
            e = material1Fraction*obj.E_material1+(1-material1Fraction)*obj.E_material2;
        end
        
        %------------------------------------
        % Calculate element stiffness matrix. 
        %------------------------------------
        function [ke, kexpansionBar] = effectiveElasticKEmatrix(obj, material1Fraction, settings)
            % Calculate E, then calculate the K element matrix
            E = effectiveElasticProperties(obj, material1Fraction, settings);
            [ke, kexpansionBar, ~]=elK_elastic(E,obj.v, obj.G);
        end
        
        % Calculate thermal expansion coefficient
        function e =  effectiveThermalExpansionCoefficient(obj,material1Fraction)
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
        function e =  effectiveHeatProperties(obj,material1Fraction)
            e = material1Fraction*obj.K_material1+(1-material1Fraction)*obj.K_material2;
        end
        
        % ---------------------------------------
        % Calculate the heat conduction matrix
        % ---------------------------------------
        function  kheat = effectiveHeatKEmatrix(obj, material1Fraction, settings)
            heatCoefficient = obj.effectiveHeatProperties(material1Fraction);
            [kheat]=elementK_heat(heatCoefficient);
        end
    end
end