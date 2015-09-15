classdef MaterialProperties
    
    
    properties 
        E_material1 = 4; % The elastic mod of material 1
        E_material2 = 2; % The elastic mod of material 2
        K_material1 = 2; % heat conduction of material 1
        K_material2 = 4; % heat conduction of material 2    
        
        v = 0.3; % Piossons ratio
        G ; % = E_material1/(2*(1+v));
    end
    
    methods
        % Constructor
        function obj = MaterialProperties
             obj.G = obj.E_material1/(2*(1+obj.v));
        end
    
        % Calculate Elastic mod
        function e =  effectiveElasticProperties(material1Fraction)
            e = material1Fraction*obj.E_material1+(1-material1Fraction)*obj.E_material2;
        end
        
        % Calculate heat transfer coefficient
        function e =  effectiveHeatProperties(material1Fraction)
             e = material1Fraction*obj.K_material1+(1-material1Fraction)*obj.K_material2;
        end
        
    end
    
    
    
end