    classdef MaterialProperties
    
    
    properties 
%         E_material1 = 4; % The elastic mod of material 1
%         E_material2 = 2; % The elastic mod of material 2
%         K_material1 = 2; % heat conduction of material 1
%         K_material2 = 4; % heat conduction of material 2    

        % ceramic Al2O3, Aluminium oxide
%         E_material2 = 380e9; % The elastic mod of material 1, 
%         K_material2 = 8.4; % heat conduction of material 1


        % TiC, Titanium carbide
%         E_material1 = 450e9; % The elastic mod of material 2
%         K_material1 = 24.28; % heat conduction of material 2    


        % http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thrcn.html
        % http://hyperphysics.phy-astr.gsu.edu/hbase/permot3.html
        % material 1, copper
         E_material1 = 10e9/1e6; %N/mm^2 The elastic mod of material 1
         K_material1 = 385/1000; % W/(Kelvin*mm) heat conduction of material 1
         
         E_material2 = 200e9/1e6; % N/mm^2 The elastic mod of material 2         
         K_material2 = 50.2/1000; % heat conduction of material 2    


        
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
    
        
        % Calculate Elastic mod
        function e =  effectiveElasticProperties(obj, material1Fraction, settings)
            
            
            if (settings.elasticMaterialInterpMethod ==1 )
                % -------------
                % Simple linear interpolation 
                %--------------
                e = material1Fraction*obj.E_material1+(1-material1Fraction)*obj.E_material2;
                
            elseif(settings.elasticMaterialInterpMethod ==2)
                % -------------
                %  Hashin–Shtrikam law
                %--------------
                w = material1Fraction ;
                E1 = obj.E_material1;
                E2 = obj.E_material2;
                
                
                lower = ((2+w)*E1+(1-w)*E2/(2*(1-w)*E1+(1+2*w)*E2))*E2;
                upper = ((w*E1)+(3-w)*E2)/((3-2*w)*E1+2*w*E2)*E1;
                e = (lower+upper)/2;
                %e = material1Fraction*obj.E_material1+(1-material1Fraction)*obj.E_material2;
            elseif(settings.elasticMaterialInterpMethod ==4)
                
                % page 3 of  Development of an advanced ceramic tool material—functionally
                
                alpha = 1/3*(1+obj.v)/(1-obj.v);
                Beta =  2/15*(4-5*obj.v)/(1-obj.v);
                c1  = material1Fraction ; % amount of ceramic in matrix
                c2 = 1-c1;
                
                E1 = obj.E_material1;
                E2 = obj.E_material2;
                
                ratio = 1+ c2*(E2-E1)/(E1+alpha*(1-c2)*(E2-E1));
                e = ratio*E1;
                
            else
                
            end
        end
        
        function ke = effectiveElasticKEmatrix(obj, material1Fraction,settings)
            % Calculate E, then calculate the K element matrix
            E = effectiveElasticProperties(obj, material1Fraction, settings);
            [ke]=elK_elastic(E,obj.v, obj.G);
            
        end
        
        
        
        % Calculate heat transfer coefficient
        function k =  effectiveHeatProperties(obj,material1Fraction,settings)
            
            if(settings.heatMaterialInterpMethod ==1)
            	k = material1Fraction*obj.K_material1+(1-material1Fraction)*obj.K_material2;
            elseif(settings.heatMaterialInterpMethod ==4)
                
                 % page 3 of  Development of an advanced ceramic tool material—functionally                
             
                c1  = material1Fraction ; % amount of ceramic in matrix
                c2 = 1-c1;
                
                K1 = obj.K_material1;
                K2 = obj.K_material2;
                
                numerator = 1+2*c2*(1-K1/K2)/(2*K1/K2+1);
                denominator = 1-c2*(1-K1/K2)/(K1/K2+1);
                k = numerator/denominator;
                
            end
        end

        % Calculate the heat conduction matrix
        function  kheat = effectiveHeatKEmatrix(obj, material1Fraction, settings)
            heatCoefficient = obj.effectiveHeatProperties(material1Fraction, settings);
            [kheat]=elementK_heat(heatCoefficient);            
        end
        
    end
    
    
    
end