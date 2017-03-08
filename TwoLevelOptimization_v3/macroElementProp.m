classdef macroElementProp
    
    properties
        % --------------------------------------
        % %% Givens
        % --------------------------------------------
        disp; % Displacement of the element, the 1x8 matrix giving the xy displacement at each node.
        elementNumber; % in the macro problem
        elementNodes; % in the macro problem
        
        yPos; % for x and w matrix arrays, y position in the array
        xPos;
        
        % Macro design variables for this element
        densitySIMP; % SIMP density of this macro element
        
        % dont' need these (at least I don't think so!!). 
        material1Fraction; % the volume fraction of material 1 at this element
        Exx;
        Eyy;
        theta;
        
        % 
        psuedoStrain;
        
      
        
        
        % --------------------------------------
        % %% Calculated
        % --------------------------------------------
        targetDensity; % The target density of this element
        B; % Displacement strain matrix
        strain; % strain
        K; % Stiffness matrix
        
        
        D_subSys; % The sub system's D
        D_sys; % the D of the system
        % eventually these should converge. 
        
        
       
        
    end
end