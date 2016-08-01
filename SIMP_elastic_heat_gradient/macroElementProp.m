classdef macroElementProp
    
    properties
        % --------------------------------------
        % %% Givens
        % --------------------------------------------
        disp; % Displacement of the element, the 1x8 matrix giving the xy displacement at each node.
        density; % density of this macro element
        material1Fraction; % the volume fraction of material 1 at this element
        elementNumber; % in the macro problem
        elementNodes; % in the macro problem
        
        yPosition; % for x and w matrix arrays, y position in the array
        xPosition;
        
        
        % --------------------------------------
        % %% Calculated
        % --------------------------------------------
        B; % Displacement strain matrix
        strain; % strain
        K; % Stiffness matrix
        
        
        D_homog;
    end
end