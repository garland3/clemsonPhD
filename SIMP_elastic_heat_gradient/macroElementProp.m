classdef macroElementProp
    
    properties
        % --------------------------------------
        % %% Settings
        % --------------------------------------------
        disp; % Displacement of the element, the 1x8 matrix giving the xy displacement at each node. 
        B; % Displacement strain matrix
        strain; % strain
        density; % density of this macro element
        K; % Stiffness matrix
        
        D_homog;
    end
end