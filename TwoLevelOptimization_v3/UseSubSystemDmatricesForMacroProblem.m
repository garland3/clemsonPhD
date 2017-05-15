function []=UseSubSystemDmatricesForMacroProblem(config,matProp)
DV = DesignVars(config);
  DV=DV.GetMacroStateVarsFromCSVFiles( config);
  
      matProp.ReadConstitutiveMatrixesFromFiles(  config)


end