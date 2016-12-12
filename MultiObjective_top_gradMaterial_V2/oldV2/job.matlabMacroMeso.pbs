#!/bin/bash
#PBS -N macroMeso1
#PBS -l select=1:ncpus=9:mem=12gb,walltime=10:00:00
#PBS -j oe

source /etc/profile.d/modules.sh
module purge
module add matlab/2015a
# module add python/2.7.6 

cd $PBS_O_WORKDIR

dos2unix *.*
rm jobP*
rm jobweight*
module add matlab/2015a
mcc -R -nodisplay  -m  main.m CalculateCheckedElements.m check.m combinedTopologyOptimization.m Configuration.m DesignVars.m elementK_heat.m elK_elastic.m FE_elasticV2.m FindAvergeChangeOfLastValue.m freezeColors.m GenerateDesignVarsForMesoProblem.m GetMacroElementPropertiesFromCSV.m GetMesoUnitCellDesignFromCSV.m Homgenization.m HomgenizationV2.m macroElementProp.m MaterialProperties.m mesoAddAdjcentCellDataObject.m  MesoDesignWrapper.m MesoStructureDesignV2.m OC.m parfor_progress.m plotResults.m  plotSingleMatCompare.m plotStrainField.m plotStructGradFromCSVfile.m ReadXMacroFromCSV.m SaveMacroProblemStateToCSV.m SaveMesoUnitCellDesignToCSV.m SingleMesoStuctureWrappe.m stlwrite.m temperatureFEA_V3.m TestMesoDesign.m TileMesoStructure.m unfreezeColors.m

./main




# use 10 cpus and 10 gb, becasue palmetto keeps killing my jobs becasue of not enough cpus
# rm *.csv
# bump up the resolution on the analysis code. 40 x20 at least.  (80 x40 would be better
# make a folder to store the .csv files
# clear previous data
# compile the matlab code
# loop from 1 1o 10
# - make a new folder with the iterationNumber to export the .csv files to
# - run the matlab code, accept in put args for the weights of the objectives. 
# 
# Note, request 12 cpus, so that each run can have its own cpu and we can run in parallel. 
#mcc -R -nodisplay  -m  combinedTopologyOptimization.m Configuration.m DesignVars.m elastic_top.m  elementK_heat.m elK_elastic.m FE_elasticV2.m  MaterialProperties.m plotResults.m temperatureFEA_V3.m  

#combinedTopologyOptimization('1','0.5','1')







