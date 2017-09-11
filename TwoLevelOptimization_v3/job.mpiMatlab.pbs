#!/bin/bash
#PBS -N mpi1
#PBS -l select=16:ncpus=8:mpiprocs=4:mem=8gb:interconnect=10g,walltime=3:00:00
#PBS -l place=scatter
#PBS -j oe

source /etc/profile.d/modules.sh
module purge
module load gcc/5.3.0 openmpi/1.10.3 matlab/2016a
# module add python/2.7.6 

cd $PBS_O_WORKDIR

rm jobP*
rm jobweight*
rm *.csv
rm *.stl
rm *.png
rm mpi*
#rm ./out0/*.csv
rm ./bank
rm ./combinedTopologyOptimization

dos2unix *.m
dos2unix ./MPI_C_program/bankElementJobs.c


# compile the c code
mpicc ./MPI_C_program/bankElementJobs.c -o ./bank
echo 'Finished compiling MPI program'

# compile the matlab code
mcc -R -nodisplay -C  -m  combinedTopologyOptimization.m    Configuration.m DesignVars.m elementK_heat.m elK_elastic.m FE_elasticV2.m FindAvergeChangeOfLastValue.m freezeColors.m GenerateDesignVarsForMesoProblem.m GetMacroElementPropertiesFromCSV.m GetMesoUnitCellDesignFromCSV.m Homgenization.m macroElementProp.m MaterialProperties.m mesoAddAdjcentCellDataObject.m  MesoDesignWrapper.m MesoStructureDesignV2.m OC.m Optimizer.m parfor_progress.m plotResults.m  SaveMacroProblemStateToCSV.m SaveMesoUnitCellDesignToCSV.m   temperatureFEA_V3.m TestMesoDesign.m TileMesoStructure.m unfreezeColors.m  VideoManager.m GenerateMesoValidationTargets.m annOutput_matUpdateV2.m annOutput_matUpdateV1.m GeneratePsuedoStrainsAndDensityTargets.m annTest.m plotMode90Data.m annOutput_RandomMesoInitialLookUpTable.m EvalutePseudoStrainAndDensityForFit.m InterpolatePseudoStrainsAndDensity.m annOutput_lookupTable_withFmincon.m CheckForConerElements.m CheckRemoveOrphanedSegments.m GenerateCompleteStructureV2Improved.m ShowOptimizerProgress.m TestForTermaination.m ConvertCSVToSTL.m csvToStlConfig.m plotAnElement.m

# weird that sometimes this folder causes problems
#rm -R combinedTopologyOptimization_mcr

echo 'Finished compiling Matlab'
mpirun -np 64 --mca mpi_cuda_support 0 ./bank

echo 'Finished mpirun'
#./bank

#////   

# mpirun -np 460 ./bank





# Single
# ///#PBS -l select=3:ncpus=8:mpiprocs=4:mem=8gb:interconnect=10g,walltime=3:00:00   # 2 GB ram per processor and 4 cpus per 2 mpi processes (mpirun 12)

# Tiny Queue
# ///#PBS -l select=16:ncpus=8:mpiprocs=4:mem=8gb:interconnect=10g,walltime=3:00:00   # 2 GB ram per processor and 4 cpus per 2 mpi processes (mpirun 64)

# Small Queue
#//#PBS -l select=64:ncpus=8:mpiprocs=4:mem=8gb:interconnect=10g,walltime=3:00:00  # 2 GB ram per processor and 4 cpus per 2 mpi processes (mpirun 256)

#medium
#//#PBS -l select=256:ncpus=8:mpiprocs=4:mem=8gb:interconnect=10g,walltime=12:00:00  # 2 GB ram per processor and 4 cpus per 2 mpi processes (mpirun 1024)












