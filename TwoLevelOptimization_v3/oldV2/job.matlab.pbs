#!/bin/bash
#PBS -N {NAME}
#PBS -l select=1:ncpus=10:mem=10gb,walltime=2:00:00
#PBS -j oe

source /etc/profile.d/modules.sh
module purge
module add matlab/2015a
# module add python/2.7.6 

cd $PBS_O_WORKDIR

./combinedTopologyOptimization {ARGS_SCRPT_CALL}




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







