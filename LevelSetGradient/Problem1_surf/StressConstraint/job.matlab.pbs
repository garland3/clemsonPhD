#!/bin/bash
#PBS -N GA
#PBS -l select=1:ncpus=24:mem=24gb,walltime=8:00:00
#PBS -j oe

source /etc/profile.d/modules.sh
module purge
module add matlab/2014a

cd $PBS_O_WORKDIR

rm *.csv

mcc -R -nodisplay  -m  GAOptimizationMain_v2_Single.m FEALevelSetWrapperGA_v2.m FEALevelSet_2D_v8.m constr2_v2.m objfun2_v2.m outPutFunction.m

./GAOptimizationMain_v2

