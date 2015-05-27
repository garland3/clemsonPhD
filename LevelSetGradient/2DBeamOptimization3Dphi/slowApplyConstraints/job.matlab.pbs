#!/bin/bash
#PBS -N slowApply
#PBS -l select=1:ncpus=16:mem=16gb,walltime=8:00:00
#PBS -j oe

source /etc/profile.d/modules.sh
module purge
module add matlab/2014a

cd $PBS_O_WORKDIR

rm *.csv

mcc -R -nodisplay  -m  GAOptimizationMain_v2_Single.m FEALevelSetWrapperGA_v2.m FEALevelSet_2D_v8.m  objfun2_v2.m outPutFunction.m freezeColors.m constr2list1.m constr2list2.m constr2list3.m constr2list4.m constr2list5.m constr2list6.m constr2list7.m

./GAOptimizationMain_v2_Single 1

