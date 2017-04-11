#!/bin/bash
#PBS -N test
#PBS -l select=1:ncpus=2,walltime=01:00:00
#PBS -j oe

source /etc/profile.d/modules.sh
module purge
module load gcc/5.3.0 matlab/2016a


cd $PBS_O_WORKDIR



dos2unix *.m


# compile the matlab code
mcc -R -nodisplay -C  -m  FitRuntimeErrorExample.m  

./FitRuntimeErrorExampletest







