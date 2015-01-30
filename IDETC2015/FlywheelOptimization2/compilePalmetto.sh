 # load an interactive session
 # qsub -I
 
 # load the matlab
module load matlab/2013a

 # make sure everything is ok for linux
    
dos2unix *.*
 
 # compile
 # -N don't load all the toolbox components. This makes compiling take for ever
 # -p add the path to the optimization toolbox module so we can use fmincon
 # -nodisplay. There is no display
 # -R singleCompThread : Only take one thread or else it tries to take everyone's sources on the partiular node
 # Compile optimize.m and include all the other files. 
 
 # single thread
 #  mcc -v -N -p /software/matlab/R2013a/toolbox/optim -R -nodisplay -R -singleCompThread -R -nojvm -m optimize.m bezierInter.m calculateDensity.m CalculateEnergy.m CalculateMass.m CalculateMaxStress.m  convertToT.m DerivBezierInter.m energyDifferential.m InputsToHpointsAndVolFractPoints.m massDifferential.m  nonlinearConstraints.m Objective.m
   
   # multithread
mcc -v -N -p /software/matlab/R2013a/toolbox/optim -R -nodisplay  -R -nojvm -m optimize.m bezierInter.m calculateDensity.m CalculateEnergy.m CalculateMass.m CalculateMaxStress.m  convertToT.m DerivBezierInter.m energyDifferential.m InputsToHpointsAndVolFractPoints.m massDifferential.m  nonlinearConstraints.m Objective.m
  
  #mkdir /newscratch/apg/flywheel3
  #cp ./* /newscratch/apg/flywheel3
  
#  ./run_my_optimize.sh $MATLAB 