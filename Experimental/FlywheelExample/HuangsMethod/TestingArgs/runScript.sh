# Run qsub -I before running this script. 
dos2unix *.*
module load matlab/2013a



mcc -v -N -R -nodisplay  -R  -nojvm -m testArgs.m

./run_testArgs.sh $MATLAB 11