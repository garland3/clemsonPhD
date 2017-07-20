# OPTIMAL DESIGN OF GRADIENT MATERIALS and BI-LEVEL OPTIMIZATION OF TOPOLOGY USING TARGETS (BOTT)

This code is part of my Ph.D. research effort at  Clemson University in the ME department. The objective is best take advantage of additive manufacturing (AM) (often called 3D printing) by generating optimal designs. These optimal designs have complex topologies (and shapes) and have gradient material distributions. Gradient materials have varying material properties within a single part. By varying the mateiral composition or mesostructure at each location, these gradient materials can be manfauctured. Gradient mateirals allow the engineer to customize each region of the part to have an optimal design. The output of this research effort is algorithms which generate optimal topologies and have gradeint materials. 

This research effort had several sections. 
* An open loop control for the Big Builder dual extruder 3D printer was developed to allow printing of gradient materials. The code reads some gcode, some material gradient distribution is specified in the python program, and the code modifes the existing gcode to have the material gradient specified in the python program. (clemsonPhD\OpenLoopGcodeController\postProcessAddGradient_v2.py)
* Optimal design of gradeint materails and topoly. This section of the research finds the optimal design of an object when given loads and boundary conditions. 
** The gradient materail can be isotropic material
** The gradient material can be varying mesostructures (within the BOTT algorithm)

## Getting Started

*Download the repository. 
*The Feedback loop is in \OpenLoopGcodeController
*The main topology and gradient material optimization code is in 'TwoLevelOptimization_v3'



### Prerequisites

What things you need to install the software and how to install them

* matlab  matlab/2016a
* Some of the code is in C for a cluster computer. These versions were used for compiling in Linux gcc/5.3.0 openmpi/1.10.3
* Python 2.7 for the feedback loop. 

```
Give examples
```

### Installing

For running the topology and gradient material optimization code on a desktop windows computer, you will need to 
* Open everything in Matlab
* Open the Configuration file and set the configurations to how you like. 
* Open the FE_elastic and temperatureFEA_V3 and setup the boundary conditions that you want. The elastic FEA has several boundary conditions possible. Specify all the loading conditions you want in the Configuration file by updating the 'loadingCase' array. 
* If running on the cluster computer, then you will need to upload the 'TwoLevelOptimization_v3\MPI_C_program'
** Within the C program, modify the #define to the correct mode . 4 modes exist for 1. Bott, 2. Validating the meosdesign (ie, run a lot of examples and compare inputs and outputs), 3. Run Coelho's algorithm 4. Run the psuedo strain and target density experiment
** To submit a job to the cluster, you can upload job.mpiMatlab.pbs. Run Dos2unix on the file, and the submit it (qsub job.mpiMatlab.pbs). 
** You can change the queue that the job will enter based on the resource request. Several options are given at the bottom of the job.mpiMatlab.pbs file. Copy and paste a resource request into the 3rd line of the file, and also change the 'mpirun -np 256 ' number. 
* The main.m file has several prebuilt input args for the combinedTopologyOptimization.m file. The combinedTopologyOptimization.m file is the true main file. Specify a mode and let it run. 
* The C program will call the compiled combinedTopologyOptimization.m file in different modes in order to run the complete BOTT algorithm. 
* On a desktop computer it might be hard to run the whole BOTT algorithm, but you can easily run the macro and individual mesostructure optimization algorithms by specifiying the mode and/or the element number. 
* Monitor progress on the cluster computer with 
** qstat -u username
** qpeek jobid | tail -n 100




## Acknowledgments

Let me know if you need help.
If you use the code in research please cite me!

