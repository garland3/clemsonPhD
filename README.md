# OPTIMAL DESIGN OF GRADIENT MATERIALS and BI-LEVEL OPTIMIZATION OF TOPOLOGY USING TARGETS (BOTT)

This code is part of my Ph.D. research effort at  Clemson University in the ME department. The objective is best take advantage of additive manufacturing (AM) (often called 3D printing) by generating optimal designs. These optimal designs have complex topologies (and shapes) and have gradient material distributions. Gradient materials have varying material properties within a single part. By varying the mateiral composition or mesostructure at each location, these gradient materials can be manfauctured. Gradient mateirals allow the engineer to customize each region of the part to have an optimal design. The output of this research effort is algorithms which generate optimal topologies and have gradeint materials. 

This research effort had several sections. 
1. An open loop control for the Big Builder dual extruder 3D printer was developed to allow printing of gradient materials. The code reads some gcode, some material gradient distribution is specified in the python program, and the code modifes the existing gcode to have the material gradient specified in the python program. (clemsonPhD\OpenLoopGcodeController\postProcessAddGradient_v2.py)
2. Optimal design of gradeint materails and topoly. This section of the research finds the optimal design of an object when given loads and boundary conditions. 
    1. The gradient materail can be isotropic material
        1. First, a 3D spline is used to model topology and gradient materrial. All locations of the surface above the zero-th isocurve have material. From 0 to 1 determines the composition from material 1 transitioning to material 2 where material 2 is more stiff than material 1. 
        2. A discritized model is used. SIMP topology optimization finds the best topology. A gradient material optimizer finds the best material distribution within the object. The gradient material optimizer uses a Augmented Lagrangian Multiplier appoarach or an optimal criteria approach in order to meet the ratio of the two mateirals constraint. 
        3. Using the ideas in method 2, the research explores considering conflicting objectives. The objects are to maximize stiffness and minimize steady-state heat transfer of the object. A weighted objcective controls the influence of the two objectives. 
    2. The gradient material can be varying mesostructures (within the BOTT algorithm). 
        1. First a macroscopic topolgoy and gradient material distribution are generated. SIMP finds the optimal topology. The gradeint material is modeled as an orthortopic material where the Exx and Ey are design variables in addition to the rotation. 
        2. A mesostructure design technique generates mesostructures at each location that are equal to the macro-level's material properties at each location. D_sys=D_meso^h. 
        3. Coordination of the two levels is facilitated through Analytical Target Cascacading (ATC) which adds constraints to the macro problem so that it must agree with the outputs of the meso structure design problems. 
        4. Additional code is added to help show the progress of the optimizaiton and generate the final complete structure. 

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
1. Open everything in Matlab
2. Open the Configuration file and set the configurations to how you like. 
3 Open the FE_elastic and temperatureFEA_V3 and setup the boundary conditions that you want. The elastic FEA has several boundary conditions possible. Specify all the loading conditions you want in the Configuration file by updating the 'loadingCase' array. 
4. If running on the cluster computer, then you will need to upload the 'TwoLevelOptimization_v3\MPI_C_program'
    1. Within the C program, modify the #define to the correct mode . 4 modes exist for 1. Bott, 2. Validating the meosdesign (ie, run a lot of examples and compare inputs and outputs), 3. Run Coelho's algorithm 4. Run the psuedo strain and target density experiment
    2. To submit a job to the cluster, you can upload job.mpiMatlab.pbs. Run Dos2unix on the file, and the submit it (qsub job.mpiMatlab.pbs). 
    3. You can change the queue that the job will enter based on the resource request. Several options are given at the bottom of the job.mpiMatlab.pbs file. Copy and paste a resource request into the 3rd line of the file, and also change the 'mpirun -np 256 ' number. 
5. The main.m file has several prebuilt input args for the combinedTopologyOptimization.m file. The combinedTopologyOptimization.m file is the true main file. Specify a mode and let it run. 
6.  The C program will call the compiled combinedTopologyOptimization.m file in different modes in order to run the complete BOTT algorithm. 
7. On a desktop computer it might be hard to run the whole BOTT algorithm, but you can easily run the macro and individual mesostructure optimization algorithms by specifiying the mode and/or the element number. 
8.  Monitor progress on the cluster computer with 
    1.  qstat -u username
    2.    qpeek jobid | tail -n 100
    
 

### Examples

#### Run Just topology optimization
Configure Configuration.m so that v1=Target density . In MaterialProperties.m specify the elastic modulus of the base material. IN temperatureFEA_V3.m and FE_elasticV2.m setup the boundary conditions and loading conditions. This must be done manually using the node numbers. 

The second arg to the  combinedTopologyOptimization specifies the weight of the dual objective weighing function.

The 4th argument is the mode number. Mode 1 is only topology optimization. Mode 200 plots the output design metrics.  

Run the code
```
combinedTopologyOptimization('1', '0.5', '1','1', 'na');
combinedTopologyOptimization('1', '0.5', '1','200', 'na');
```

#### Run Topology optimization and Gradient material optimization
1. Specify the target density of each material by setting v1=density of mateiral 1 and v2= density of material 2 in the Configuration.m file. 
2. Set useExxEyy=0 and useRotation=0 in the Configuation file. 

Run the code

```
combinedTopologyOptimization('1', '0.5', '1','1', 'na');
combinedTopologyOptimization('1', '0.5', '1','200', 'na');
```


#### Run Topology optimization and Gradeint Orthotropic material (BOTT macro optimization)
1. Specify the target total topology density be setting v1=target density. v2=0
2. Set useExxEyy=1 and useRotation=1 in the Configuation file. 
3. Specify to use a target average mesostructure density or a target avg elastic modulus useTargetMesoDensity in the Configuation.m file. 
    1. Target avg E is targetAvgExxEyy
    2. Target density is targetExxEyyDensity
4. Specify the minimum summed elastic modulus at each location (Exx+Eyy>Emin), minEallowed=

Run the code

```
 combinedTopologyOptimization('1', '1', '1','60', 'na');
 combinedTopologyOptimization('1', '1', '1','200', 'na');
 ```
 
 
 #### Run Meosstructure Optimization
 1. Make sure you have design targets for the mesostructure optimization by running the macro optimization first. 
 2. Select the method of mesostructure design that you want. 
 
 ```
UseLookUpTableForPsuedoStrain=1; %0 = feedback loop, 1 = use look up.
mesoVolumeUpdateMethod=1; % 1 = average, 2 = Target the larger
lookupSearchScheme=2; % 2 = search table and scale eta, 4 = linear interpolation with particle swarm
 ```
 
 Then run the code for all the mesostructures. N is the number of macro elements
 
 ```
for i =1:N
    combinedTopologyOptimization('1', '1', '1','100', int2str(i)); 
end
 ```
 
 #### Run Pseudo strain and density target Experiment
 1. Turn on the test configuration. strainAndTargetTest=1;
 2. Specify the length of the vector of combinations (step=1/targetTestVectorLen)

 The code generates the targets, then runs the mesostructure optimization. 
 
 ```
combinedTopologyOptimization('1', '1', '1','113', 'na'); % genrate psuedo strain and density targets
for i =1:N
    combinedTopologyOptimization('1', '1', '1','100', int2str(i)); 
end
combinedTopologyOptimization('1', '1', '1','114', 'na'); % save psuedo strain and density targets results in data files
combinedTopologyOptimization('1', '1', '1','115', 'na'); % plot psuedo strain and density targets results
```

In the C program, you specify mode 4 to run the same code, but the mesostructure designs are run in parallel. 

#### Run the mesostructure validation test
1. Turn on the validation mode validationModeOn=1
2. Specify the number of problems to run validationGridSizeNelx^3 is the number of problems. 

The code systematially changes Exx, Eyy, and theta targets for the mesostructure design problem. 

```
combinedTopologyOptimization('1', '1', '1','111', 'na');
for i =1:N
    combinedTopologyOptimization('1', '1', '1','100', int2str(i)); 
end
combinedTopologyOptimization('1', '1', '1','112', 'na');
```


In the C program, you specify mode 2





## Acknowledgments

Let me know if you need help.
If you use the code in research please cite me!

