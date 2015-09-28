import os
import shutil
import stat
import subprocess

# precompile the matlab code

# dos2unix *.*
# module add python/2.7.6 
# module add matlab/2015a
# mcc -R -nodisplay  -m  combinedTopologyOptimization.m Configuration.m DesignVars.m  elementK_heat.m elK_elastic.m FE_elasticV2.m  MaterialProperties.m plotResults.m temperatureFEA_V3.m  


# bump up the resolution on the analysis code. 40 x20 at least.  (80 x40 would be better
# make a folder to store the .csv files
# clear previous data
# compile the matlab code
# loop from 1 1o 10
# - make a new folder with the iterationNumber to export the .csv files to
# - run the matlab code, accept in put args for the weights of the objectives. 
# -- change the csv output location based off an input arg
# 
# Note, request 12 cpus, so that each run can have its own cpu and we can run in parallel. 

# For testing, use interactive connection to palmetto
#  qsub -I -l select=1:ncpus=2,walltime=3:00:00



pathHome = os.path.dirname(os.path.realpath(__file__))

programName = 'main'
programFullPath = os.path.join(pathHome,programName)
print programFullPath




for x in range(0, 11):
	directoryName = "out"+str(x)
	directoryNew = os.path.join(pathHome,directoryName)
	
	# remove the previous data
	if  os.path.exists(directoryNew):
		shutil.rmtree(directoryNew) 
	
	# make a new directory
	if not os.path.exists(directoryNew):
		os.makedirs(directoryNew)
		
	
	name = "jobP%i" % x
	pbsname = "jobweight%imatlab.pbs"  % x
   
    # Make different weightings of the objective functions
	w1 = x/10.
	#w2 = 1-w1
	args = "  1 %f %i" % (w1,  x)
                
                    
	with open(pbsname, "wt") as fout:
		with open("job.matlab.pbs", "rt") as fin:
			for line in fin:             
				if (line.find("{NAME}")  != -1):                    
					fout.write(line.replace('{NAME}', name))
				elif (line.find("{ARGS_SCRPT_CALL}")  != -1):                    
					fout.write(line.replace('{ARGS_SCRPT_CALL}', args))
				else:
					fout.write(line)

	cmd = "qsub "+ pbsname
	print cmd
	subprocess.call(["qsub",pbsname])
		
	
	#newProgramLocation = os.path.join(directoryNew,programName)
	#shutil.copyfile(programFullPath, newProgramLocation)
	#os.chmod(newProgramLocation,stat.S_IEXEC)
		
	
	
# import subprocess

# write several .pbs scripts,
# where each one has a slightly different name
# and where each on sends different args to the matblab code. 
# Submit the .pbs file to job queue on the cluster using 'qsub'

# for x in range(0, 11):
    # pbsname = "job%imatlab.pbs"  % x
  
    # print pbsname
    # name = "job%i" % x
   
    # Make different weightings of the objective functions
    # w1 = x/10.
    # w2 = 1-w1
    # args = "  %f %f %i" % (w1, w2, x)
                    
                    
    # with open(pbsname, "wt") as fout:
         # with open("job.matlab.pbs", "rt") as fin:
             # for line in fin:             
                # if (line.find("{NAME}")  != -1):                    
                    # fout.write(line.replace('{NAME}', name))
                # elif (line.find("{ARGS_SCRPT_CALL}")  != -1):                    
                    # fout.write(line.replace('{ARGS_SCRPT_CALL}', args))
                # else:
                 # fout.write(line)
    
    # cmd = "qsub "+ pbsname
    # subprocess.call(["qsub",pbsname])