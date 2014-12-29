import subprocess
for x in range(0, 11):
    pbsname = "job%imatlab.pbs"  % x
  
    print pbsname
    name = "job%i" % x
   
    w1 = x/10.
    w2 = 1-w1
    args = "  %f %f %i" % (w1, w2, x)
                    
                    
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
    subprocess.call(["qsub",pbsname])

            # write several .pbs scripts,
            # where each one has a slightly different name
            # and where each on sends differnet args to the matblab code. 
        
