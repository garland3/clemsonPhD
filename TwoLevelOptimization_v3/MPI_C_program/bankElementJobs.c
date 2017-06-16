#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include<stdlib.h>
#include <time.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
#define TERMINATION_TAG 123456789 
#define MASTER_RANK 0
#define SUB_MASTER_RANK 8
#define SUB_SUB_MASTER_RANK 2
#define NUM_MACRO_MESO_ITER 5 // 3,5,11
#define START_MACRO_MESO_ITER 1 // 1

//#define NELX_MACRO  15 // 30 // 1331 
//#define NELY_MACRO 15 //15 // 1

#define MODE 4 //1 =BOTT Multiscale optimization 2. Meso Validation 3.  multiscaleMethodCompareToCoelho 4. PseudoStrain and Density Target Test


int callMatlab(int mode, int macro_meso_iteration, int element);
int MultiscaleMasterNode( int macro_meso_iteration , int pool_size , int my_rank);
int MesoValidationMasterNode(int macro_meso_iteration,int pool_size,int my_rank );

 int DetermineNumberOfElements(int macro_meso_iteration,int folderNumber );



int main ( int argc, char **argv )
{    
   int pool_size, my_rank, destination;
   int i_am_the_master = FALSE; 
   int nelm; //Number of elements 
   MPI_Status status;
    FILE *log_file;
	int sender;
    int maxMacroMesoIterationCount;   
   
   // For each macro loop 
   // run macro level on master
   // broad cast (which will force everyone to wait)
   // Begin sending out jobs untill finished
   // run the combine together version of the matlab code. 
   // Get the Exx, eyy, theta values extracted from the D_h sub systems
   
   // Remove the log file before starting to append to it. 
   remove("./log.txt"); 
   
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   
 
  
   
     if (my_rank == MASTER_RANK){
		i_am_the_master = TRUE;
     }
		 
	   log_file = fopen ("./log.txt", "a");
	
    
    //------------
    // if Mode is validation, then only do 1 iteration. 
    //------------
    if(MODE==1 ){
        maxMacroMesoIterationCount=NUM_MACRO_MESO_ITER;
    }else if(MODE==2 || MODE==4){
         maxMacroMesoIterationCount=1;
    } else if( MODE==3){
        maxMacroMesoIterationCount=100;
    }
   
   for (int k = START_MACRO_MESO_ITER;k<maxMacroMesoIterationCount+1;k++){
     
        // -------------------------------------------------------
        //
        // MASTER code!! I am the master, you are the padawan!
        //
        // -------------------------------------------------------
        if (i_am_the_master) {
             fprintf(log_file, "\n--------\nMaster Node Starint Iteration%d\n-----------\n", k);
             fflush(log_file);   
            if(MODE==1 || MODE==3){
                MultiscaleMasterNode(k,pool_size,my_rank);
            }
            else if(MODE==2 || MODE==4){
                 MesoValidationMasterNode(k,pool_size,my_rank);
            }     

                  
         }         
         else
        {           
            // -------------------------------------------------------
            //
            // Padawan code!! Obey the master
            //
            // -------------------------------------------------------  
                if (my_rank == SUB_MASTER_RANK && (MODE==1 || MODE==3)){
                    fprintf(log_file, "\n--------\nNew Macro Meso Iteration %d\n-----------\n", k);
                    fflush(log_file);     
                }            
               
                // ------------------------
                // 1. Pause till master code sends out a message that he is done with macro level
                // ------------------------
                MPI_Barrier(MPI_COMM_WORLD);
                
                 // If we have more nodes than elements, then ignore nelm+1 nodes and have them do nothing
                 // they must still participate in the Barriers. 
                //if(my_rank<nelm+1)
                //{
                    
                    //MPI_Bcast(b, COLS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
                   // fprintf(log_file, "Got past the barrier as  work\n");
                    //fflush(log_file);
                      
                    //-------------------------------------
                    // 2. Get a job from the master
                    //-------------------------------------
                    int elementNumber;
                    int okMessage = 1;
                    MPI_Recv(&elementNumber, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    fprintf(log_file, "received a message from %d, tag %d\n", status.MPI_SOURCE, status.MPI_TAG);
                    fflush(log_file);
                      
                    // While we keep getting new jobs. 
                    while (status.MPI_TAG != TERMINATION_TAG) { /* The job is not finished */
                        
                         
                        // Run matlab code for element = elementNumber
                        // callMatlab(int mode, int macro_meso_iteration, int element)
                        callMatlab(100,k, elementNumber);
                        fprintf(log_file, "Finished Matlab Meso Structure Design. ElementNum = %d, MacroMesoIteration %d\n", elementNumber, k);
                        fflush(log_file);
                      //  printf("Finished Element Optimization for e = %d withstatus %d\n", elementNumber,status);
                         
                        MPI_Send (&okMessage, 1, MPI_INT, MASTER_RANK, elementNumber, MPI_COMM_WORLD);
                        // fprintf(log_file, "sent row %d to %d\n", row, MASTER_RANK);
                        // fflush(log_file);
                        MPI_Recv(&elementNumber, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        fprintf(log_file, "received a message from %d, tag %d\n", status.MPI_SOURCE, status.MPI_TAG);
                        fflush(log_file);
                    }
                    printf( "process with rank %d exiting on  tag %d\n", my_rank,status.MPI_TAG);
                    //fflush(log_file);
              //  } else{
                //      printf( "My rank is higher than max jobs needed. Rank = %d\n", my_rank);
                //}// end  if(my_rank<nel+1)
                  
                MPI_Barrier(MPI_COMM_WORLD); // force all to wait till all elements are finished
                   
                   
                   
                if (my_rank == SUB_MASTER_RANK && (MODE==1 || MODE==3)){
                     // --------------------------------------------
                    // 4. run the code to combine everything together
                    // --------------------------------------------
                    fprintf(log_file, "SUB_MASTER is now running Mode 202 Recombine everything into complete structure\n");
                    fflush(log_file);
                    printf("SUB_MASTER --> All elements done. Recombine the entire matrix");
                    // callMatlab(int mode, int macro_meso_iteration, int element)
                    callMatlab(202,k, 1);
                    fprintf(log_file, "SUB_MASTER Finished");
                    fflush(log_file);
                    printf("SUB_MASTER --> Finished Recombine the entire matrix");
                    
                    if(MODE==1){
                          printf("SUBMASTER --> Calculate Objective Using D_sub values\n");			
                          callMatlab(90,k, 1);   
                            printf("SUBMASTER ---------> Finished\n");			
                    }
                   
                } 
              
        } // end if  (i_am_the_master)
        
        
      
       
   } // End the master for loop over the macro-meso iterations
   
   if (i_am_the_master && MODE==1) {
        // Call to finilize the design and track metrics. 
       // callMatlab(int mode, int macro_meso_iteration, int element)
       callMatlab(60,NUM_MACRO_MESO_ITER+1, 1);
       callMatlab(200,NUM_MACRO_MESO_ITER+1, 1);
   }

   MPI_Finalize ();   
}


int callMatlab(int mode, int macro_meso_iteration, int element){
	
	char command[200];
	char commandHuman[200];
	int useCommandLineArgs = 1;
	//double w1 = 1;
	int w1 = 1;
    int returnValue=0;
    
    char name[200];
    int length;
    
    MPI_Get_processor_name(name, &length);

   sprintf(command, "./combinedTopologyOptimization %d %d %d %d %d ",useCommandLineArgs, w1, macro_meso_iteration, mode, element);
   sprintf(commandHuman, "combinedTopologyOptimization , macroMeso =  %d , mode = %d, element =  %d , Processor %s", macro_meso_iteration, mode, element,name);
   puts(commandHuman);
  // sleep(1);
  // Try 5 times to get the command right with a return of 0
  // Matlab sometimes randomely fails. 
   for (int i = 0; i < 5; i++) {
        returnValue=system(command);
       // printf( "Finished Command With Exit value %d \n", returnValue);
        if(returnValue!=0){
               printf( "Command Failed with code %d try again: %s \n", returnValue, commandHuman);
              //  sleep(5);
        } else {
            break;
        }
   }
   
     if(returnValue!=0){
         exit(5);
     }
   
  return returnValue;		
}

// ----------------------------------------------
// =============================================
//              MULTISCALE OPTIMIZATION  
//              MASTER NODE
// =============================================
// ----------------------------------------------
int MultiscaleMasterNode(int macro_meso_iteration,int pool_size,int my_rank ){
    int  nelm ;//= NELX_MACRO*NELY_MACRO;	
    struct timespec tim, tim2;
    tim.tv_sec = 0;
    tim.tv_nsec = 100000000; // 0.1 second
    int sender, destination; 
     MPI_Status status;    

      
    // ----------------------------------------------
    // 1. run macro level on master
    // ----------------------------------------------
    // callMatlab(int mode, int macro_meso_iteration, int element)
    printf("Master code running Macro Gradient Materail Optimization for Interation %d\n",macro_meso_iteration);
    
    if( MODE==1){
          callMatlab(60,macro_meso_iteration, 1);
    } else if(MODE==3){
         callMatlab(1,macro_meso_iteration, 1);
    }
  

    callMatlab(200,macro_meso_iteration, 1);
    //sleep(1);
    
     nelm= DetermineNumberOfElements(macro_meso_iteration,0 );
    
    
    // -----------------
    // 2 Pause everyone until we are ready for the element iterations. 
    // -----------------
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Bcast(b, ROWS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    
    // --------------------------------------------
    // 3. Send out jobs. 
    // --------------------------------------------
    //int count = 0;
    int elementNumber=1;
     
   
    for (int destination = 0; destination < pool_size; destination++) {
         if (destination != my_rank) {
             if(elementNumber<nelm+1){
                     
                //for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];
                //MPI_Send(int_buffer, COLS, MPI_INT, destination, count,
                //		 MPI_COMM_WORLD);
                MPI_Send(&elementNumber, 1, MPI_INT, destination, elementNumber, MPI_COMM_WORLD);
                printf("sent Element Job %d of %d to %d\n", elementNumber, nelm,destination);
                //count = count + 1;
                elementNumber=elementNumber+1;
                nanosleep(&tim , &tim2) ;// Matlab is giving me some weird Bus error. Maybe trying to read the same files all at the same time is
                  // causing issues? I'll slow down the process by making them 1 second apart. 
                  // https://stackoverflow.com/questions/7684359/how-to-use-nanosleep-in-c-what-are-tim-tv-sec-and-tim-tv-nsec
             } else {
                  MPI_Send(NULL, 0, MPI_INT, destination, TERMINATION_TAG, MPI_COMM_WORLD);
                   printf("No more jobs. Terminated process %d with tag %d\n", destination, TERMINATION_TAG);
                 
             }
         }
      }
      
    printf("\nSent Intial jobs to workers\n");

    int valueReturned;
    int rec_elementNumber;
    
    
      for (int i = 0; i < nelm; i++) {
          // ------------------------------------
          // Wait for a message saying someone finished their job. 
          // ------------------------------------
         MPI_Recv (&valueReturned, 1, MPI_INT, MPI_ANY_SOURCE,   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
         sender = status.MPI_SOURCE;
         rec_elementNumber = status.MPI_TAG;
         //c[row] = int_buffer[0];
        // printf("received  results from element %d from %d\n", rec_elementNumber, sender);
         destination = sender;
         
          // ------------------------------------
          // If more jobs, then send a new one. , plus 1 becuase we started at 1
          // ------------------------------------
         if (elementNumber < nelm+1) {
            //for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];					
            MPI_Send(&elementNumber, 1, MPI_INT, destination, elementNumber, MPI_COMM_WORLD);
            printf("sent Element Job %d of %d to %d\n", elementNumber, nelm,destination);
            //count = count + 1;
            elementNumber=elementNumber+1;
            
         }
         else {
               // ------------------------------------
              // If NO jobs, then send a termination messsage
              // ------------------------------------
            MPI_Send(NULL, 0, MPI_INT, destination, TERMINATION_TAG, MPI_COMM_WORLD);
            printf("terminated process %d with tag %d, i value is %d of %i \n", sender, TERMINATION_TAG, i,nelm);
         }
      }
      
       
        // Still need to wait, since all process must reach barrier
        // --------------------------------------------
        printf( "At position before Master Barrier\n");
        MPI_Barrier(MPI_COMM_WORLD); // force all to wait till all elements are finished
         
        
        // --------------------------------------------
        // 5. Get all the Exx, Eyy, theta values from the D_sub matrixes. Save all in csv files. 
        // --------------------------------------------
        if( MODE==1){
            printf("MASTER --> Get all the Exx, Eyy, theta values from the D_sub matrixes\n");			
            callMatlab(203,macro_meso_iteration, 1);   
             printf("MASTER -->Finished extracting Macro Vars from D_sub values. Finished Mode 203\n");			
             
           
        }
        
        return 1;

}


// ----------------------------------------------
// =============================================
//              MESO VALIDATION 
//              MASTER NODE
// =============================================
// ----------------------------------------------
int MesoValidationMasterNode(int macro_meso_iteration,int pool_size ,int my_rank){
	int k = macro_meso_iteration;
    int  nelm ;//= NELX_MACRO*NELY_MACRO;
    struct timespec tim, tim2;
   tim.tv_sec = 0;
   tim.tv_nsec = 500000;
   int sender,destination;
   // int terminateTag;
     MPI_Status status;
    // terminateTag =nelm+10;
    
      
    // ----------------------------------------------
    // 1. run macro level on master
    // ----------------------------------------------
    // callMatlab(int mode, int macro_meso_iteration, int element)
    if(MODE ==2){
        printf("Generate Targets for Meso Validation%d\n",k);
        callMatlab(111,k, 1); // Manully upload the files. 
    } else if (MODE ==4){
         printf("Generate Targets Pseudo Strian and Density : mode = 113%d\n",k);
        callMatlab(113,k, 1); // Manully upload the files. 
    }
    
    nelm= DetermineNumberOfElements(macro_meso_iteration,0 );
   // nelm=8000;
   
    
    
    // -----------------
    // 2 Pause everyone until we are ready for the element iterations. 
    // -----------------
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Bcast(b, ROWS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    
    // --------------------------------------------
    // 3. Send out jobs. 
    // --------------------------------------------
    //int count = 0;
    int elementNumber=1;
    //int startJob=
   
    for (int destination = 0; destination < pool_size; destination++) {
         if (destination != my_rank) {
             if(elementNumber<nelm+1){
                     
                //for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];
                //MPI_Send(int_buffer, COLS, MPI_INT, destination, count,
                //		 MPI_COMM_WORLD);
                MPI_Send(&elementNumber, 1, MPI_INT, destination, elementNumber, MPI_COMM_WORLD);
                printf("sent Element Job %d of %d to %d\n", elementNumber, nelm,destination);
                //count = count + 1;
                elementNumber=elementNumber+1;
                nanosleep(&tim , &tim2) ;// Matlab is giving me some weird Bus error. Maybe trying to read the same files all at the same time is
                  // causing issues? I'll slow down the process by making them 1 second apart. 
                  // https://stackoverflow.com/questions/7684359/how-to-use-nanosleep-in-c-what-are-tim-tv-sec-and-tim-tv-nsec
             } else{
                  MPI_Send(NULL, 0, MPI_INT, destination, TERMINATION_TAG, MPI_COMM_WORLD);
                  printf("terminated process %d with tag %d\n", sender, TERMINATION_TAG);
             }
         }
      }
      
    printf("\nSent Intial jobs to workers\n");

    int valueReturned;
    int rec_elementNumber;
    
    
      for (int i = 0; i < nelm; i++) {
          // ------------------------------------
          // Wait for a message saying someone finished their job. 
          // ------------------------------------
         MPI_Recv (&valueReturned, 1, MPI_INT, MPI_ANY_SOURCE,   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
         sender = status.MPI_SOURCE;
         rec_elementNumber = status.MPI_TAG;
         //c[row] = int_buffer[0];
        // printf("received  results from element %d from %d\n", rec_elementNumber, sender);
         destination = sender;
         
          // ------------------------------------
          // If more jobs, then send a new one. , plus 1 becuase we started at 1
          // ------------------------------------
         if (elementNumber < nelm+1) {
            //for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];					
            MPI_Send(&elementNumber, 1, MPI_INT, destination, elementNumber, MPI_COMM_WORLD);
            printf("sent Element Job %d of %d to %d\n", elementNumber, nelm,destination);
            //count = count + 1;
            elementNumber=elementNumber+1;
            
         }
         else {
               // ------------------------------------
              // If NO jobs, then send a termination messsage
              // ------------------------------------
            MPI_Send(NULL, 0, MPI_INT, destination, TERMINATION_TAG, MPI_COMM_WORLD);
             printf("terminated process %d with tag %d, i value is %d of %i \n", sender, TERMINATION_TAG, i,nelm);
         }
      }
      
       
        // Still need to wait, since all process must reach barrier
        // --------------------------------------------
        printf( "At position before Master Barrier\n");
         MPI_Barrier(MPI_COMM_WORLD); // force all to wait till all elements are finished
         
        
        // --------------------------------------------
        // 5. REad validation results
        // --------------------------------------------
        if(MODE ==2){
            printf("MASTER --> Interprete the Meso validation results, ");			
            callMatlab(112,k, 1);   
        } else {
            printf("MASTER --> Compile PSeudo strain and density information. MODE 114 ");			
            callMatlab(114,k, 1);   
        }
        
        return 1;

}

int DetermineNumberOfElements(int macro_meso_iteration,int folderNumber ){
   int bufSize = 1000000;
    char str[bufSize];
    FILE * file;

    char outname[200] ;
   
  
    sprintf(outname,"./out%i/SIMPdensityfield%i.csv",folderNumber,macro_meso_iteration);
     printf("Attempting to DetermineNumberOfElements from %s\n",outname);
   //  sprintf(outname,"../out%i/SIMPdensityfield%i.csv",folderNumber,macro_meso_iteration);
    file = fopen( outname, "r");
    int rowCount = 0;
    int commaCount=0;
    
    if (file) {
        while (fscanf(file, "%s", str)!=EOF){
            int charCount = strlen(str);
            rowCount++;
             printf("rowCount %i : %s\n",rowCount, str);
             
             if(rowCount==1){
                for(int m=0; m<charCount; m++) {
                    char temp = str[m];
                    if(temp== ',') {
                        commaCount++;
                    } 
                }
             }
        }
           
        fclose(file);
    } else {
        printf("Error reading file\n");
    }
    
    printf("Row Count = %i, Comma Count = %i, total Values = %i\n",rowCount,commaCount,rowCount*(commaCount+1));
    
    return commaCount,rowCount*(commaCount+1);
	
}  