#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include<stdlib.h>

#define TRUE 1
#define FALSE 0
#define MASTER_RANK 0
#define SUB_MASTER_RANK 1
#define NUM_MACRO_MESO_ITER 2
#define NELX_MACRO  39
#define NELY_MACRO 21

void callMatlab(int mode, int macro_meso_iteration, int element);

int main ( int argc, char **argv )
{
   int pool_size, my_rank, destination;
   int i_am_the_master = FALSE; 
   int nelm; //Number of elements
      
 
  
   MPI_Status status;
   int terminateTag;
    FILE *log_file;
	
	int sender;
   
   // For each macro loop 
   // run macro level on master
   // broad cast (which will force everyone to wait)
   // Begin sending out jobs untill finished
   // run the combine together version of the matlab code. 
   
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   nelm = NELX_MACRO*NELY_MACRO;
   terminateTag = nelm +10;
   
    if (my_rank == MASTER_RANK){
		i_am_the_master = TRUE;
	}else {
		 

	   log_file = fopen ("/tmp/gustav_log", "w");
	}
   
   for (int k = 1;k<NUM_MACRO_MESO_ITER+1;k++){
	 
		// -------------------------------------------------------
		//
		// MASTER code!! I am the master, you are the padawan!
		//
		// -------------------------------------------------------
		if (i_am_the_master) {
			// ----------------------------------------------
			// 1. run macro level on master
			// ----------------------------------------------
			// callMatlab(int mode, int macro_meso_iteration, int element)
			printf("Master code running Macro Gradient Materail Optimization for Interation %d\n",k);
			callMatlab(60,k, 1);
            callMatlab(200,k, 1);
			//sleep(1);
			
			
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
			for (destination = 0; destination < pool_size; destination++) {
				 if (destination != my_rank) {
					//for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];
					//MPI_Send(int_buffer, COLS, MPI_INT, destination, count,
					//		 MPI_COMM_WORLD);
					MPI_Send(&elementNumber, 1, MPI_INT, destination, elementNumber, MPI_COMM_WORLD);
					printf("sent Element Job %d to %d\n", elementNumber, destination);
					//count = count + 1;
					elementNumber=elementNumber+1;
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
					printf("sent Element Job %d to %d\n", elementNumber, destination);
					//count = count + 1;
					elementNumber=elementNumber+1;
					
				 }
				 else {
					   // ------------------------------------
					  // If NO jobs, then send a termination messsage
					  // ------------------------------------
					MPI_Send(NULL, 0, MPI_INT, destination, terminateTag, MPI_COMM_WORLD);
					printf("terminated process %d with tag %d, i value is %d\n", sender, terminateTag, i);
				 }
			  }
			  
			    // --------------------------------------------
				// 4. MOved to SubMaster Process
				// run the code to combine everything together
				//
				// Still need to wait, since all process must reach barrier
				// --------------------------------------------
			
				 MPI_Barrier(MPI_COMM_WORLD); // force all to wait till all elements are finished
			  	
			  
		 }
		 
		 else
		{
		   
			// -------------------------------------------------------
			//
			// Padawan code!! Obey the master
			//
			// -------------------------------------------------------
		   
		  			
		   
		   // ------------------------
		   // 1. Pause till master code sends out a message that he is done with macro level
		   // ------------------------
			  MPI_Barrier(MPI_COMM_WORLD);
			  //MPI_Bcast(b, COLS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
			  fprintf(log_file, "Got past the barrier as  work\n");
			  fflush(log_file);
			  
			  //-------------------------------------
			  // 2. Get a job from the master
			  //-------------------------------------
			  int elementNumber;
			  int okMessage = 1;
			  MPI_Recv(&elementNumber, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			  fprintf(log_file, "received a message from %d, tag %d\n", status.MPI_SOURCE, status.MPI_TAG);
			  fflush(log_file);
			  
			  // While we keep getting new jobs. 
			  while (status.MPI_TAG != terminateTag) { /* The job is not finished */
				
				 
				// Run matlab code for element = elementNumber
				// callMatlab(int mode, int macro_meso_iteration, int element)
				callMatlab(100,k, elementNumber);
				//printf("Finished Element Optimization for e = %d\n", elementNumber);
				 
				 MPI_Send (&okMessage, 1, MPI_INT, MASTER_RANK, elementNumber, MPI_COMM_WORLD);
				// fprintf(log_file, "sent row %d to %d\n", row, MASTER_RANK);
				// fflush(log_file);
				 MPI_Recv(&elementNumber, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				 //fprintf(log_file, "received a message from %d, tag %d\n", status.MPI_SOURCE, status.MPI_TAG);
				 //fflush(log_file);
			  }
			  printf( "exiting on  tag %d\n", status.MPI_TAG);
			  //fflush(log_file);
			  
			   MPI_Barrier(MPI_COMM_WORLD); // force all to wait till all elements are finished
			   
			   
			   
			   if (my_rank == SUB_MASTER_RANK){
				    // --------------------------------------------
					// 4. run the code to combine everything together
					// --------------------------------------------
					printf("SUB_MASTER --> All elements done. Recombine the entire matrix");
					// callMatlab(int mode, int macro_meso_iteration, int element)
					callMatlab(202,k, 1);
			   }
			   
			
		}
		
		
	  
	   
   } // End the master for loop over the macro-meso iterations

  

   MPI_Finalize ();
   
   // Call to finilize the design and track metrics. 
   // callMatlab(int mode, int macro_meso_iteration, int element)
	callMatlab(12,NUM_MACRO_MESO_ITER, 1);
}


void callMatlab(int mode, int macro_meso_iteration, int element){
	
	char command[200];
	char commandHuman[200];
	int useCommandLineArgs = 1;
	//double w1 = 1;
	int w1 = 1;

   sprintf(command, "./combinedTopologyOptimization %d %d %d %d %d ",useCommandLineArgs, w1, macro_meso_iteration, mode, element);
   sprintf(commandHuman, "combinedTopologyOptimization , macroMeso =  %d , mode = %d, element =  %d ", macro_meso_iteration, mode, element);
   puts(commandHuman);
  // sleep(1);
  system(command);
	
	
	
	
}