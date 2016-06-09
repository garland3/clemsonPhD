#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include <stddef.h>

#define WORKTAG 1
#define DIETAG 2

#define COMPLETE 3
#define NOT_COMPLETE 4

#define NO_MORE_WORK 345

#define MAX_ITERATIONS 10
#define NUM_METHODS 6



typedef struct 
{
   int iteration;
   int method;
} unit_of_work_t;

/* Local functions */

static void master(void);
static void slave(void);

static unit_of_work_t get_next_work_item(int workArray[NUM_METHODS][MAX_ITERATIONS]);
static void process_results(unit_of_work_t result, int workArray[NUM_METHODS][MAX_ITERATIONS]);
static unit_of_work_t do_work(unit_of_work_t work);



extern void projectmain_(int *iteration, int *method,int *nx,int *ny);



int main(int argc, char **argv)
{
  int myrank;

  /* Initialize MPI */

  MPI_Init(&argc, &argv);
  
    
    
    
    
  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    master();
  } else {
    slave();
  }

  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


static void master(void)
{
  int ntasks, rank,i,j;
  unit_of_work_t work;
  unit_of_work_t result;
  MPI_Status status;
  
  int workArray[NUM_METHODS][MAX_ITERATIONS];
  
  /* create a type for struct work */
    // http://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-sending-over-mpi
    const int nitems=2;
    int          blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype mpi_work_type;
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(unit_of_work_t, iteration);
    offsets[1] = offsetof(unit_of_work_t, method);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_work_type);
    MPI_Type_commit(&mpi_work_type);
  
  
  // set everything in the array to be not complete
  for(i=0;i<NUM_METHODS;i++){
    for(j=0;j<MAX_ITERATIONS;j++){            
        workArray[i][j]= NOT_COMPLETE;
        }
   }
   

  /* Find out how many processes there are in the default
     communicator */

  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  /* Seed the slaves; send one unit of work to each slave. */

  for (rank = 1; rank < ntasks; ++rank) {

    /* Find the next item of work to do */

    work = get_next_work_item(workArray);

    /* Send it to each rank */

    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             mpi_work_type,           /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
  }

  /* Loop over getting new work requests until there is no more work
     to be done */

  work = get_next_work_item(workArray);
  while (work.method != NO_MORE_WORK) {

    /* Receive results from a slave */

    MPI_Recv(&result,           /* message buffer */
             1,                 /* one data item */
             mpi_work_type,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

    /* Send the slave a new work unit */
    

    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             mpi_work_type,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    /* Get the next unit of work to be done */

    
    work = get_next_work_item(workArray);
  }

  /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */

  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Recv(&result, 1, mpi_work_type, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  /* Tell all the slaves to exit by sending an empty message with the
     DIETAG. */

  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  }
}


static void slave(void)
{
  unit_of_work_t work;
  unit_of_work_t results;
  MPI_Status status;
  
  /* create a type for struct work */
    // http://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-sending-over-mpi
    const int nitems=2;
    int          blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype mpi_work_type;
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(unit_of_work_t, iteration);
    offsets[1] = offsetof(unit_of_work_t, method);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_work_type);
    MPI_Type_commit(&mpi_work_type);
  

  while (1) {

    /* Receive a message from the master */

    MPI_Recv(&work, 1, mpi_work_type, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    /* Check the tag of the received message. */

    if (status.MPI_TAG == DIETAG) {
      return;
    }

    /* Do the work */

    results = do_work(work);

    /* Send the result back */

    MPI_Send(&results, 1, mpi_work_type, 0, 0, MPI_COMM_WORLD);
  }
}

static int iteration_g =0;
static int method_g=0;

static unit_of_work_t get_next_work_item(int workArray[NUM_METHODS][MAX_ITERATIONS])
{
  /* Fill in with whatever is relevant to obtain a new unit of work
     suitable to be given to a slave. */
   
   unit_of_work_t work;
   
   if(iteration_g>=MAX_ITERATIONS){
        work.method = NO_MORE_WORK;
        return work;
   }

   work.method = method_g;
   work.iteration  = iteration_g;
   
   method_g++;
   if(method_g>=NUM_METHODS){
        iteration_g++;
        method_g=0;  
   }
   
   
   
   
  
   return work;
}

/*
Fill in the workArray at the location of the results to indicate that it is done. 
*/
static void process_results(unit_of_work_t result, int workArray[NUM_METHODS][MAX_ITERATIONS])
{
   // workArray[result.method][result.iteration]= COMPLETE;
  /* Fill in with whatever is relevant to process the results returned
     by the slave */
}


static unit_of_work_t do_work(unit_of_work_t work)
{
  /* Fill in with whatever is necessary to process the work and
     generate a result */
     int myrank;
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
     
     // The method numbers are 1,2,3,4,5,6 so that we need to add 1 to make it work correctly
     int temp_method = work.method +1;
     
     printf("\nSlave # %d starting new work  [iter,method] = %d, %d\n", myrank, work.iteration, temp_method);
     
     int nx = 25;
     int ny = 25;
     
    
     projectmain_(&(work.iteration), &temp_method,&nx,&ny);
     // projectmain_(&iteration,&method,&nx,&ny);
}