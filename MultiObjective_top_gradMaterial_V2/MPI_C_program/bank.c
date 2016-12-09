#include <stdio.h>
#include <mpi.h>
#define COLS 100
#define ROWS 100
#define TRUE 1
#define FALSE 0
#define MASTER_RANK 0

int main ( int argc, char **argv )
{
   int pool_size, my_rank, destination;
   int i_am_the_master = FALSE; 
   int a[ROWS][COLS], b[ROWS], c[ROWS], i, j;
   int int_buffer[BUFSIZ];
   MPI_Status status;
   

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   if (my_rank == MASTER_RANK) i_am_the_master = TRUE;

   if (i_am_the_master) {

      int row, count, sender;

      for (i = 0; i < COLS; i++) {
         b[i] = 1;
         for (j = 0; j < ROWS; j++){
			 a[i][j] = i;
		 }
      }

      MPI_Bcast(b, ROWS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

      count = 0;
      for (destination = 0; destination < pool_size; destination++) {
         if (destination != my_rank) {
            for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];
            MPI_Send(int_buffer, COLS, MPI_INT, destination, count,
                     MPI_COMM_WORLD);
            printf("sent row %d to %d\n", count, destination);
            count = count + 1;
         }
      }

      for (i = 0; i < ROWS; i++) {
         MPI_Recv (int_buffer, BUFSIZ, MPI_INT, MPI_ANY_SOURCE, 
                   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
         sender = status.MPI_SOURCE;
         row = status.MPI_TAG;
         c[row] = int_buffer[0];
         printf("\treceived row %d from %d\n", row, sender);
         if (count < ROWS) {
            for (j = 0; j < COLS; j++) int_buffer[j] = a[count][j];
            MPI_Send(int_buffer, COLS, MPI_INT, sender, count,
                     MPI_COMM_WORLD);
            printf("sent row %d to %d\n", count, sender);
            count = count + 1;
         }
         else {
            MPI_Send(NULL, 0, MPI_INT, sender, ROWS, MPI_COMM_WORLD);
            printf("terminated process %d with tag %d\n", sender, ROWS);
         }
      }
   }
   else { /* I am not the master */

      int sum, row;
      FILE *log_file;

      log_file = fopen ("/tmp/gustav_log", "w");

      MPI_Bcast(b, COLS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
      fprintf(log_file, "received broadcast from %d\n", MASTER_RANK);
      fflush(log_file);
      MPI_Recv(int_buffer, COLS, MPI_INT, MASTER_RANK, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &status);
      fprintf(log_file, "received a message from %d, tag %d\n",
                   status.MPI_SOURCE, status.MPI_TAG);
      fflush(log_file);
      while (status.MPI_TAG != ROWS) { /* The job is not finished */
         row = status.MPI_TAG; sum = 0;
         for (i = 0; i < COLS; i++) sum = sum + int_buffer[i] * b[i];
         int_buffer[0] = sum;
         MPI_Send (int_buffer, 1, MPI_INT, MASTER_RANK, row, MPI_COMM_WORLD);
         fprintf(log_file, "sent row %d to %d\n", row, MASTER_RANK);
         fflush(log_file);
         MPI_Recv (int_buffer, COLS, MPI_INT, MASTER_RANK, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &status);
         fprintf(log_file, "received a message from %d, tag %d\n",
                 status.MPI_SOURCE, status.MPI_TAG);
         fflush(log_file);
      }
      fprintf(log_file, "exiting on  tag %d\n", status.MPI_TAG);
      fflush(log_file);
   }

   MPI_Finalize ();
}