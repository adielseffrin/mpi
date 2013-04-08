#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mpe.h"
#include <math.h>
#include <string.h>

void
printMatrix (double *M, int m, int n)
{
  int lin, col;
  for (lin = 0; lin < m; lin++)
    {
      for (col = 0; col < n; col++)
	printf ("%.2f \t", M[((lin * n) + col)]);
      printf ("\n");
    }
}

void
imprimeMatriz_all (const char *const message, double *M, int m, int n)
{
  int rank, size, i;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Barrier (MPI_COMM_WORLD);
  if (!rank)
    printf ("%s\n", message);
  MPI_Barrier (MPI_COMM_WORLD);
  for (i = 0; i < size; i++)
    {
      if (rank == i)
	{
	  printf ("rank = %d\n", i);
	  printMatrix (M, m, n);
	}
      MPI_Barrier (MPI_COMM_WORLD);
      printf ("\n");
      MPI_Barrier (MPI_COMM_WORLD);
    }
  MPI_Barrier (MPI_COMM_WORLD);
}




double *
allocateMatrix (int m, int n)
{
  double *M;
  M = (double *) malloc (m * n * sizeof (double));

  if (!M)
    {
      fprintf (stderr, "Error: Not enough memory\n");
      exit (1);
    }

  bzero (M, (m * n) * sizeof (double));
  return M;
}




int
main (int argc, char *argv[])
{
  FILE *fp;
  char *nome;
  int rank, size;
  int m1, n1, m2, n2;
  int row, col, i, k, lines;
  double startwtime = 0.0, endwtime;
  double *M1, *M2, *M3, *M1_tmp, *M3_tmp, *M3_exp;
  int event1a, event1b, event2a, event2b,
    event3a, event3b, event4a, event4b, event5a, event5b, event6a, event6b;
  int event1, event2, event3, event4, event5, event6, event7;
  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  m1 = m2 = n1 = n2 = 1000;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

 // MPI_Get_processor_name (processor_name, &namelen);
//  fprintf (stderr, "Process %d running on %s\n", rank, processor_name);

#if defined( NO_MPI_LOGGING )
  MPE_Init_log ();
#endif

  MPE_Log_get_state_eventIDs (&event1a, &event1b);
  MPE_Log_get_state_eventIDs (&event2a, &event2b);
  MPE_Log_get_state_eventIDs (&event3a, &event3b);
  MPE_Log_get_state_eventIDs (&event4a, &event4b);
  MPE_Log_get_state_eventIDs (&event5a, &event5b);
  MPE_Log_get_state_eventIDs (&event6a, &event6b);


  if (rank == 0)
    {
      MPE_Describe_state (event1a, event1b, "Broadcast M3", "red");
      MPE_Describe_state (event2a, event2b, "Broadcast M2", "orange");
      MPE_Describe_state (event3a, event3b, "Scatter M1", "blue");
      MPE_Describe_state (event4a, event4b, "Compute", "green");
      MPE_Describe_state (event5a, event5b, "Gather", "brown");
      MPE_Describe_state (event6a, event6b, "Populate M3", "black");
    }

  /* Get event ID for Solo-Event(single timestamp object) from MPE */
  MPE_Log_get_solo_eventID (&event1);
  MPE_Log_get_solo_eventID (&event2);
  MPE_Log_get_solo_eventID (&event3);
  MPE_Log_get_solo_eventID (&event4);
  MPE_Log_get_solo_eventID (&event5);
  MPE_Log_get_solo_eventID (&event6);
  MPE_Log_get_solo_eventID (&event7);

  if (rank == 0)
    {
      MPE_Describe_event (event1, "Broadcast Post", "white");
      MPE_Describe_event (event2, "Scatter Post", "purple");
      MPE_Describe_event (event3, "Compute Start", "navy");
      MPE_Describe_event (event4, "Compute End", "gray");
      MPE_Describe_event (event5, "Gather Request", "pink");
      MPE_Describe_event (event6, "Populate Start", "yellow");
      MPE_Describe_event (event7, "Populate End", "olive");
    }



  lines = (int) ceil (m1 / (float) size);

  
  M1 = allocateMatrix (m1, n1);
  M3 = allocateMatrix (m1, n2);
  M2 = allocateMatrix (m2, n2);
  M1_tmp = allocateMatrix (lines, n1);
  M3_tmp = allocateMatrix (lines, n2);
  M3_exp = allocateMatrix (lines * size, n2);


  srand (100);

  for (col = 0; col < n1; col++)
    {
      for (row = 0; row < m1; row++)
	{
	  M1[((row * m1) + col)] = rand () % 50;
	  M2[(row * m2) + col] = rand () % 50;
	 // M3[(row * m2) + col] = 0;
	 // M3_tmp[(row * m2) + col] = 0;
	}
    //  for (row = 0; row < lines; row++)
	//M1_tmp[(row * m2) + col] = 0;
  //    for (row = 0; row < lines * size; row++)
//	M3_exp[(row * m2) + col] = 0;
    }

  if (rank == 0)
    startwtime = MPI_Wtime ();

  MPE_Log_event (event3a, 0, NULL);
  MPI_Scatter (M1, lines * n1, MPI_DOUBLE, M1_tmp, lines * n1, MPI_DOUBLE, 0,
	       MPI_COMM_WORLD);
  MPE_Log_event (event3b, 0, NULL);

  MPE_Log_event (event2, 0, NULL);

  MPE_Log_event (event1a, 0, NULL);
  MPI_Bcast (M2, m2 * n2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPE_Log_event (event1b, 0, NULL);

  MPE_Log_event (event1, 0, NULL);

  MPE_Log_event (event2a, 0, NULL);
  MPI_Bcast (M3_tmp, lines * n2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPE_Log_event (event2b, 0, NULL);

  MPE_Log_event (event1, 0, NULL);


 /*
  if (rank == size - 1)
  {
    i = (n1 % lines) * n1;
    if (i != 0)
      for (col = i; col <  lines * n1; col++)
	      M1_tmp[col] = 0;
  }
 
*/
/*
char mystr[200];

snprintf(mystr,200,   "arq_%d%d.out",rank,rank);

if((fp = fopen(mystr,"w"))!= NULL){
	for (i=0; i<lines; i++) {
	  for (k=0; k<n1; k++)
	    fprintf(fp,"%.2f \t", M1_tmp[((i*n1)+k)]);
	  fprintf(fp,"\n"); 
   }	
	fclose(fp);
}

*/
  MPE_Log_event (event3, 0, NULL);
  MPE_Log_event (event4a, 0, NULL);

  for (row = 0; row < lines; row++)
    {
      for (col = 0; col < n1; col++)
	    {
        double val = 0.0;
    	  for (k = 0; k < n1; k++)
        {
	        val += M1_tmp[(row * n1) + k] * M2[(k * n1) + col];
  	    }
	    M3_tmp[(row * n1) + col] = val;
    	}
    }
    
    
    
  MPE_Log_event (event4b, 0, NULL);
  MPE_Log_event (event4, 0, NULL);

  MPE_Log_event (event5a, 0, NULL);


  MPI_Gather (M3_tmp, lines * n1, MPI_DOUBLE, M3_exp, lines * n1, MPI_DOUBLE,
	      0, MPI_COMM_WORLD);



  MPE_Log_event (event5b, 0, NULL);

  MPE_Log_event (event5, 0, NULL);
  MPE_Log_event (event6, 0, NULL);
  MPE_Log_event (event6a, 0, NULL);


  for (col = 0; col < n1; col++)
  {
    for (row = 0; row < m1; row++)
  	{
  	  M3[((row * m1) + col)] = M3_exp[((row * m1) + col)];
  	}
  }
  
  
  MPE_Log_event (event6a, 0, NULL);
  MPE_Log_event (event7, 0, NULL);

  MPE_Log_sync_clocks ();

#if defined( NO_MPI_LOGGING )
  if (argv != NULL)
    MPE_Finish_log (argv[0]);
  else
    MPE_Finish_log ("matrizlog");
#endif
/*
if(rank==0){
printMatrix(M3_exp,lines*size,n1);
printf("\n---------------------------------------\n");
printMatrix(M3,m1,n1);
}


if(rank == 0){  
printf("matrix 1------------------------\n");
printMatrix(M1,m1,n1); 
printf("matrix 2------------------------\n");
printMatrix(M2,m2,n2);
printf("matrix 3------------------------\n");
printMatrix(M3,m1,n2);
}

if((fp = fopen("out1","w"))!= NULL){
	for (i=0; i<m1; i++) {
	  for (k=0; k<n1; k++)
	    fprintf(fp,"%.2f \t", M1[((i*n1)+k)]);
	  fprintf(fp,"\n"); 
   }	
	fclose(fp);
}
if((fp = fopen("out2","w"))!= NULL){
	for (i=0; i<m1; i++) {
	  for (k=0; k<n1; k++)
	    fprintf(fp,"%.2f \t", M2[((i*n1)+k)]);
	  fprintf(fp,"\n"); 
   }	
	fclose(fp);
}
*/
  if (rank == 0)
    {
      endwtime = MPI_Wtime ();
      printf ("wall clock time = %f\n", endwtime - startwtime);
    }
  MPI_Finalize ();
  return 0;
}
