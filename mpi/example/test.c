#include "mpi.h"
#include <stdio.h>

int main( int argc, char * argv[ ] )
{
int processId; /* rank dos processos */
int noProcesses; /* Número de processos */
int nameSize; /* Tamanho do nome */
char computerName[MPI_MAX_PROCESSOR_NAME];
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);
MPI_Comm_rank(MPI_COMM_WORLD, &processId);
MPI_Get_processor_name(computerName, &nameSize);
fprintf(stderr,"Hello from process %d on %s\n. Total process: %d", processId,
computerName, noProcesses);
MPI_Finalize();
return 0; 
}