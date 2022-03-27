#include "mpi.h"
#include <stdio.h>

/* problem parameters */
#define f(x) ((x) * (x))
#define numberRects 50
#define lowerLimit 2.0
#define upperLimit 5.0
int main(int argc, char *argv[])
{
    /*Variáveis MPI */
    int dest, noProcesses, processId, src, tag;
    MPI_Status status;
    /* Variáveis do problema */
    int i, k;
    double area, at, height, lower, width, total, range;
    /* Inicializacção do MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    
    /*Ajustando o tamanho para o subproblema*/
    range = (upperLimit - lowerLimit) / noProcesses;
    width = range / numberRects;
    lower = lowerLimit + range * processId;
    /*calculando área para o subproblema*/
    area = 0.0;
    if (processId == 0) scanf("%d", &k);
    printf("meu processo e %d \n", processId);
    for (i = 0; i < numberRects; i++)
    {
        at = lower + i * width + width / 2.0;
        height = f(at);
        area = area + width * height;
    }

    /* coletando informações e imprimindo resultado */
    tag = 0;
    if (processId == 0)
    /* Se o rank é zero ele coleta os resultados */
    {
        total = area;
        for (src = 1; src < noProcesses; src++)
        {
            MPI_Recv(&area, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD,
                     &status);
            total = total + area;
        }
        fprintf(stderr, "A area entre %f e %f e: %f\n", lowerLimit, upperLimit,
                total);
    }

    else
    /* Todos os outros processos somente enviam os
    resultados */
    {
        dest = 0;
        MPI_Send(&area, 1, MPI_DOUBLE, dest, tag,
                 MPI_COMM_WORLD);
    };
    /* Finalizando */
    MPI_Finalize();
    return 0;
}