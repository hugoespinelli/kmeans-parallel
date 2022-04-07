#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

#define TAG_COMMOM 1
#define DIM 3
#define MAIN_PROCESS 0

int main(int argc, char *argv[])
{
	int dest, totalSizeProcess, processId, start, end;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalSizeProcess);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	int i, j, k, n, c, chunk_size, slot;
	double dmin, dx;
	double *x, *mean, *sum;
	int *cluster, *cp_cluster, *count, color;
	int flips, total_flips;

	// Distribuição dos inputs da master com os slaves
	if (processId == MAIN_PROCESS)
	{
		scanf("%d", &k);
		scanf("%d", &n);
		x = (double *)malloc(sizeof(double) * DIM * n);
		mean = (double *)malloc(sizeof(double) * DIM * k);
		sum = (double *)malloc(sizeof(double) * DIM * k);
		cluster = (int *)malloc(sizeof(int) * n);
		cp_cluster = (int *)malloc(sizeof(int) * n);
		count = (int *)malloc(sizeof(int) * k);
		for (i = 0; i < k; i++)
			scanf("%lf %lf %lf", mean + i * DIM, mean + i * DIM + 1, mean + i * DIM + 2);
		for (i = 0; i < n; i++)
			scanf("%lf %lf %lf", x + i * DIM, x + i * DIM + 1, x + i * DIM + 2);
		for (int i = 1; i < totalSizeProcess; i++)
		{
			MPI_Send(&k, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(&n, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(mean, DIM * k, MPI_DOUBLE, i, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(x, DIM * n, MPI_DOUBLE, i, TAG_COMMOM, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(&k, 1, MPI_INT, 0, TAG_COMMOM, MPI_COMM_WORLD, &status);
		MPI_Recv(&n, 1, MPI_INT, 0, TAG_COMMOM, MPI_COMM_WORLD, &status);
		x = (double *)malloc(sizeof(double) * DIM * n);
		mean = (double *)malloc(sizeof(double) * DIM * k);
		sum = (double *)malloc(sizeof(double) * DIM * k);
		cluster = (int *)malloc(sizeof(int) * n);
		count = (int *)malloc(sizeof(int) * k);
		MPI_Recv(mean, DIM * k, MPI_DOUBLE, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD, &status);
		MPI_Recv(x, DIM * n, MPI_DOUBLE, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD, &status);
	}

	// Divisao de pedaços de elementos entre os slaves
	chunk_size = n / (totalSizeProcess - 1);
	start = chunk_size * (processId - 1);
	end = chunk_size * processId;

	for (i = 0; i < n; i++)
		cluster[i] = 0;

	// printf("Initialize variables centroids and numbers took %.4lf seconds to run.\n", (float)(clock() - start) / CLOCKS_PER_SEC );
	flips = n;
	while (flips > 0)
	{
		// printf("flips %d process %d \n", flips, processId);
		flips = 0;
		for (j = 0; j < k; j++)
		{
			count[j] = 0;
			for (i = 0; i < DIM; i++)
				sum[j * DIM + i] = 0.0;
		}

		if (processId == MAIN_PROCESS)
		{
			// Envia os centroides da master para o slave
			for (int i = 1; i < totalSizeProcess; i++)
			{
				MPI_Send(mean, DIM * k, MPI_DOUBLE, i, TAG_COMMOM, MPI_COMM_WORLD);
			}

			// Recebe os flips de atualização de cores dos clusters e as cores do cluster na master
			total_flips = 0;
			for (int i = 1; i < totalSizeProcess; i++)
			{
				MPI_Recv(&total_flips, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD, &status);
				flips += total_flips;

				// Atualização dos elementos nos clusters
				MPI_Recv(cp_cluster, DIM * n, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD, &status);
				for (int j = 0; j < chunk_size; j++)
				{
					slot = ((i - 1) * chunk_size) + j;
					cluster[slot] = cp_cluster[slot];
				}
			}

			// Manda os flips de atualização para os slaves ficarem sincronizados com a master
			for (int i = 1; i < totalSizeProcess; i++)
			{
				MPI_Send(&flips, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD);
			}
		}
		else
		{
			// Recebe os centroides recalculados da master
			MPI_Recv(mean, DIM * k, MPI_DOUBLE, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD, &status);

			// Inicia as contas de distancias euclidianas dos pontos ate os centroides
			for (i = start; i < end; i++)
			{
				dmin = -1;
				color = cluster[i];
				{
					for (c = 0; c < k; c++)
					{
						dx = 0.0;
						for (j = 0; j < DIM; j++)
						{
							dx += (x[i * DIM + j] - mean[c * DIM + j]) * (x[i * DIM + j] - mean[c * DIM + j]);
						}

						if (dx < dmin || dmin == -1)
						{
							color = c;
							dmin = dx;
						}
					}
				}

				if (cluster[i] != color)
				{
					flips += 1;
					cluster[i] = color;
				}
			}

			// Envia os flips de atualização e cores dos clusters para a master
			MPI_Send(&flips, 1, MPI_INT, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(cluster, n, MPI_INT, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD);
			
			// Recebe os flips de atualização da master
			MPI_Recv(&flips, 1, MPI_INT, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD, &status);
			continue;
		}


		// A master que efetua os calculos de atualização dos centroides
		for (i = 0; i < n; i++)
		{
			count[cluster[i]]++;
			for (j = 0; j < DIM; j++)
				sum[cluster[i] * DIM + j] += x[i * DIM + j];
		}
		for (i = 0; i < k; i++)
		{
			for (j = 0; j < DIM; j++)
			{
				mean[i * DIM + j] = sum[i * DIM + j] / count[i];
			}
		}
	}


	for (i = 0; i < k; i++)
	{
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", mean[i * DIM + j]);
		printf("\n");
	}

#ifdef DEBUG
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", x[i * DIM + j]);
		printf("%d\n", cluster[i]);
	}
#endif
	MPI_Finalize();
	return (0);
}
