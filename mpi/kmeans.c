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
	}

	// all_script = clock();
	printf("Process id %d has started! \n", processId);

	if (processId == MAIN_PROCESS)
	{
		for (i = 0; i < k; i++)
			scanf("%lf %lf %lf", mean + i * DIM, mean + i * DIM + 1, mean + i * DIM + 2);
		for (i = 0; i < n; i++)
			scanf("%lf %lf %lf", x + i * DIM, x + i * DIM + 1, x + i * DIM + 2);
		for (int i = 1; i < totalSizeProcess; i++)
		{
			printf("Main process sending to %d \n", i);
			MPI_Send(&k, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(&n, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(mean, DIM * k, MPI_DOUBLE, i, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(x, DIM * n, MPI_DOUBLE, i, TAG_COMMOM, MPI_COMM_WORLD);
		}
	}
	else
	{
		printf("Side process %d receiving... \n", processId);
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

	// 100
	// 5
	// 100 / 4
	// 25
	// start 1 - 0
	// start 2 - 25
	// end 1 - 25
	// end 2 - 50
	chunk_size = n / (totalSizeProcess - 1);
	start = chunk_size * (processId - 1);
	end = chunk_size * processId;

	for (i = 0; i < n; i++)
		cluster[i] = 0;

	// Dividir os chunks de cada processamento

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
			printf("Envia os clusters para os slaves \n");
			for (int i = 1; i < totalSizeProcess; i++)
			{
				MPI_Send(mean, DIM * k, MPI_DOUBLE, i, TAG_COMMOM, MPI_COMM_WORLD);
			}
		} else {
			// printf("cluster[0] = %d", cluster[0]);
			MPI_Recv(mean, DIM * k, MPI_DOUBLE, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD, &status);
			// printf("Processo %d recebe os clusters da master \n", processId);
			// printf("cluster[0] = %d", cluster[0]);
		}

		// Se for processo nao principal, processa as distancias euclidiana
		if (processId != MAIN_PROCESS)
		{
			// printf("Processa as distancias \n");
			for (i = start; i < end; i++)
			{
				dmin = -1;
				// printf("Meu processo e %d e o indice e %d end e %d \n", processId, i, end);
				color = cluster[i];
				// printf("1 \n");
				{
					for (c = 0; c < k; c++)
					{
						dx = 0.0;
						// printf("2 \n");
						for (j = 0; j < DIM; j++)
						{
							dx += (x[i * DIM + j] - mean[c * DIM + j]) * (x[i * DIM + j] - mean[c * DIM + j]);
							// printf("3 \n");
						}

						if (dx < dmin || dmin == -1)
						{
							// printf("4 \n");
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
			printf("Processo %d envia flips %d \n", processId, flips);
			MPI_Send(&flips, 1, MPI_INT, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD);
			MPI_Send(cluster, n, MPI_INT, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD);
			printf("Envia as cores calculadas para a master \n");
		}
		else  // Se for processo principal, pega todos os flips e atualiza as cores dos clusters
		{
			total_flips = 0;
			for (int i = 1; i < totalSizeProcess; i++)
			{
				MPI_Recv(&total_flips, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD, &status);
				flips += total_flips;
				// printf("total flips %d \n", total_flips);
				MPI_Recv(cp_cluster, DIM * n, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD, &status);
				printf("chunk size %d\n", chunk_size);
				for (int j = 0; j < chunk_size; j++)
				{
					slot = ((i-1) * chunk_size) + j;
					// printf("slot %d \n", slot);
					cluster[slot] = cp_cluster[slot];
				}
			}
			printf("Atualiza master com clusters e flips %d\n", flips);
		}


		// Atualiza flips de todos os processos
		if (processId != MAIN_PROCESS)
		{
			MPI_Recv(&flips, 1, MPI_INT, MAIN_PROCESS, TAG_COMMOM, MPI_COMM_WORLD, &status);
			continue;
		} else {
			for (int i = 1; i < totalSizeProcess; i++)
			{
				MPI_Send(&flips, 1, MPI_INT, i, TAG_COMMOM, MPI_COMM_WORLD);
			}
			// printf("Flips do processo principal %d \n", flips);
		}

		// Pegar os resultados de cada clusters

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
		// printf("Mean distance took %.4lf seconds to run.\n", (float)(clock() - start) / CLOCKS_PER_SEC );
		// acc_mean += clock() - start;
	}
	for (i = 0; i < k; i++)
	{
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", mean[i * DIM + j]);
		printf("\n");
	}
// printf("Acc distance %.4lf , Acc sum %.4lf Acc mean %.4lf.\n",
// 	(float)(acc_distance/CLOCKS_PER_SEC),
// 	(float)(acc_sum/CLOCKS_PER_SEC),
// 	(float)(acc_mean/CLOCKS_PER_SEC)
// );
// printf("All script took %.4lf seconds to run.\n", (float)(clock() - all_script) / CLOCKS_PER_SEC );
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
