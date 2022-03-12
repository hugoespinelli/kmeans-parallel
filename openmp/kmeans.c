#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define DIM 3
#define NUM_THREADS INPUT_THREADS

int main(void) {
	clock_t start, all_script, acc_distance, acc_sum, acc_mean;
	acc_distance = 0.0;
	acc_sum = 0.0;
	acc_mean = 0.0;
	int i, j, k, n, c;
	double dmin, dx;
	double *x, *mean, *sum;
	int *cluster, *count, color;
	int flips;
	scanf("%d", &k);
	scanf("%d", &n);
	#ifdef _OPENMP
	omp_set_num_threads(NUM_THREADS);
	#endif
	x = (double *)malloc(sizeof(double)*DIM*n);
	mean = (double *)malloc(sizeof(double)*DIM*k);
	sum= (double *)malloc(sizeof(double)*DIM*k);
	cluster = (int *)malloc(sizeof(int)*n);
	count = (int *)malloc(sizeof(int)*k);
	all_script = clock();
	start = clock();
	for (i = 0; i<n; i++) 
		cluster[i] = 0;
	for (i = 0; i<k; i++)
		scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2);
	for (i = 0; i<n; i++)
		scanf("%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2);
	// printf("Initialize variables centroids and numbers took %.4lf seconds to run.\n", (float)(clock() - start) / CLOCKS_PER_SEC );
	flips = n;
	while (flips>0) {
		flips = 0;
		start = clock();
		for (j = 0; j < k; j++) {
			count[j] = 0; 
			for (i = 0; i < DIM; i++) 
				sum[j*DIM+i] = 0.0;
		}
		#pragma omp parallel private(i, c, j, color, dmin, dx)
		{
			#pragma omp for
			for (i = 0; i < n; i++) {
				dmin = -1; color = cluster[i];
				{
					for (c = 0; c < k; c++) {
						dx = 0.0;
						for (j = 0; j < DIM; j++) {
							dx +=  (x[i*DIM+j] - mean[c*DIM+j])*(x[i*DIM+j] - mean[c*DIM+j]);
						} 
						
						if (dx < dmin || dmin == -1) {
							color = c;
							dmin = dx;
						}
					}
				}

				if (cluster[i] != color) {
					# pragma omp atomic
					flips += 1;
					cluster[i] = color;
				}
			}

		}
		
		// printf("Sum all centroids distances took %.4lf seconds to run.\n", (float)(clock() - start) / CLOCKS_PER_SEC );
		// acc_distance += clock() - start;
		// start = clock();
	    for (i = 0; i < n; i++) {
			count[cluster[i]]++;
			for (j = 0; j < DIM; j++) 
				sum[cluster[i]*DIM+j] += x[i*DIM+j];
		}
		// printf("Sum all distances took %.4lf seconds to run.\n", (float)(clock() - start) / CLOCKS_PER_SEC );
		// acc_sum += clock() - start;
		// start = clock();
		for (i = 0; i < k; i++) {
			for (j = 0; j < DIM; j++) {
				mean[i*DIM+j] = sum[i*DIM+j]/count[i];
  			}
		}
		// printf("Mean distance took %.4lf seconds to run.\n", (float)(clock() - start) / CLOCKS_PER_SEC );
		// acc_mean += clock() - start;
	}
	for (i = 0; i < k; i++) {
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", mean[i*DIM+j]);
		printf("\n");
	}
	// printf("Acc distance %.4lf , Acc sum %.4lf Acc mean %.4lf.\n", 
	// 	(float)(acc_distance/CLOCKS_PER_SEC),
	// 	(float)(acc_sum/CLOCKS_PER_SEC), 
	// 	(float)(acc_mean/CLOCKS_PER_SEC)
	// );
	printf("All script took %.4lf seconds to run.\n", (float)(clock() - all_script) / CLOCKS_PER_SEC );
	#ifdef DEBUG
	for (i = 0; i < n; i++) {
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", x[i*DIM+j]);
		printf("%d\n", cluster[i]);
	}
	#endif
	return(0);
}
