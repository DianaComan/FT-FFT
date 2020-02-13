#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <complex.h>
#include <math.h>

int N;
FILE* f1;
FILE* f2;
double* v;
int b;
double complex* out;
int P;
pthread_barrier_t barrier;


void* threadFunction(void *var)
{
	int k;
	int thread_id = *(int*)var;
	if (thread_id <= N % P - 1) {
		for (k = thread_id * (N / P + 1); k < (thread_id + 1) * (N / P + 1); ++k) {
			complex double sum = 0.0;
			for (int t = 0; t < N; t++) {  // For each input element
			double angle = 2 * M_PI * t * k / N;
			sum += v[t] * cexp(-angle * I);
		}
		out[k] = sum;
		}
	}
	else {
		for (k = (N % P) + thread_id * (N / P); k < (N % P) + (thread_id + 1) * (N / P); ++k) {
			complex double sum = 0.0;
		for (int t = 0; t < N; t++) {  // For each input element
			double angle = 2 * M_PI * t * k / N;
			sum += v[t] * cexp(-angle * I);
		}
		out[k] = sum;
		}
	}
	/*for (int k = 0; k < N; k++) {  // For each output element
		complex double sum = 0.0;
		for (int t = 0; t < N; t++) {  // For each input element
			double angle = 2 * M_PI * t * k / N;
			sum += v[t] * cexp(-angle * I);
		}
		out[k] = sum;
	}*/
	
}

int main(int argc, char *argv[])
{
	//int i, j;
	//int P = 4;
	
	double a;
	f1 = fopen(argv[1],"r");
	P = atoi(argv[3]);
	b = fscanf(f1,"%d", &N);
	v = malloc(sizeof(double) * N);
	out = malloc(sizeof(double complex) * N);
	
	for (int i = 0; i < N; i++) {
		b = fscanf(f1,"%lf", &a);
		v[i] = a;
	}
	fclose(f1);

	pthread_barrier_init(&barrier, NULL, P);

	pthread_t tid[P];
	int thread_id[P];
	for(int i = 0;i < P; i++)
		thread_id[i] = i;

	for(int i = 0; i < P; i++) {
		pthread_create(&(tid[i]), NULL, threadFunction, &(thread_id[i]));
	}

	for(int i = 0; i < P; i++) {
		pthread_join(tid[i], NULL);
	}

	pthread_barrier_destroy(&barrier);

	//printare in f2
	f2 = fopen(argv[2],"w");
	b = fprintf(f2,"%d\n", N);
	
	for (int i = 0; i < N; i++) {
		b = fprintf(f2,"%lf %lf\n", creal(out[i]), cimag(out[i]));
	}
	fclose(f2);

	free(v);
	free(out);
	return 0;
}
