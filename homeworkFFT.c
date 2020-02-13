#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <complex.h>
#include <math.h>
#include <time.h>
typedef double complex cd;
int N;
FILE* f1;
FILE* f2;
double complex* v;
int b;
int log2n;
double complex* out;
int P;
pthread_barrier_t barrier;
 

const double PI = 3.1415926536; 
  
unsigned int bitReverse(unsigned int x, int log2n) 
{ 
	int n = 0; 
    for (int i = 0; i < log2n; i++) 
    { 
	    n <<= 1; 
        n |= (x & 1); 
        x >>= 1; 
    } 
	return n; 
} 

void* threadFunction(void *var) 
{ 
    int n = N; 
	//int s;
	int thread_id = *(int*)var;
	if ((P==1) || (P == 0)) {
    for (int s = 1; s <= log2n; ++s) { 
        int m = 1 << s; 
        int m2 = m >> 1;
        cd w = 1; 
		cd wm = cos(PI/m2) + I*sin(PI/m2);  
		for (int j = 0; j < m2; ++j) { 
		
            for (int k = j; k < n; k += m) {  
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
            } 
            w *= wm; 
        } 
	}
	}
	
	if (P == 2) {
		//if (thread_id == 0) {
		for (int s = 1; s <= log2n; ++s) { 
		if (thread_id == 0) {
		
        int m = 1 << s; 
        int m2 = m >> 1;
        cd w = 1; 
		cd wm = cos(PI/m2) + I*sin(PI/m2);  
		for (int j = 0; j < m2; ++j) { 
			
            for (int k = j; k < n; k += m) {  
				if ((((m2 + k) %2 == 0) || (k%2==0)) || ((((m2 + k) %2 == 1) && (k%2==1))) ) {
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				}
            } 
            w *= wm; 
        } 
		} else {			
		
        int m = 1 << s; 
        int m2 = m >> 1;
        cd w = 1; 
		cd wm = cos(PI/m2) + I*sin(PI/m2);  
		for (int j = 0; j < m2; ++j) { 
            for (int k = j; k < n; k += m) {  
				if (!(((m2 + k) %2 == 0) || (k%2==0)) || ((((m2 + k) %2 == 1) && (k%2==1))) ) 
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				}
            } 
            w *= wm; 
        }
	}
	
	}
} 

	if (P == 4) {
		//if (thread_id == 0) {
		for (int s = 1; s <= log2n; ++s) { 
		if (thread_id % 4 == 0) {
		
        int m = 1 << s; 
        int m2 = m >> 1;
        cd w = 1; 
		cd wm = cos(PI/m2) + I*sin(PI/m2);  
		for (int j = 0; j < m2; ++j) { 
            for (int k = j; k < n; k += m) {  
				if (((m2 + k) %4 == thread_id) && (k%4==thread_id)) {
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				} 
				else {
				if (((m2 + k + 1) %4 == thread_id ) && ((k%4 + 1)==thread_id )) {
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				}
				else {
				if ((((m2 + k + 2)) %4 == thread_id ) && ((k%4 + 2 )==thread_id + 1)) {
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				}
				else {
				if (((m2 + k + 3) %4 == thread_id + 1) && ((k%4 + 3)==thread_id + 1)) {
                cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				}
				else {
					cd t = w * out[k + m2];  
                cd u = out[k]; 
                out[k] = u + t;  
                out[k + m2] = u - t;  
				}
				}
				}
				
				
				}
				} 
            w *= wm; 
        } 
		}
	
	
		}

	}

}
int main(int argc, char *argv[])
{
	double a;
	f1 = fopen(argv[1],"r");
	P = atoi(argv[3]);
	b = fscanf(f1,"%d", &N);
	v = malloc(sizeof(double complex) * N);
	out = malloc(sizeof(double complex) * N);
	
	for (int i = 0; i < N; i++) {
		b = fscanf(f1,"%lf", &a);
		v[i] = a;
	}
	fclose(f1);
	log2n = log2(N);
	for (unsigned int i = 0; i < N; ++i) {
        int rev = bitReverse(i, log2n); 
        out[i] = v[rev]; 
    } 
	clock_t begin = clock();

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

	
/* here, do your time-consuming job */

clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
///printf("%lf\n", time_spent);

	//printare in f2
	f2 = fopen(argv[2],"w");
	b = fprintf(f2,"%d\n", N);
	
	for (int i = 0; i < N; i++) {
		b = fprintf(f2,"%lf %lf\n", creal(out[i]), -cimag(out[i]));
	}
	fclose(f2);
	
	free(v);
	free(out);
	return 0;

}