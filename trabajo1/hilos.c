#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

/* Global variable: accessible to all threads */
int thread_count;

void *Hello(void* rank); /* Thread function */

int main(int argc, char* argv[]) { 
	long thread;
	pthread_t* thread_handles;
	srand(time(0));
	
	/* Get number of threads from command line */ 
	thread_count = strtol(argv[1], NULL, 10);

	thread_handles = malloc (thread_count*sizeof(pthread_t));

	for (thread = 0; thread < thread_count; thread++)
	pthread_create(&thread_handles[thread], NULL, 
		Hello, (void*) thread);
	
	printf("Hello from the main thread\n");
	long* val[thread_count];
	for (thread = 0; thread < thread_count; thread++) 
		pthread_join(thread_handles[thread], (void*)&val[thread]);

	puts("valores aleatorios");
	for (int i = 0; i < thread_count; i++) {
		printf("%ld\n", val[i]);
	}

	free(thread_handles);
	return 0; 
} /* main*/

void* Hello(void* rank) { 
	long my_rank = (long) rank;
		/* Use long in case of 64−bit system */ 
	printf("Hello from thread %ld of %d\n", my_rank, 
		thread_count);
	return (void*) ((long) rand() % 200); 
}  /*Hello*/