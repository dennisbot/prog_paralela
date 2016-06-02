#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <pthread.h>

void test (void* rank);

int main(int argc, char* argv[]) { 
	srand(time(0));
	time_t t = time(0);
	
	struct tm * timeinfo = localtime(&t);
	
	long num = 100 + rand() % 200;
	printf ( "Current local time and date: %s", asctime (timeinfo) );
	printf("%ld y puntero es : %p\n", num, &num);
	test((void*) num);
	return 0; 
} 
/* main */
void test(void* rank) { 
	long num = (long) rank;
	/* Use long in case of 64−bit system */ 
	printf("Hello from thread %ld of %p\n", num, rank);
} 
/* test */