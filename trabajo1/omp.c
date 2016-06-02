#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
void mymethod(int it) {
	printf("%d : %d of %d\n", it, omp_get_thread_num(), omp_get_num_threads());
}
int main() {
	
	int thread_count = 8;
	#pragma omp parallel for num_threads(thread_count)
	for (int i = 0; i < 7; i++)
	mymethod(i);
	return 0;
}