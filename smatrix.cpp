#include <bits/stdc++.h>
#define N 10000
using namespace std;

int i, j, k;
int mul1[N][N], mul2[N][N], res[N][N], temp[N][N], sum;

void resetMatrix() {
	sum = 0;
	memset(res, sizeof(int), 0);
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			res[i][j] = 0;
		}
	}
}
void printMatrix(int matrix[N][N], string message) {
	printf("%s\n", message.c_str());
	 for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			printf("%d\t", matrix[i][j]);
		}
		puts("");
	}
	puts("");
}
// three nested loop version (with a acumulator to verify that we are working with the same
// matrix when comparing with other implementations)
void basicWay() {
	resetMatrix();
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j) {
			for (k = 0; k < N; ++k)
				res[i][j] += mul1[i][k] * mul2[k][j];
			sum += res[i][j];
		}
}
void transposeWay() {
	resetMatrix();
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
				temp[j][i] = mul2[i][j];
	
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j) {
			for (k = 0; k < N; ++k) {
				res[i][j] += mul1[i][k] * temp[j][k];
			}
			sum += res[i][j];
		}
} 
void performOne() {
	sum = 0;
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			res[i][j] = mul1[i][j] * mul2[i][j], sum += res[i][j];
}

void performTwo() {
	sum = 0;
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j)
			res[j][i] = mul1[j][i] * mul2[j][i], sum += res[j][i];
}

int main() {
	
	// memset(mul1, N / sizeof(int), 0);
	// memset(mul2, N / sizeof(int), 0);
	srand(time(0));
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			mul1[i][j] = 1;
			mul2[i][j] = rand() % 10;
		}
	}
	mul1[0][0] = 1; mul1[0][1] = 2;
	mul1[1][0] = 3; mul1[1][1] = 4;

	mul2[0][0] = 4; mul2[0][1] = 3;
	mul2[1][0] = 2; mul2[1][1] = 1;
	
	// printMatrix(mul1);
	// printMatrix(mul2);

	// clock_t t = clock();
	// basicWay();
	// printf("basicWay time: %.3f\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// printf("sum is: %d\n", sum);
	// t = clock();
	// transposeWay();
	// printf("transposeWay time: %.4f\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// printf("sum is: %d\n", sum);
	
	// perform one and two
	clock_t t = clock();
	performOne();
	printf("performOne time: %.2f\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("SUM: %d\n", sum);
	t = clock();
	performTwo();
	printf("performOne time: %.2f\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("SUM: %d\n", sum);

	
	// for (i = 0; i < N; i++) {
	// 	for (j = 0; j < N; j++) {
	// 		printf("%d\t", res[i][j]);
	// 	}
	// 	puts("");
	// }
	return 0;
}