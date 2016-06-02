#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cmath>
// #include <mpi.h>

using namespace std;

#define DEBUG

//global variables
int nTasks, id;
//
#include "helper.h"
#include "implementaciones.h"

int main(int argc, char** argv)
{
	Input in;
	MPI_Input mpi_in;

	in = readInput(argc, argv);
	// nBodySecuenciaV1(in);
	// nBodySecuenciaV2(in);
	// NBodyOMP_V1(in);
	NBodyOMP_V2(in);
	// NBodyOMP_V3(in);

	//MPI_Init(&argc,&argv);

	//MPI_Comm_size(MPI_COMM_WORLD,&nTasks);
	//MPI_Comm_rank(MPI_COMM_WORLD,&id);
	//
	//if (id == 0)
	//	mpi_in = read_MPI_Input(argc, argv);

	//n_body_mpi_v1(mpi_in);
	//n_body_mpi_v2(mpi_in);

	//MPI_Finalize();

	puts("terminado ...");
	return 0;
}
