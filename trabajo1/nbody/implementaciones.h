#pragma once
#include "helper.h"

void nBodySecuenciaV1(const Input &in)
{
	vector<Pair2d> fuerza, posicion, velocidad;
	vector<double> masa;
	double T_end, dt;
	int steps, N, q;
	Pair2d p_diff;
	double distancia, distancia_3;
	N = in.N;
	dt = in.dt;
	T_end = in.T;
	fuerza.resize(N);
	posicion = in.bodiesInfo.posicion;
	velocidad = in.bodiesInfo.velocidad;
	masa = in.bodiesInfo.masa;

#ifdef DEBUG
	cout << "------------Condiciones iniciales:" << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG

	steps = (int)(T_end / dt);

	for (int i = 0; i < steps; i++)
	{
		//para cada particula q, calculamos la fuerza total sobre q
		for (int q = 0; q < N; q++)
		{
			for (int k = 0; k < q; k++)
			{
				p_diff = posicion[k] - posicion[q];
				distancia = p_diff.norma();
				distancia_3 = distancia * distancia * distancia;
				fuerza[q] = fuerza[q] - p_diff * ((G * masa[k] * masa[q]) / distancia_3);
			}

			for (int k = q + 1; k < N; k++)
			{
				p_diff = posicion[k] - posicion[q];
				distancia = p_diff.norma();
				distancia_3 = distancia * distancia * distancia;
				fuerza[q] = fuerza[q] - p_diff * ((G * masa[k] * masa[q]) / distancia_3);
			}
		}
		
		// para cada particula q calcular la posición y velocidad de q
		for (int q = 0; q < N; q++)
		{
			posicion[q] = posicion[q] + velocidad[q] * dt;
			velocidad[q] = velocidad[q] + fuerza[q] * (dt / masa[q]);
		}
	}

#ifdef DEBUG
	cout << "------------Resume:" << endl;
	cout << "Steps = " << steps << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG
}
void nBodySecuenciaV2(const Input &in)
{
	vector<Pair2d> fuerza, posicion, velocidad;
	vector<double> masa;
	double T_end, dt;
	int steps, N, q;
	Pair2d f, p_diff;
	double distancia, distancia_3;

	N = in.N;
	dt = in.dt;
	T_end = in.T;
	fuerza.resize(N);
	posicion = in.bodiesInfo.posicion;
	velocidad = in.bodiesInfo.velocidad;
	masa = in.bodiesInfo.masa;

#ifdef DEBUG
	cout << "------------Condiciones iniciales:" << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG

	steps = (int)(T_end / dt);

	for (int i = 0; i < steps; i++)
	{
		//for each particle q, Compute total fuerza on q
		for (int q = 0; q < N; q++)
		{
			for (int k = q + 1; k < N; k++)
			{
				p_diff = posicion[q] - posicion[k];
				distancia = p_diff.norma();
				distancia_3 = distancia * distancia * distancia;
		
				f = p_diff * ((G * masa[q] * masa[k]) / distancia_3);
				fuerza[q] = fuerza[q] + f;
				fuerza[k] = fuerza[k] - f;
			}
		}
		
		//for each particle q Compute posicion and velocidad of q
		for (int q = 0; q < N; q++)
		{
			posicion[q] = posicion[q] + velocidad[q] * dt;
			velocidad[q] = velocidad[q] + fuerza[q] * (dt / masa[q]);
		}
	}

#ifdef DEBUG
	cout << "------------Resume:" << endl;
	cout << "Steps = " << steps << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG
}
void NBodyOMP_V1(const Input &in)
{
	vector<Pair2d> fuerza, posicion, velocidad;
	vector<double> masa;
	double T_end, dt;
	int threads, steps, N, q;
	Pair2d p_diff;
	double distancia, distancia_3;

	N = in.N;
	dt = in.dt;
	T_end = in.T;
	fuerza.resize(N);
	posicion = in.bodiesInfo.posicion;
	velocidad = in.bodiesInfo.velocidad;
	masa = in.bodiesInfo.masa;

#ifdef DEBUG
	cout << "------------Condiciones iniciales:" << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG

	steps = (int)(T_end / dt);
	#pragma omp parallel
	{
		threads = omp_get_num_threads();

		for (int i = 0; i < steps; i++)
		{
			//for each particle q, Compute total fuerza on q
			#pragma omp parallel for schedule(static, N/threads) private(p_diff, distancia, distancia_3)
			for (int q = 0; q < N; q++)
			{
				for (int k = 0; k < q; k++)
				{
					p_diff = posicion[k] - posicion[q];
					distancia = p_diff.norma();
					distancia_3 = distancia * distancia * distancia;
					fuerza[q] = fuerza[q] - p_diff * ((G * masa[k] * masa[q]) / distancia_3);
				}

				for (int k = q + 1; k < N; k++)
				{
					p_diff = posicion[k] - posicion[q];
					distancia = p_diff.norma();
					distancia_3 = distancia * distancia * distancia;
					fuerza[q] = fuerza[q] - p_diff * ((G * masa[k] * masa[q]) / distancia_3);
				}
			}
		
			//for each particle q Compute posicion and velocidad of q
			#pragma omp parallel for schedule(static, N/threads)
			for (int q = 0; q < N; q++)
			{
				posicion[q] = posicion[q] + velocidad[q] * dt;
				velocidad[q] = velocidad[q] + fuerza[q] * (dt / masa[q]);
			}
		}
	}

#ifdef DEBUG
	cout << "------------Resume:" << endl;
	cout << "Steps = " << steps << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG
}
void NBodyOMP_V2(const Input &in)
{
	vector<Pair2d> fuerza, posicion, velocidad;
	vector<double> masa;
	double T_end, dt;
	int threads, steps, N, q;
	Pair2d f, p_diff;
	double distancia, distancia_3;

	N = in.N;
	dt = in.dt;
	T_end = in.T;
	fuerza.resize(N);
	posicion = in.bodiesInfo.posicion;
	velocidad = in.bodiesInfo.velocidad;
	masa = in.bodiesInfo.masa;

	vector<omp_lock_t> my_locks;
	my_locks.resize(N);
	for (int i = 0; i < N; i++)
		omp_init_lock(&(my_locks[i]));


#ifdef DEBUG
	cout << "------------Condiciones iniciales:" << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG

	steps = (int)(T_end / dt);

	#pragma omp parallel
	{
		threads = omp_get_num_threads();

		for (int i = 0; i < steps; i++)
		{
			//for each particle q, Compute total fuerza on q
			#pragma omp parallel for schedule(static, N/threads) private(p_diff, distancia, distancia_3, f)
			for (int q = 0; q < N; q++)
			{
				for (int k = q + 1; k < N; k++)
				{
					p_diff = posicion[q] - posicion[k];
					distancia = p_diff.norma();
					distancia_3 = distancia * distancia * distancia;
		
					f = p_diff * ((G * masa[q] * masa[k]) / distancia_3);
				
					omp_set_lock(&(my_locks[q]));
					fuerza[q] = fuerza[q] + f;
					omp_unset_lock(&(my_locks[q]));

					omp_set_lock(&(my_locks[k]));
					fuerza[k] = fuerza[k] - f;
					omp_unset_lock(&(my_locks[k]));
				}
			}
		
			//for each particle q Compute posicion and velocidad of q
			#pragma omp parallel for
			for (int q = 0; q < N; q++)
			{
				posicion[q] = posicion[q] + velocidad[q] * dt;
				velocidad[q] = velocidad[q] + fuerza[q] * (dt / masa[q]);
			}
		}
	}

#ifdef DEBUG
	cout << "------------Resume:" << endl;
	cout << "Steps = " << steps << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG
}
void NBodyOMP_V3(const Input &in)
{
	vector<Pair2d> *loc_fuerzas, fuerza, posicion, velocidad;
	vector<double> masa;
	double T_end, dt;
	int thread_id, threads, steps, N, q;
	Pair2d f, p_diff;
	double distancia, distancia_3;

	N = in.N;
	dt = in.dt;
	T_end = in.T;
	fuerza.resize(N);
	posicion = in.bodiesInfo.posicion;
	velocidad = in.bodiesInfo.velocidad;
	masa = in.bodiesInfo.masa;

	vector<omp_lock_t> my_locks;
	my_locks.resize(N);
	for (int i = 0; i < N; i++)
		omp_init_lock(&(my_locks[i]));


#ifdef DEBUG
	cout << "------------Condiciones iniciales:" << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG

	steps = (int)(T_end / dt);

	#pragma omp parallel private(thread_id)
	{
		#pragma omp master
		{
			threads = omp_get_num_threads();
			loc_fuerzas = new vector<Pair2d>[threads];
			for (int i = 0; i < threads; i++)
				loc_fuerzas[i].resize(N);
		}

		thread_id = omp_get_thread_num();

		#pragma omp barrier

		for (int i = 0; i < steps; i++)
		{
			//for each particle q, Compute total fuerza on q
			//remove static schedule 
			#pragma omp parallel for schedule(static, N/threads) private(p_diff, distancia, distancia_3, f)
			for (int q = 0; q < N; q++)
			{
				for (int k = q + 1; k < N; k++)
				{
					p_diff = posicion[q] - posicion[k];
					distancia = p_diff.norma();
					distancia_3 = distancia * distancia * distancia;
		
					f = p_diff * ((G * masa[q] * masa[k]) / distancia_3);
				
					//omp_set_lock(&(my_locks[q]));
					//fuerza[q] = fuerza[q] + f;
					//omp_unset_lock(&(my_locks[q]));

					//omp_set_lock(&(my_locks[k]));
					//fuerza[k] = fuerza[k] - f;
					//omp_unset_lock(&(my_locks[k]));
					loc_fuerzas[thread_id][q] = loc_fuerzas[thread_id][q] + f;
					loc_fuerzas[thread_id][k] = loc_fuerzas[thread_id][q] - f;
				}
			}

			#pragma omp parallel for
			for (int q = 0; q < N; q++)
			{
				for (int id = 0; id < threads; id++)
				{
					fuerza[q] = fuerza[q] + loc_fuerzas[thread_id][q];
				}
			}
			
		
			//for each particle q Compute posicion and velocidad of q
			#pragma omp parallel for
			for (int q = 0; q < N; q++)
			{
				posicion[q] = posicion[q] + velocidad[q] * dt;
				velocidad[q] = velocidad[q] + fuerza[q] * (dt / masa[q]);
			}
		}
	}

#ifdef DEBUG
	cout << "------------Resume:" << endl;
	cout << "Steps = " << steps << endl;
	print(posicion, velocidad, fuerza, masa);
#endif // DEBUG
}
// void n_body_mpi_v1(const MPI_Input &in)
// {
// 	vector<T>
// 		l_fuerza,
// 		masa,
// 		posicion,
// 		g_velocidad,
// 		l_velocidad;
// 	T
// 		T_end,
// 		dt,
// 		p_diff[2],
// 		distancia,
// 		distancia_3;
// 	int
// 		steps,
// 		l_N,
// 		g_N,
// 		q;
// 	double
// 		in_buf[3];

// 	if (id == 0)
// 	{
// 		in_buf[0] = in.N;
// 		in_buf[1] = in.dt;
// 		in_buf[2] = in.T;

// 		masa = in.mpi_bodiesInfo.masa;
// 		posicion = in.mpi_bodiesInfo.posicion;

// 		g_velocidad = in.mpi_bodiesInfo.velocidad;
// 	}

// 	MPI_Bcast(&in_buf[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	g_N = (int) (in_buf[0]);
// 	dt = in_buf[1];
// 	T_end = in_buf[2];

// 	masa.resize(g_N);
// 	posicion.resize(2 * g_N);

// 	l_N = g_N / nTasks;

// 	if (id == 0)
// 	{
// 		cout << "g_N = " << g_N << endl;
// 		cout << "l_N = " << l_N << endl;
// 		cout << "dt = " << dt << endl;
// 		cout << "T_end = " << T_end << endl;
// 		cout << "masa.size() = " << masa.size() << endl;
// 		cout << "posicion.size() = " << posicion.size() << endl;
// 	}

// 	MPI_Bcast(	&masa[0], 
// 				g_N, 
// 				MPI_DOUBLE, 
// 				0, 
// 				MPI_COMM_WORLD);
// 	MPI_Bcast(	&posicion[0], 
// 				2 * g_N, 
// 				MPI_DOUBLE, 
// 				0, 
// 				MPI_COMM_WORLD);

	

// 	l_velocidad.resize(2 * l_N);
// 	MPI_Scatter(&g_velocidad[0],
// 				2 * l_N,
// 				MPI_DOUBLE,
// 				&l_velocidad[0],
// 				2 * l_N,
// 				MPI_DOUBLE,
// 				0,
// 				MPI_COMM_WORLD);

// 	l_fuerza.resize(2 * l_N);

// #ifdef DEBUG
// 	cout << "------------Condiciones iniciales:" << endl;
// 	//print(posicion, velocidad, fuerza, masa);
// #endif // DEBUG

// 	steps = (int)(T_end / dt);

// 	for (int i = 0; i < steps; i++)
// 	{
// 		//for each particle q, Compute total fuerza on q
// 		for (int l_q = id*l_N; l_q < (id + 1)*l_N; l_q++)
// 		{
// 			for (int k = 0; k < l_q; k++)
// 			{
// 				p_diff[0] = posicion[2*k] - posicion[2*l_q];
// 				p_diff[1] = posicion[2*k+1] - posicion[2*l_q+1];
// 				distancia = sqrt(p_diff[0]*p_diff[0] + p_diff[1]*p_diff[1]);
// 				distancia_3 = distancia * distancia * distancia;
// 				l_fuerza[2*(l_q % l_N)] -= p_diff[0] * ((G * masa[k] * masa[l_q]) / distancia_3);
// 				l_fuerza[2*(l_q % l_N) + 1] -= p_diff[1] * ((G * masa[k] * masa[l_q]) / distancia_3);
// 			}

// 			for (int k = l_q + 1; k < g_N; k++)
// 			{
// 				p_diff[0] = posicion[2*k] - posicion[2*l_q];
// 				p_diff[1] = posicion[2*k+1] - posicion[2*l_q+1];
// 				distancia = sqrt(p_diff[0]*p_diff[0] + p_diff[1]*p_diff[1]);
// 				distancia_3 = distancia * distancia * distancia;
// 				l_fuerza[2*(l_q % l_N)] -= p_diff[0] * ((G * masa[k] * masa[l_q]) / distancia_3);
// 				l_fuerza[2*(l_q % l_N) + 1] -= p_diff[1] * ((G * masa[k] * masa[l_q]) / distancia_3);
// 			}
// 		}
		
// 		//for each particle q Compute posicion and velocidad of q
// 		for (int l_q = id*l_N; l_q < (id + 1)*l_N; l_q++)
// 		{
// 			posicion[2*l_q] += l_velocidad[2*(l_q % l_N)] * dt;
// 			posicion[2*l_q] += l_velocidad[2*(l_q % l_N) + 1] * dt;

// 			l_velocidad[2*(l_q % l_N)] += l_fuerza[2*(l_q % l_N)] * (dt / masa[l_q]);
// 			l_velocidad[2*(l_q % l_N) + 1] += l_fuerza[2*(l_q % l_N) + 1] * (dt / masa[l_q]);
// 		}

// 		MPI_Allgather(	&(posicion[2*id*l_N]),
// 						2 * l_N, 
// 						MPI_DOUBLE, 
// 						&posicion[0], 
// 						2 * l_N, 
// 						MPI_DOUBLE, 
// 						MPI_COMM_WORLD); 
// 	}

// 	MPI_Gather(	&(l_velocidad[2*id*l_N]),
// 					2 * l_N, 
// 					MPI_DOUBLE, 
// 					&g_velocidad[0], 
// 					2 * l_N, 
// 					MPI_DOUBLE,
// 					0,
// 					MPI_COMM_WORLD); 

// 	if (id == 0)
// 	{
// #ifdef DEBUG
// 		cout << "------------Resume:" << endl;
// 		cout << "Steps = " << steps << endl;
// 		//MPI_print(posicion, g_velocidad);
// #endif // DEBUG		
// 	}
// }
// void n_body_mpi_v2(const MPI_Input &in)
// {
// 	cout << "Pendiente!" << endl;
// }
