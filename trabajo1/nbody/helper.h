#pragma once
#define G 6.673e-11 //N m2/kg2

struct Pair2d
{
	double x, y;
	Pair2d() : x(.0f), y(.0f) {};
	Pair2d(double _x, double _y) : x(_x), y(_y) {};

	double norma()
	{
		return sqrt(this->x * this->x + this->y * this->y);
	}

	double norma2()
	{
		return this->x * this->x + this->y * this->y;
	}

	Pair2d operator+(const Pair2d &p2)
	{
		Pair2d out;
		out.x = this->x + p2.x;
		out.y = this->y + p2.y;
		return out;
	}

	Pair2d operator-(const Pair2d &p2)
	{
		Pair2d out;
		out.x = this->x - p2.x;
		out.y = this->y - p2.y;
		return out;
	}

	Pair2d operator*(const double &val)
	{
		Pair2d out;
		out.x = this->x * val;
		out.y = this->y * val;
		return out;
	}

	Pair2d operator/(const double &val)
	{
		Pair2d out;
		out.x = this->x / val;
		out.y = this->y / val;
		return out;
	}

	friend ostream& operator<<(ostream& os, const Pair2d& dt);
};

ostream& operator<<(ostream& os, const Pair2d& par)
{
    os << "[" << par.x << " , " << par.y << "]";
    return os;
}

struct MPI_BodiesInfo
{
	vector<double> posicion; // 2*N
	vector<double> velocidad; // 2*N
	vector<double> masa; // N
};

struct BodiesInfo
{
	vector<Pair2d> posicion;
	vector<Pair2d> velocidad;
	vector<double> masa;
};

struct MPI_Input
{
	int N; // N-bodies
	double dt;
	double TT;
	MPI_BodiesInfo mpi_bodiesInfo;
};

struct Input
{
	int N; // N-bodies
	double dt;
	double T;
	BodiesInfo bodiesInfo;
};



Input readInput(int argc, char** argv)
{
	Input in;

	if (argc < 4)
	{
		cout << "Use: parametros-> N dt T" << endl;
		cout << "- parametros: n bodies->" << endl;
		cout << "- N: Número de cuerpos a simular" << endl;
		cout << "- dt: time step" << endl;
		cout << "- T: time limit" << endl;
		exit(-1);
	}

	in.N = atoi(argv[1]);
	in.dt = atof(argv[2]);
	in.T = atof(argv[3]);
	in.bodiesInfo.masa.resize(in.N);
	in.bodiesInfo.posicion.resize(in.N);
	in.bodiesInfo.velocidad.resize(in.N);

	for (int i = 0; i < in.N; i++)
	{
		in.bodiesInfo.masa[i] = 100 + rand() % 200;
		in.bodiesInfo.posicion[i] = Pair2d(rand() % 10000, rand() % 10000);
		in.bodiesInfo.velocidad[i] = Pair2d(rand() % 100, rand() % 100);
	}

	return in;
}

MPI_Input read_MPI_Input(int argc, char** argv)
{
	MPI_Input in;

	if (argc < 4)
	{
		cout << "Use: parametros-> N dt T" << endl;
		cout << "- parametros: n bodies->" << endl;
		cout << "- N: Número de cuerpos a simular" << endl;
		cout << "- dt: time step" << endl;
		cout << "- T: time limit" << endl;
		exit(-1);
	}

	in.N = atoi(argv[1]);
	in.dt = atof(argv[2]);
	in.TT = atof(argv[3]);
	in.mpi_bodiesInfo.masa.resize(in.N);
	in.mpi_bodiesInfo.posicion.resize(2 * in.N);
	in.mpi_bodiesInfo.velocidad.resize(2 * in.N);

	for (int i = 0; i < in.N; i++)
	{
		in.mpi_bodiesInfo.masa[i] = 100 + rand() % 200;

		in.mpi_bodiesInfo.posicion[2 * i] = rand() % 10000;
		in.mpi_bodiesInfo.posicion[2 * i + 1] = rand() % 10000;

		in.mpi_bodiesInfo.velocidad[2 * i] = rand() % 100;
		in.mpi_bodiesInfo.velocidad[2 * i + 1] = rand() % 100;
	}

	return in;
}

Input readInput(int N)
{
	Input in;

	in.N = N;
	in.dt = 0.01;
	in.T = 100;
	in.bodiesInfo.masa.resize(N);
	in.bodiesInfo.posicion.resize(N);
	in.bodiesInfo.velocidad.resize(N);

	for (int i = 0; i < N; i++)
	{
		in.bodiesInfo.masa[i] = 100 + rand() % 200;
		in.bodiesInfo.posicion[i] = Pair2d(rand() % 10000, rand() % 10000);
		in.bodiesInfo.velocidad[i] = Pair2d(rand() % 100, rand() % 100);
	}

	return in;
}
void print(vector<Pair2d> posicion, vector<Pair2d> velocidad, vector<Pair2d> fuerza, vector<double> masa)
{
	int N;
	N = posicion.size();
	for (int i = 0; i < N; i++)
	{
		cout << "----Cuerpo " << (i + 1) << ":" << endl;
		cout << "masa: " << masa[i] << endl;
		cout << "posicion: " << posicion[i] << endl;
		cout << "velocidad: " << velocidad[i] << endl;
		cout << "fuerza: " << fuerza[i] << endl;
	}
}




