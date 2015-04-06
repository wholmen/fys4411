#include <iostream>
#include <armadillo>
#include <WaveFunctions.h>
#include <VMCSolver.h>
#include <mpi.h>

using namespace arma;
using namespace std;




int main(int nargs, char* args[])
{

    WaveFunctions helium = WaveFunctions(0.5,1,4);

    mat r = zeros<mat>(4,3);
    cout << helium.WaveFunction(r);


    VMCSolver solve = VMCSolver(helium);
    int numprocs, my_rank, N, Nlocal;
    double Etotal, Elocal;

    Etotal = 0;
    N = 10000;
    Nlocal = N/4;
    MPI_Init(&nargs, &args);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    solve.MCintegration(Nlocal);

    Elocal = solve.Energy;
    Etotal += Elocal;

    cout << "Processor with rank " << my_rank << " calculated energy " << Elocal << endl;

    MPI_Finalize();

    cout << Etotal / 4.0;
    //cout << endl << solve.LocalEnergy(r) << endl;

    //cout << solve.Energy;
}


