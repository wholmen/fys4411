#include <iostream>
#include <fstream>
#include <armadillo>
#include <WaveFunctions.h>
#include <VMCSolver.h>
#include <string>
#include <sstream>
#include <mpi.h>

using namespace std;
using namespace arma;

ofstream myfile, blocking_file;

int numprocs, my_rank, NstepTotal, NstepSteepestDescent;
double Etotal(0), Elocal(0), E2total(0), E2local(0), variance(0);
double beta_new, beta_old, Steepest_step, precision;
double time_start, time_stop, total_time;


int main(int argc, char* argv[])
{
    // Define what particle we are solving
    //double Nparticles = 10;// Neon
    //double Nparticles = 4; // Beryllium
    double Nparticles = 10;

    int Ncycles = 1; int Nsteps = 2e3; // Ncycles: How many times to run MCIntegration for each each core. Nsteps: Iterations in MCIntegration
    NstepTotal = Nsteps*Nparticles;
    vec CycleE(Ncycles), CycleE2(Ncycles), TotalCycleE(Ncycles), TotalCycleE2(Ncycles), AllEnergies(NstepTotal); mat AllPositions(NstepTotal,Nparticles);

    // Set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();

    // Opening files for writing
    if (my_rank == 0) myfile.open("Neon_Final.txt");
    ostringstream oss;
    oss << "Neon_my_rank_" << my_rank << ".txt";
    string var = oss.str();
    blocking_file.open(var.c_str());

    // Setting up the solver class
    WaveFunctions atom = WaveFunctions(beta_old,Nparticles);
    VMCSolver solve = VMCSolver(atom, my_rank);
    solve.AnalyticalEnergy = false;

    for (int nc=0; nc<Ncycles; nc++){
        // Steepest descent method to find optimal beta.
        beta_old = 0; beta_new = 0.6;
        Steepest_step = 0.2; // Around 0.1 is ok if beta_new is really close to optimal beta (~0.7)
        precision = 0.0001;

        int n = 0; NstepSteepestDescent = 100;
        while (abs(beta_new - beta_old) > precision && n < 1000){
            beta_old = beta_new;

            solve.MCintegration_FindVariables(NstepSteepestDescent,beta_old);

            beta_new = beta_old - Steepest_step*solve.DiffEbeta;

            if (my_rank==1) cout << "Beta: " << beta_new << endl;
            n++;
        }
        // Steepest descent method over. beta_new is the optimal beta.
        //solve.FindStepLength(beta_new);
        solve.step = 1.3;
        // Running the solver. Choose either between MCintegration or Importance sampling
        solve.MCintegration(CycleE, CycleE2, nc, AllEnergies, AllPositions, Ncycles, Nsteps, beta_new);
        //solve.ImportanceSampling(CycleE, CycleE2, nc, AllEnergies, AllPositions, Ncycles, Nsteps, beta_new, 0.005);
    }

    for (int i=0; i<Ncycles; i++){
        MPI_Reduce(&CycleE(i), &TotalCycleE(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&CycleE2(i), &TotalCycleE2(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    time_stop = MPI_Wtime();
    total_time = time_stop - time_start;

    if (my_rank == 0) {
        for (int i=0; i<Ncycles; i++){
            Etotal += TotalCycleE(i); Elocal = TotalCycleE(i) / numprocs;
            E2total += TotalCycleE2(i); E2local = TotalCycleE2(i) / numprocs;
            variance = E2local - Elocal*Elocal;
            myfile << "Mean value of all " << numprocs << " cores. Run number " << i << endl;
            myfile << "Beta:     " << beta_new << endl;
            //myfile << "Step:     " << solve.step << endl;
            myfile << "Energy:   " << Elocal << endl;
            myfile << "Variance: " << variance << endl;
            myfile << "Error:    " << sqrt(variance) << endl << endl;
            //myfile << "Total time needed: " << total_time;
        }
        myfile << "Mean values for all runs." << endl;
        myfile << "Mean Energy:   " << Etotal / numprocs / Ncycles << endl;
        variance = E2total / numprocs / Ncycles - Etotal*Etotal / numprocs / numprocs / Ncycles / Ncycles;
        myfile << "Mean Variance: " << variance  << endl;
        myfile << "Mean Error:    " << sqrt(variance) << endl;
        myfile << "Total time needed: " << total_time;
        myfile.close();
    }

    //blocking_file << Nsteps << " " << NstepTotal << endl;
    for (int i=0; i<NstepTotal; i++){
        blocking_file << AllEnergies(i) << " ";
        for (int p=0; p<Nparticles; p++) blocking_file << AllPositions(i,p) << " ";
        blocking_file << endl;
    }

    blocking_file.close();

    MPI_Finalize();
}

