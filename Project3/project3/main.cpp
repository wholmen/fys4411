#include <iostream>
#include <fstream>
#include <armadillo>
#include <WaveFunctions.h>
#include <VMCSolver.h>
#include <mpi.h>
#include <string>
#include <sstream>


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
    double Nparticles = 10; double alpha = 8.0; double alpha2; // Helium
    //double Nparticles = 4; double alpha = 3.5; // Beryllium

    int Ncycles = 1; int Nsteps = 1e3; // Ncycles: How many times to run MCIntegration for each each core. Nsteps: Iterations in MCIntegration
    NstepTotal = Nsteps*Nparticles;
    vec CycleE(Ncycles), CycleE2(Ncycles), TotalCycleE(Ncycles), TotalCycleE2(Ncycles), AllEnergies(NstepTotal); mat AllPositions(NstepTotal,Nparticles);

    // Set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();

    // Opening files for writing
    if (my_rank == 0) myfile.open("Neon.txt");
    ostringstream oss;
    oss << "Beryllium_my_rank_" << my_rank << ".txt";
    string var = oss.str();
    blocking_file.open(var.c_str());

    // Setting up the solver class
    WaveFunctions helium = WaveFunctions(alpha,beta_old,Nparticles);
    VMCSolver solve = VMCSolver(helium, my_rank);
    solve.AnalyticalEnergy = true;

    for (int nc=0; nc<Ncycles; nc++){
        // Steepest descent method to find optimal beta.
        beta_old = 0; beta_new = 0.6;
        Steepest_step = 0.1; // Around 0.1 is ok if beta_new is really close to optimal beta (~0.7)
        precision = 0.00001;
        alpha2 = 8.0;
        //alpha2 = alpha + nc*0.05;

        int n = 0; NstepSteepestDescent = 1000;
        while (abs(beta_new - beta_old) > precision && n < 1000){
            beta_old = beta_new;
            solve.MCintegration_FindVariables(NstepSteepestDescent,beta_old, alpha2);
            beta_new = beta_old - Steepest_step*solve.DiffEbeta;

            if (my_rank==1) cout << beta_new << endl;
            n++;
        }
        // Steepest descent method over. beta_new is the optimal beta.
        solve.FindStepLength(beta_new, alpha2);
        // Running the solver. Choose either between MCintegration or Importance sampling
        solve.MCintegration(CycleE, CycleE2, nc, AllEnergies, AllPositions, Ncycles, Nsteps, beta_new, alpha2);
        //solve.ImportanceSampling(CycleE, CycleE2, nc, AllEnergies, AllPositions, Ncycles, Nsteps, beta_new, alpha, 0.005);
    }

    for (int i=0; i<Ncycles; i++){
        MPI_Reduce(&CycleE(i), &TotalCycleE(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&CycleE2(i), &TotalCycleE2(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    time_stop = MPI_Wtime();
    total_time = time_stop - time_start;

    if (my_rank == 0) {
        for (int i=0; i<Ncycles; i++){
            //alpha2 = alpha + i * 0.05;
            alpha2 = 3.65;
            Etotal += TotalCycleE(i); Elocal = TotalCycleE(i) / numprocs;
            E2total += TotalCycleE2(i); E2local = TotalCycleE2(i) / numprocs;
            variance = E2local - Elocal*Elocal;
            myfile << "Mean value of all " << numprocs << " cores. Run number " << i << endl;
            myfile << "Alpha:    " << alpha2 << endl;
            myfile << "Beta:     " << beta_new << endl;
            myfile << "Step:     " << solve.step << endl;
            myfile << "Energy:   " << Elocal << endl;
            myfile << "Variance: " << variance << endl;
            myfile << "Error:    " << sqrt(variance) << endl << endl;
        }
        myfile << "Mean values for all runs." << endl;
        myfile << "Mean Energy:   " << Etotal / numprocs / Ncycles << endl;
        variance = E2total / numprocs / Ncycles - Etotal*Etotal / numprocs / numprocs / Ncycles / Ncycles;
        myfile << "Mean Variance: " << variance  << endl;
        myfile << "Mean Error:    " << sqrt(variance) << endl;
        myfile << "Total time needed: " << total_time;
        myfile.close();
    }

    blocking_file << Nsteps << " " << NstepTotal << endl;
    for (int i=0; i<NstepTotal; i++){
        blocking_file << AllEnergies(i) << " ";
        for (int p=0; p<Nparticles; p++) blocking_file << AllPositions(i,p) << " ";
        blocking_file << endl;
    }

    blocking_file.close();
    //cout << Etotal / 4.0 << endl;

    MPI_Finalize();

    return 0;
}


