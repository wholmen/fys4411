#include <iostream>
#include <fstream>
#include <armadillo>
#include <WaveFunctions.h>
#include <VMCSolver.h>
#include <string>
#include <sstream>


using namespace std;
using namespace arma;

int main()
{
    vec Cycle_Energy, Cycle_Energy2, AllEnergies;
    int ncycle;
    mat AllPositions;
    int Ncycles, Niterations;
    double beta;

    Ncycles = 1; ncycle = 1; Niterations = 1000000; beta = 0.9;
    Cycle_Energy = zeros<vec>(Ncycles); Cycle_Energy2 = zeros<vec>(Ncycles); AllEnergies = zeros<vec>(Ncycles);

    WaveFunctions helium = WaveFunctions(beta,2);
    VMCSolver solve = VMCSolver(helium,0);
    solve.AnalyticalEnergy = true;
    solve.FindStepLength(0.5);
    //solve.MCintegration(Cycle_Energy, Cycle_Energy2, ncycle, AllEnergies, AllPositions, Ncycles, Niterations, beta);
    solve.ImportanceSampling(Cycle_Energy, Cycle_Energy2, ncycle, AllEnergies, AllPositions, Ncycles, Niterations, beta, 0.005);
    cout << "Ground state energy: " << solve.Energy << endl;

    helium = WaveFunctions(beta,2);
    solve = VMCSolver(helium,0);
    solve.AnalyticalEnergy = true;
    //solve.FindStepLength(0.5);
    //solve.MCintegration(Cycle_Energy, Cycle_Energy2, ncycle, AllEnergies, AllPositions, Ncycles, Niterations, beta);
    //solve.ImportanceSampling(Cycle_Energy, Cycle_Energy2, ncycle, AllEnergies, AllPositions, Ncycles, Niterations, beta, 0.005);
    //cout << "Ground state energy: " << solve.Energy << endl;
    //GaussianOrbitals particle = GaussianOrbitals(1,0,0,0,2);
    //cout << particle.psi(0.2) << endl;
    //cout << exp(-1.6*0.2);
    //cout << 0;
}

