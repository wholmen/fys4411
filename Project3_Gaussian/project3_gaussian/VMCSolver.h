#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <iostream>
#include <WaveFunctions.h>
#include <lib.h>

using namespace std;
using namespace arma;

// This is the integrator class. The time steps and metropolis algorithm is in this class.

class VMCSolver
{
private:
    long idum;
    double Nparticles, Ndimensions;

public:
    double step, Energy, Variance, AcceptRate, h, charge, DiffEbeta;
    bool AnalyticalEnergy;
    int my_rank;
    WaveFunctions atom;

    VMCSolver(): idum(-1) {}
    VMCSolver(WaveFunctions a, int my_rank1): step(2.45), AcceptRate(0), h(0.001), AnalyticalEnergy(false), my_rank(my_rank1)
    {atom = a; charge = atom.Nparticles; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions; idum = -1-my_rank;}

    // Functions to find parametres
    void FindStepLength(double b);
    void MCintegration_FindVariables(int Nsteps, double beta);

    // Functions to calculate
    void MCintegration(vec &Cycle_Energy, vec &Cycle_Energy2, int ncycle, vec &AllEnergies, mat &AllPositions, int Ncycles, int Niterations, double b);
    void ImportanceSampling(vec &Cycle_Energy, vec &Cycle_Energy2, int ncycle, vec &AllEnergies, mat &AllPositions, int Ncycles, int Niterations, double b, double DT);

    // Functions to calculate energy
    double LocalEnergy(mat R);
    double LocalEnergyAnalytical(mat R);

};



#endif // VMCSOLVER_H
