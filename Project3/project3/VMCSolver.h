#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <iostream>
#include <WaveFunctions.h>
#include <lib.h>
#include <GaussianOrbitals.h>

using namespace std;
using namespace arma;

class VMCSolver
{

private:
    long idum;
    double Nparticles, Ndimensions;

public:
    double step, Energy, Variance, AcceptRate, h, charge, DiffEbeta;
    bool AnalyticalEnergy, HFBasis;
    int my_rank;

    WaveFunctions atom;
    GaussianOrbitals atom2;
    void FindStepLength(double b, double a);

    VMCSolver(): idum(-1) {}
    VMCSolver(WaveFunctions a, int myrank): step(2.45), AcceptRate(0), h(0.001), AnalyticalEnergy(false), my_rank(myrank)
    {atom = a; charge = atom.Nparticles; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions; HFBasis = false; idum = -1-my_rank;}

    VMCSolver(GaussianOrbitals a):  idum(-1), step(1.5), AcceptRate(0), h(0.001), AnalyticalEnergy(false)
    {atom2 = a; charge = atom2.Nparticles; Nparticles = atom2.Nparticles; Ndimensions = atom2.Ndimensions; HFBasis = true;}

    void MCintegration(vec &Cycle_Energy, vec &Cycle_Energy2, int ncycle, vec &AllEnergies, mat &AllPositions, int Ncycles, int Niterations, double b, double a);
    void MCintegration_FindVariables(int Nsteps, double beta, double a);
    void ImportanceSampling(vec &Cycle_Energy, vec &Cycle_Energy2, int ncycle, vec &AllEnergies, mat &AllPositions, int Ncycles, int Niterations, double b, double a, double DT);
    double LocalEnergy(mat R);
    double LocalEnergyAnalytical(mat R);

};


#endif // VMCSOLVER_H






