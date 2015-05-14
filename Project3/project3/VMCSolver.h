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
    double Nparticles;
    double Ndimensions;

public:
    double step;
    double Energy;
    double Variance;
    double AcceptRate;
    double h;
    double charge;
    bool AnalyticalEnergy;
    bool HFBasis;
    WaveFunctions atom;
    GaussianOrbitals atom2;
    void FindStepLength();

    VMCSolver(): idum(-1) {}
    VMCSolver(WaveFunctions a): idum(-1), step(2.45), AcceptRate(0), h(0.001), AnalyticalEnergy(false)
    {atom = a; charge = atom.Nparticles; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions; HFBasis = false;}

    VMCSolver(GaussianOrbitals a):  idum(-1), step(1.5), AcceptRate(0), h(0.001), AnalyticalEnergy(false)
    {atom2 = a; charge = atom2.Nparticles; Nparticles = atom2.Nparticles; Ndimensions = atom2.Ndimensions; HFBasis = true;}

    void MCintegration(int N);
    void ImportanceSampling(double DT, int Nstep);
    double LocalEnergy(mat R);
    double LocalEnergyAnalytical(mat R);

};


#endif // VMCSOLVER_H






