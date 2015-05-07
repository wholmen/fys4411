#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <iostream>
#include <WaveFunctions.h>
#include <lib.h>

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
    WaveFunctions atom;

    VMCSolver(): idum(-1) {}
    VMCSolver(WaveFunctions a): idum(-1), step(1.5), AcceptRate(0), h(0.001), AnalyticalEnergy(false)
    {atom = a; charge = atom.Nparticles; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions;}

    void MCintegration(int N);
    double LocalEnergy(mat R);
    double LocalEnergyAnalytical(mat R);
};


#endif // VMCSOLVER_H
