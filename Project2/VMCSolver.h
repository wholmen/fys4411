#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <iostream>
#include <WaveFunctions.h>
#include <../Project1/project1/lib.h>

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
    WaveFunctions atom;

    VMCSolver() {}
    VMCSolver(WaveFunctions a): idum(-1), step(2.55), AcceptRate(0), h(0.001)
    {atom = a; charge = atom.Nparticles; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions;}

    void Integration();
    void MCintegration(int N);
    double LocalEnergy(mat R);
};


#endif // VMCSOLVER_H
