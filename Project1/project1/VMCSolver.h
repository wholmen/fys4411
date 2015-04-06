#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <iostream>
#include <armadillo>
#include <lib.h>
#include <Atom.h>

using namespace std;
using namespace arma;

class VMCSolver
{

private:
    int Nstep;
    long idum;

    double Nparticles;
    double Ndimensions;

public:
    double step;
    double Energy;
    double Variance;
    double AcceptRate;

    Atom atom;

    VMCSolver() {}
    VMCSolver(Atom a): Nstep(100000), idum(-1), step(2.55), AcceptRate(0)
    {atom = a; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions;}

    void MCintegration();
    void FindStepLength();
    void ImportanceSampling();
};


#endif // VMCSOLVER_H
