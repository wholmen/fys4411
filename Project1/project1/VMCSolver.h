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
    bool closed;

public:
    double step;
    double Energy;
    double Variance;
    double AcceptRate;
    double Phitotal;

    Atom atom;

    VMCSolver() {}
    VMCSolver(Atom a, bool CLOSED): Nstep(10000), idum(-1), step(1.5), AcceptRate(0)
    {atom = a; Nparticles = atom.Nparticles; Ndimensions = atom.Ndimensions; closed = CLOSED;}

    void MCintegration();
    void FindStepLength();
    void ImportanceSampling(double DT);
};


#endif // VMCSOLVER_H
