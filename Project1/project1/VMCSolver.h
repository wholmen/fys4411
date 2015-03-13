#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <iostream>
#include <armadillo>
#include <lib.h>

using namespace std;
using namespace arma;

class VMCSolver
{

private:
    int Nstep;
    int Nparticles;
    int Ndimensions;
    double charge;
    double h;
    double sqh;
    long idum;
    double LocalEnergy(mat R);
    double WaveFunction(mat R);

public:
    double alpha;
    double step;
    double beta;
    double Energy;
    double Variance;
    double AcceptRate;
    int WFnumber;
    VMCSolver();
    void MCintegration();
    void FindStepLength();
    void ImportanceSampling();
    void QuantumForce(mat r, mat &Qforce);
};


#endif // VMCSOLVER_H
