#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H

#include <armadillo>
#include <iostream>
#include <HydrogenOrbitals.h>

using namespace std;
using namespace arma;

class WaveFunctions
{
private:
    double distance(rowvec r1, rowvec r2);
    double length(rowvec r);

public:
    int Nparticles, Ndimensions;
    double alpha, beta;
    WaveFunctions(){}
    WaveFunctions(double Alpha, double Beta, int nparticles): Nparticles(nparticles), Ndimensions(3), alpha(Alpha), beta(Beta) {}
    double WaveFunction(mat r);
};


#endif // WAVEFUNCTIONS_H
