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
    double factorial(int a);

public:
    int Nparticles, Ndimensions;
    double alpha, beta;
    bool GaussianBasis;
    mat D, DInverse;

    // Initializing the class
    WaveFunctions(){}
    WaveFunctions(double Alpha, double Beta, int nparticles): Nparticles(nparticles), Ndimensions(3), alpha(Alpha), beta(Beta), GaussianBasis(false) {}

    // This Function returns the wave function
    double WaveFunction(mat r);

    // To calculate kinetic energy
    double diff2_Slater(mat R);
    double diff_Slater_diff_Corrolation(mat R);
    double diff2_Corrolation(mat R);
};


#endif // WAVEFUNCTIONS_H
