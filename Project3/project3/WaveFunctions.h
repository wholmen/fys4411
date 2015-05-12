#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H

#include <armadillo>
#include <iostream>
#include <HydrogenOrbitals.h>
#include <GaussianOrbitals.h>
//#include <Orbitals.h>

using namespace std;
using namespace arma;


class WaveFunctions
{
private:
    double distance(rowvec r1, rowvec r2);
    double length(rowvec r);

public:
    int Nparticles, Ndimensions, Np2;
    double alpha, beta;
    bool GaussianBasis;
    mat DinvUp, DinvDown, Dup, Ddown;
    //HydrogenOrbitals *particles;
    GaussianOrbitals *particles;

    // Initializing the class
    WaveFunctions(){}
    WaveFunctions(double Alpha, double Beta, int nparticles): Nparticles(nparticles), Ndimensions(3), alpha(Alpha), beta(Beta), GaussianBasis(false)
    {
    Np2 = Nparticles/2;
    particles = new GaussianOrbitals[Nparticles]; GaussianBasis = true;
    //particles = new HydrogenOrbitals[Nparticles];
    Dup = zeros<mat>(Np2,Np2); Ddown = zeros<mat>(Np2,Np2);
    }

    //(GaussianBasis == false) ? HydrogenOrbitals *particles : GaussianOrbitals *particles;

    // Initializing the system
    void Initialize_System(mat r);
    void UpdateInverseSD(mat R, int i);

    // This Function returns the wave function
    double WaveFunction(mat r);

    // To calculate kinetic energy
    double diff2_Slater(mat R);
    double diff_Slater_diff_Corrolation(mat R);
    double diff2_Corrolation(mat R);
};


#endif // WAVEFUNCTIONS_H
