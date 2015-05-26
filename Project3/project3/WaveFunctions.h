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
    int Nparticles, Ndimensions, Np2;
    double alpha, beta, h;
    mat DinvUp, DinvDown, DinvUpQF, DinvDownQF, Dup, Ddown;
    HydrogenOrbitals *particles;

    // Initializing the class
    WaveFunctions(){}
    WaveFunctions(double Alpha, double Beta, int nparticles): Nparticles(nparticles), Ndimensions(3), alpha(Alpha), beta(Beta), h(0.001)
    {
    Np2 = Nparticles/2;
    particles = new HydrogenOrbitals[Nparticles];
    Dup = zeros<mat>(Np2,Np2); Ddown = zeros<mat>(Np2,Np2);
    }

    //(GaussianBasis == false) ? HydrogenOrbitals *particles : GaussianOrbitals *particles;

    // Initializing the system
    void Initialize_System(mat r);
    void UpdateInverseSD(int i, bool QF);
    void UpdateInverseQF(bool Move);
    void FactorizeDeterminant(mat &SDup, mat &SDdown, mat R, int d);

    // This Function returns the wave function
    double WaveFunction(mat r, bool DifferentiateBeta);

    // Calculating quantum force
    void QuantumForce(mat r, mat &Qforce);

    // To calculate kinetic energy
    double KineticEnergy(mat R);
};


#endif // WAVEFUNCTIONS_H
