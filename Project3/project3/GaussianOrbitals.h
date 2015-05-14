#ifndef GAUSSIANORBITALS_H
#define GAUSSIANORBITALS_H

#include <armadillo>
#include <iostream>

#define pi 3.14159265359

using namespace std;
using namespace arma;

class GaussianOrbitals
{
private:
    double Normalization(double exp, int i, int j, int k);
    double factorial(int a);
    double Xi(double w, double exp, double x, double y, double z, int i, int j, int k);

public:
    int Nparticles, Ndimensions;

    GaussianOrbitals() {}
    GaussianOrbitals(int nparticles): Nparticles(nparticles), Ndimensions(3) {

    }
    double WaveFunction(mat R);
    double dpsi(double r);
    double d2psi(double r);

    double KineticEnergy(mat R);
};

#endif // GAUSSIANORBITALS_H
