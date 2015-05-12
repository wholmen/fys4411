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
    double GTO(double r, double xi);
    double dGTO(double r, double xi);
    double d2GTO(double r, double xi);
public:
    int n,l,ml,ms,Nparticles;
    double alpha;
    GaussianOrbitals() {}
    GaussianOrbitals(int n1, int l1, int ml1, int ms1, double a, int Nparticles1): n(n1), l(l1), ml(ml1), ms(ms1), Nparticles(Nparticles1), alpha(a) {

    }
    double psi(double r);
    double dpsi(double r);
    double d2psi(double r);
};

#endif // GAUSSIANORBITALS_H
