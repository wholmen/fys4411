#ifndef GAUSSIANORBITALS_H
#define GAUSSIANORBITALS_H

#include <armadillo>
#include <iostream>
#define pi 3.1415

using namespace std;
using namespace arma;

class GaussianOrbitals
{
private:
public:
    int n,l,ml,ms, Nparticles;

    GaussianOrbitals() {}
    GaussianOrbitals(int n1, int l1, int ml1, int ms1, int np): n(n1), l(l1), ml(ml1), ms(ms1), Nparticles(np) { }

    // Wavefunctions
    double psi(rowvec r, int row);
    rowvec dpsi(rowvec r, int row);
    double d2psi(rowvec r, int row);

    // Help functions
    double N(double a, int i, int j, int k);
    double factorial(int a);
    double Length(rowvec r);
    double Xi(rowvec r, int i, int j, int k, double a, double c);
    rowvec XiD(rowvec r, int i, int j, int k, double a, double c);
    double Xi2D(rowvec r, int i, int j, int k, double a, double c);
};


#endif // GAUSSIANORBITALS_H
