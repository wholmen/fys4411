#ifndef GAUSSIANORBITALS_H
#define GAUSSIANORBITALS_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class GaussianOrbitals
{
private:
    double psi1s(double r);
    double psi2s(double r);
    double psi2p(double r);
public:
    int n,l,ml,ms;
    double alpha;
    GaussianOrbitals() {}
    GaussianOrbitals(int n1, int l1, int ml1, int ms1, double a): n(n1), l(l1), ml(ml1), ms(ms1), alpha(a) { }
    double psi(double r);
};

#endif // GAUSSIANORBITALS_H
