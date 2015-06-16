#ifndef HYDROGENORBITALS_H
#define HYDROGENORBITALS_H

#include <armadillo>
#include <iostream>
#include <complex>

using namespace std;
using namespace arma;

class HydrogenOrbitals
{
private:
    double psi1s(double r, double d);
    double psi2s(double r, double d);
    double psi2p(rowvec R, double r, double d);
public:
    int n,l,ml,ms;
    double alpha;
    HydrogenOrbitals() {}
    HydrogenOrbitals(int n1, int l1, int ml1, int ms1, double a): n(n1), l(l1), ml(ml1), ms(ms1), alpha(a) { }
    double psi(double r, rowvec R);
    double dpsi(double r, rowvec R);
    double d2psi(double r, rowvec R);
};


#endif // HYDROGENORBITALS_H
