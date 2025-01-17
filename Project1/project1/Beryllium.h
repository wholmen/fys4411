#ifndef BERYLLIUM_H
#define BERYLLIUM_H

#include<lib.h>
#include<iostream>
#include<armadillo>

using namespace std;
using namespace arma;


class Beryllium
{
private:
    double phi1s(double r);
    double phi2s(double r);
    double corrolation(mat R);
    double sign(double a);

public:
    double alpha;
    double beta;
    double charge;
    double Z;
    double Nparticles;
    double Ndimensions;
    double h;
    Beryllium() {}
    Beryllium(double alpha, double beta): alpha(alpha), beta(beta), charge(4), Z(4), Nparticles(4), Ndimensions(3), h(0.01) {}
    double WaveFunction(mat R);
    double LocalEnergy(mat R);
    double LocalEnergyClosed(mat R);
};


#endif // BERYLLIUM_H

