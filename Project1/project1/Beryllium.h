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

public:
    double alpha;
    double beta;
    double charge;
    double Z;
    double Nparticles;
    double Ndimensions;
    Beryllium() {}
    Beryllium(double alpha, double beta): alpha(alpha), beta(beta), charge(4), Z(4), Nparticles(4), Ndimensions(3) {}
    double WaveFunction(mat R);
    double LocalEnergy(mat R);
};


#endif // BERYLLIUM_H

