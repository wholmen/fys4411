#ifndef HELIUM_H
#define HELIUM_H

#include<iostream>
#include<armadillo>
#include<lib.h>

using namespace arma;
using namespace std;

class Helium
{
private:
public:
    double alpha;
    double beta;
    double charge;
    double Nparticles;
    double Ndimensions;
    int WFnumber;
    double h;
    double sqh;

    Helium() {}
    Helium(double Alpha, double Beta, int WFNumber): alpha(Alpha), beta(Beta), charge(2), Nparticles(2), Ndimensions(3), WFnumber(WFNumber), h(0.001), sqh(sqrt(h)) {}

    double WaveFunction(mat R);
    double LocalEnergy(mat R);
    void QuantumForce(mat r, mat &Qforce);
};

#endif // HELIUM_H
