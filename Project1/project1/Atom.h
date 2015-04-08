#ifndef ATOM_H
#define ATOM_H
#include <Beryllium.h>
#include <Helium.h>

class Atom{
public:
    double Nparticles;
    double Ndimensions;
    double charge;
    double alpha;
    double beta;
    Helium atom;
    Beryllium atom2;

    Atom(Helium H) {atom = H; Nparticles=H.Nparticles; Ndimensions=H.Ndimensions; charge=H.charge; alpha=H.alpha; beta=H.beta;}
    Atom(Beryllium B) {atom2 = B; Nparticles=B.Nparticles; Ndimensions=B.Ndimensions; charge=B.charge; alpha=B.alpha; beta=B.beta;}
    Atom() {}

    double WaveFunction(mat R) {return (charge==2) ? atom.WaveFunction(R) : atom2.WaveFunction(R);}
    double LocalEnergy(mat R) {return (charge==2) ? atom.LocalEnergy(R) : atom2.LocalEnergy(R);}
    double LocalEnergyClosed(mat R) {return (charge == 2) ? atom.LocalEnergyClosed(R) : atom2.LocalEnergyClosed(R);}
    void QuantumForce(mat r, mat &Qforce) {atom.QuantumForce(r, Qforce);}
};


#endif // ATOM_H
