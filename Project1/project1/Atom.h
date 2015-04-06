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

    double WaveFunction(mat R) {if (charge==2) return atom.WaveFunction(R); else return atom2.WaveFunction(R);}
    double LocalEnergy(mat R) {if (charge==2) return atom.LocalEnergy(R); else return atom2.LocalEnergy(R);}
    void QuantumForce(mat r, mat &Qforce) {atom.QuantumForce(r, Qforce);}
};


#endif // ATOM_H
