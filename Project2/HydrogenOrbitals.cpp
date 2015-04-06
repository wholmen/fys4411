#include <HydrogenOrbitals.h>

double HydrogenOrbitals::psi1s(double r){return exp(-alpha*r);}

double HydrogenOrbitals::psi2s(double r){return (1-alpha*r/2.0)*exp(-alpha*r/2.0);}

double HydrogenOrbitals::psi2p(double r){return alpha*r*exp(-alpha*r/2.0);}

double HydrogenOrbitals::psi(double r){
    if (n==1) {return psi1s(r);}
    else if (n==2) {
        if (l==0) {return psi2s(r);}
        if (l==1) {return psi2p(r);}
        else {return false;}
    }
    else {return false;}
}
