#include <HydrogenOrbitals.h>
#define pi 3.1415
typedef complex<double> dcomp;

double HydrogenOrbitals::psi1s(double r, double d){
    if (d == 0) return exp(-alpha*r);
    if (d == 1) return exp(-alpha*r)*(-alpha);
    if (d == 2) return exp(-alpha*r)*alpha* (alpha - 2/r);
    else return false;
    }

double HydrogenOrbitals::psi2s(double r, double d){
    if (d == 0) return (1-alpha*r/2.0)*exp(-alpha*r/2.0);
    if (d == 1) return alpha*(0.25*alpha*r - 1.0)*exp(-0.5*alpha*r);
    if (d == 2) return alpha*(alpha*r*(-0.125*alpha*r + 0.75) + 0.5*alpha*r - 2.0)*exp(-0.5*alpha*r)/r;
    else return false;
}

double HydrogenOrbitals::psi2p(double r, double d){
    //double Y(1), theta;
    //dcomp i; i = sqrt(-1); dcomp phi;
    //theta = acos(z/r);
    //phi = atan(y/r);

    //if (ml == 0) Y = sqrt(3/4/pi);
    //if (ml == 1) Y = sqrt(3/8/pi)*sin(theta)*exp(-i*phi);
    //if (ml == -1)Y = sqrt(3/8/pi)*sin(theta)*exp( i*phi);

    if (d == 0 && ml == -1) return (alpha*r*exp(-alpha*r/2.0));
    if (d == 0 && ml == 0) return (alpha*r*exp(-alpha*r/4.0));
    if (d == 0 && ml == 1) return (alpha*r*exp(-alpha*r/8.0));
    if (d == 1) return (alpha*(-0.5*alpha*r + 1)*exp(-0.5*alpha*r));
    if (d == 2) return (alpha*(alpha*r*(0.25*alpha*r - 1.0) - 1.0*alpha*r + 2)*exp(-0.5*alpha*r)/r);
    else return false;
}

double HydrogenOrbitals::psi(double r){
    if (n==1) return psi1s(r,0);
    else if (n==2 && l==0) return psi2s(r,0);
    else if (n==2 && l==1) return psi2p(r,0);
    else return false;
}

double HydrogenOrbitals::dpsi(double r){
    if (n==1) return psi1s(r,1);
    else if (n==2 && l==0) return psi2s(r,1);
    else if (n==2 && l==1) return psi2p(r,1);
    else return false;
}

double HydrogenOrbitals::d2psi(double r){
    if (n==1) return psi1s(r,2);
    else if (n==2 && l==0) return psi2s(r,2);
    else if (n==2 && l==1) return psi2p(r,2);
    else return false;
}
