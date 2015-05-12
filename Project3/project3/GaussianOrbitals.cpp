#include <GaussianOrbitals.h>


double GaussianOrbitals::GTO(double r, double xi){
    return 2*xi/pow(pi,0.75) *exp(-xi*r*r);
}

double GaussianOrbitals::dGTO(double r, double xi){
    return 2*xi/pow(pi,0.75) * 2 * xi * r * exp(-xi*r*r);
}

double GaussianOrbitals::d2GTO(double r, double xi){
    return 4*xi*xi/pow(pi,0.75) * exp(-xi*r*r) * (3 + 2*xi*r*r);
}

double GaussianOrbitals::psi(double r){
    double sum(0);
    if (Nparticles == 2){
        sum += 0.17523*GTO(r,13.6267);
        sum += 0.893483*GTO(r,1.99935);
        sum += 1*GTO(r,0.382993);
    }
    return sum;
}

double GaussianOrbitals::dpsi(double r){
    double sum(0);
    if (Nparticles == 2){
        sum += 0.17523*dGTO(r,13.6267);
        sum += 0.893483*dGTO(r,1.99935);
        sum += 1*dGTO(r,0.382993);
    }
    return sum;
}

double GaussianOrbitals::d2psi(double r){
    double sum(0);
    if (Nparticles == 2){
        sum += 0.17523*d2GTO(r,13.6267);
        sum += 0.893483*d2GTO(r,1.99935);
        sum += 1*d2GTO(r,0.382993);
    }
    return sum;
}


