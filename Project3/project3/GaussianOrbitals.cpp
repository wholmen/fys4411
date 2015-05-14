#include <GaussianOrbitals.h>

double GaussianOrbitals::factorial(int a){
    double fact(1.0);
    for (int i=a; i>0; i--) fact *= i;
    return fact;
}

double GaussianOrbitals::Normalization(double a, int i, int j, int k){
    return pow(2*a/pi,0.75) * sqrt( pow(8*a,i+j+k) * factorial(i)*factorial(j)*factorial(k)/factorial(2*i)/factorial(2*j)/factorial(2*k) );
}

double GaussianOrbitals::KineticEnergy(mat R){
    return 1;
}

double GaussianOrbitals::Xi(double w, double a, double x, double y, double z, int i, int j, int k){
    double r2 = x*x + y*y + z*z;
    return Normalization(a,i,j,k) * w * pow(x,i) * pow(y,j) * pow(z,k) * exp(-a*r2);
}

double GaussianOrbitals::WaveFunction(mat R){
    double wavefunc(0.0);
    if (Nparticles == 2){
        // Only 1s-orbitals
        for (int i=0; i<Nparticles-1; i++){
            wavefunc += Xi(0.175230, 13.62670, R(i,0), R(i,1), R(i,2), 0, 0, 0);
            wavefunc += Xi(0.893483, 1.999350, R(i,0), R(i,1), R(i,2), 0, 0, 0);
            wavefunc += Xi(1.000000, 0.382993, R(i,0), R(i,1), R(i,2), 0, 0, 0);
            //cout << wavefunc << endl;
        }
    }
    else if (Nparticles == 4){

    }
    else if (Nparticles == 10){

    }
    //cout << ;
    return wavefunc;
}

double GaussianOrbitals::dpsi(double r){
    double sum(0);
    if (Nparticles == 2){

    }
    return sum;
}

double GaussianOrbitals::d2psi(double r){
    double sum(0);
    if (Nparticles == 2){

    }
    return sum;
}


