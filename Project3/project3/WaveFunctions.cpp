#include <WaveFunctions.h>

double WaveFunctions::length(rowvec r){
    double l=0;
    for (int j=0; j<Ndimensions; j++) l += r(j)*r(j);
    return sqrt(l);
}

double WaveFunctions::distance(rowvec r1, rowvec r2){
    double dist = 0;
    for (int i=0; i<Ndimensions; i++) dist += (r1(i) - r2(i))*(r1(i) - r2(i));
    return sqrt(dist);
}

double WaveFunctions::factorial(int a){
    double factor = 1;
    for (int i=a; i>0; i--) factor *= i;
    return factor;
}

double WaveFunctions::WaveFunction(mat r){
    // Making a list of all the different hydrogen orbitals with the different wave functions.
    HydrogenOrbitals *particles = new HydrogenOrbitals[Nparticles];
    int n(0), l(0), ml(0), ms(1);

    for (int i=0; i<Nparticles; i++){
        if (ms < 1){ms ++; particles[i] = HydrogenOrbitals(n,l,ml,ms,alpha);}
        else{
            if (ml < l) {ml ++; ms = 0; particles[i] = HydrogenOrbitals(n,l,ml,ms,alpha);}
            else {
                if (l < n - 1){l++; ml = -l; ms = 0; particles[i] = HydrogenOrbitals(n,l,ml,ms,alpha);}
                else {
                   n++; l = 0; ml = -l; ms = 0; particles[i] = HydrogenOrbitals(n,l,ml,ms,alpha);
                }
            }
        }
    }
    // Computing the Slater matrix using the list above. Splitting determinant into Dup and Ddown
    int Np2(Nparticles/2), y, e, i, j;
    double SD, corrolation(1), a(0), d;
    mat Dup(Np2,Np2), Ddown(Np2,Np2);

    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 0) continue;
            Dup(e,j) = particles[i].psi(length(r.row(j)));
            if (j == Np2 - 1) e++; // increasing index in Dup
        }
    }
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 1) continue;
            y = j + Np2;
            Ddown(e,j) = particles[i].psi(length(r.row(y)));
            if (j == Np2 -1) e++;
        }
    }
    SD = det(Dup) * det(Ddown);

    // Computing the Corrolation factor
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<i; j++){
            (particles[i].ms == particles[j].ms) ? a = 0.25: a = 0.5;
            d = distance(r.row(i), r.row(j));
            corrolation *= exp( (a*d) / (1.0 + beta*d));
        }
    }
    return SD;//*corrolation; // Returning the Slater determinant mulitplied with the Jastrow factor.
}

double WaveFunctions::diff2_Slater(mat R){
    double T(0);
    if (GaussianBasis == false){
        for (int i=0; i<Nparticles; i++){
            T += alpha * alpha - 2 * alpha / length( R.row(i) );
        }
    }
    return T;
}

double WaveFunctions::diff_Slater_diff_Corrolation(mat R){

}

double WaveFunctions::diff2_Corrolation(mat R){

}














/*
double WaveFunctions::WaveFunction(mat r){

    // Making a list of all the different hydrogen orbitals with the different wave functions.
    GaussianOrbitals *particles = new GaussianOrbitals[Nparticles];
    int n = 0; int l = 0; int ml = 0; int ms = 1;

    for (int i=0; i<Nparticles; i++){
        if (ms < 1){ms ++; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha);}
        else{
            if (ml < l) {ml ++; ms = 0; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha);}
            else {
                if (l < n - 1){l++; ml = -l; ms = 0; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha);}
                else {
                   n++; l = 0; ml = -l; ms = 0; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha);
                }
            }
        }
    }

    // Computing the Slater matrix using the list above.
    mat SlaterMatrix = zeros<mat>(Nparticles,Nparticles);

    for (int i=0; i<Nparticles; i+=1){
        for (int j=0; j<Nparticles; j+=1){
            SlaterMatrix(i,j) = particles[j].psi(length(r.row(i)));
        }
    }

    //cout << SlaterMatrix << det(SlaterMatrix);
    // Taking the determinant of this matrix gives the Slater Determinant.
    double SD = det(SlaterMatrix);


    // Computing the Corrolation factor
    double corrolation = 0; double a = 0;
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<i; j++){
            if (particles[i].ms == particles[j].ms) {a = 0.25;}
            else {a = 0.5; }

            double d = distance(r.row(i), r.row(j));
            corrolation *= exp( (a*d) / (1 + beta*d));
        }
    }

    // Returning the Slater determinant mulitplied with the Jastrow factor.
    return SD;//*corrolation;

}
*/
