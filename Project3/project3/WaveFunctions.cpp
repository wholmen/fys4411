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

void WaveFunctions::Initialize_System(mat r){
    int n(0), l(0), ml(0), ms(1);
    // Making a list of all the different hydrogen orbitals with the different wave functions.
    for (int i=0; i<Nparticles; i++){
        if (ms < 1){ms ++; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha,Nparticles);}
        else{
            if (ml < l) {ml ++; ms = 0; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha,Nparticles);}
            else {
                if (l < n - 1){l++; ml = -l; ms = 0; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha,Nparticles);}
                else {
                   n++; l = 0; ml = -l; ms = 0; particles[i] = GaussianOrbitals(n,l,ml,ms,alpha,Nparticles);
                }
            }
        }
    }
    // Computing the Slater matrix using the list above. Splitting determinant into Dup and Ddown
    int y, e, i, j;
    //mat Dup(Np2,Np2), Ddown(Np2,Np2);

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
    DinvUp = inv(Dup); DinvDown = inv(Ddown);
}

double WaveFunctions::WaveFunction(mat r){
    // Computing the Slater matrix using the list above. Splitting determinant into Dup and Ddown
    int y, e, i, j;
    double SD, corrolation(1), a(0), d;
    //mat Dup(Np2,Np2), Ddown(Np2,Np2);

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
    return SD*corrolation; // Returning the Slater determinant mulitplied with the Jastrow factor.
}

void WaveFunctions::UpdateInverseSD(mat r, int i){
    // Step accepted, now update the inverse matrices.
    mat DinvUpNew = zeros<mat>(Np2,Np2); mat DinvDownNew = zeros<mat>(Np2,Np2); double Sj(0), detratio(0); int j,l;

    if (i < Np2) for (l=0; l<Np2; l++) {detratio += Dup(l,i) * DinvUp(l,i);}
    else         for (l=0; l<Np2; l++) {detratio += Ddown(l,i-Np2) * DinvDown(l,i-Np2);}

    if (i < Np2) {
        for (j=0; j<Np2; j++){
            if (j!=i){
                for (l=0; l<Np2; l++) Sj += Dup(l,i) * DinvUp(l,j);
                for (l=0; l<Np2; l++) DinvUpNew(l,j) = DinvUp(l,j) - Sj*DinvUp(l,i) / detratio;
            }
        }
        for (l=0; l<Np2; l++) DinvUpNew(l,i) = DinvUp(l,i) / detratio;
        DinvUp = DinvUpNew;
    }
    else{
        for (j=0; j<Np2; j++){
            if (j!=i-Np2){
                for (l=0; l<Np2; l++) Sj += Ddown(l,i-Np2) * DinvDown(l,j);
                for (l=0; l<Np2; l++) DinvDownNew(l,j) = DinvDown(l,j) - Sj*DinvDown(l,i-Np2) / detratio;
            }
        }
        for (l=0; l<Np2; l++) DinvDownNew(l,i-Np2) = DinvDown(l,i-Np2) / detratio;
        DinvDown = DinvDownNew;
    }
}

double WaveFunctions::diff2_Slater(mat R){
    double T(0); int y, e, i, j;

    mat D2up(Np2,Np2), D2down(Np2,Np2);
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 0) continue;
            D2up(e,j) = particles[i].d2psi(length(R.row(j)));
            if (j == Np2 - 1) e++; // increasing index in Dup
        }
    }
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 1) continue;
            y = j + Np2;
            D2down(e,j) = particles[i].d2psi(length(R.row(y)));
            if (j == Np2 -1) e++;
        }
    }

    for (i=0; i<Nparticles; i++){
        if (i < Np2) for (j=0; j<Np2; j++) T += D2up(j,i) * DinvUp(j,i);
        else         for (j=0; j<Np2; j++) T += D2down(j,i-Np2) * DinvDown(j,i-Np2);
    }
    return T;
}

double WaveFunctions::diff_Slater_diff_Corrolation(mat R){
    // Computing the Slater part
    double T(0); int y, e, i, j; rowvec detpart = zeros<rowvec>(Ndimensions);

    mat D1up(Np2,Np2), D1down(Np2,Np2);
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 0) continue;
            D1up(e,j) = particles[i].dpsi(length(R.row(j)));
            if (j == Np2 - 1) e++; // increasing index in Dup
        }
    }
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 1) continue;
            y = j + Np2;
            D1down(e,j) = particles[i].dpsi(length(R.row(y)));
            if (j == Np2 -1) e++;
        }
    }


    rowvec rij; double aij, Rij; rowvec Jastrow(Ndimensions);
    for (i=0; i<Nparticles; i++){
        detpart = zeros<rowvec>(Ndimensions); Jastrow = zeros<rowvec>(Ndimensions);

        if (i < Np2) for (j=0; j<Np2; j++) detpart += D1up(j,i) * DinvUp(j,i) * R.row(i) / length(R.row(i));
        else         for (j=0; j<Np2; j++) detpart += D1down(j,i-Np2) * DinvDown(j,i-Np2) * R.row(i) / length(R.row(i));

        //cout << " detpart " << detpart << " analytical " << -alpha*R.row(i)/length(R.row(i)) << endl;

        for (int j=0; j<Nparticles; j++){
            if (i==j) continue;
            rij = R.row(j)-R.row(i); Rij = length(rij);
            aij = (particles[i].ms == particles[j].ms) ? 0.25 : 0.5;
            Jastrow += rij / Rij * aij / pow(1+beta*Rij,2);
        }
        //cout << " Jastrow " << Jastrow << " analytical " << pow(-1,i+1) * (R.row(0) - R.row(1)) / (2*length(R.row(0) - R.row(1))*pow(1+beta*length(R.row(0) - R.row(1)),2)) << endl;

        T -= dot(detpart,Jastrow);
    }
    //double R12 = length(R.row(0) - R.row(1) ); double r1 = length(R.row(0)); double r2 = length(R.row(1));

    //cout << "Numerical T is " << T;
    //T = -1.0 / (2*pow(1+beta*R12,2)) * alpha*(r1+r2)/R12 * (1 - dot(R.row(0),R.row(1))/r1*r2);
    //cout << " And Analytical T is " << T << endl;

    return T;
}

double WaveFunctions::diff2_Corrolation(mat R){
    double T(0), Rki,Rkj,aki,akj; rowvec rki, rkj;
    for (int k=0; k<Nparticles; k++){
        for (int i=0; i<Nparticles; i++){
            for (int j=0; j<Nparticles; j++){
                if (i==k || j==k) continue;
                rki = R.row(k) - R.row(i); rkj = R.row(k) - R.row(j);
                Rki = length(rki); Rkj = length(rkj);
                aki = (particles[j].ms != particles[i].ms) ? 0.25: 0.5;
                akj = (particles[i].ms != particles[j].ms) ? 0.25: 0.5;
                //T += dot(rki, rkj) / (Rki*Rkj) * aki / ((1 + beta*Rki)*(1+beta*Rki)) * akj / ((1 + beta*Rkj)*(1+beta*Rkj));
                T += dot(rki, rkj) / (Rki*Rkj) * aki / pow(1 + beta*Rki,2) * akj / pow(1 + beta*Rkj,2);
            }
        if (i == k) continue;
        //T += 2*aki / (Rki * (1+beta*Rki)*(1+beta*Rki)) - 2*aki*beta / ((1+beta*Rki)*(1+beta*Rki)*(1+beta*Rki));
        T += 2*aki / (Rki * pow(1+beta*Rki,2)) - 2*aki*beta / pow(1+beta*Rki,3);
        }
    }
    //cout << T;
    //T = 0;
    //double R12 = length(R.row(0) - R.row(1) );
    //T = 2* ( 1.0 / (R12*pow(1+beta*R12,2)) + 0.25 / pow(1 + beta*R12,4) - beta / pow(1+beta*R12,3) );
    //cout << " " << T << endl;
    return T;
}




