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
    // for (int i=0; i<10; i++){ cout << particles[i].n << particles[i].l << particles[i].ml << particles[i].ms << endl; }
    // Computing the Slater matrix using the list above. Splitting determinant into Dup and Ddown
    FactorizeDeterminant(Dup,Ddown,r,0);

    for (int i=0; i<10; i++) { cout << particles[i].l; }

    cout << 1;
    DinvUp = inv(Dup); DinvDown = inv(Ddown);
    DinvUpQF = DinvUp; DinvDownQF = DinvDown;
    cout << 2;
}

double WaveFunctions::WaveFunction(mat r, bool DifferentiateBeta){
    // Computing the Slater matrix using the list above. Splitting determinant into Dup and Ddown
    double SD, corrolation(1), a(0), d;

    FactorizeDeterminant(Dup,Ddown,r,0);
    SD = det(Dup) * det(Ddown);

    // Computing the Corrolation factor
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<i; j++){
            (particles[i].ms == particles[j].ms) ? a = 0.25: a = 0.5;
            d = distance(r.row(i), r.row(j));
            if (DifferentiateBeta) corrolation *= -a*d*d / pow((1 + beta*d),2) * exp( (a*d) / (1.0 + beta*d));
            else corrolation *= exp( (a*d) / (1.0 + beta*d));
        }
    }
    return SD*corrolation; // Returning the Slater determinant mulitplied with the Jastrow factor.
}

void WaveFunctions::UpdateInverseSD(int i, bool QF){
    // bool QF determines wether we update inverse matrix for computing the quantum force. If QF==true, the updated D-1 will not
    // be used to compute ELocal, but will be used to update QForceNew. Function UpdateInverseQF will update D-1 to compute EL. Called if move is accepted.
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
        (not QF) ? DinvUp = DinvUpNew : DinvUpQF = DinvUpNew;
    }
    else{
        for (j=0; j<Np2; j++){
            if (j!=i-Np2){
                for (l=0; l<Np2; l++) Sj += Ddown(l,i-Np2) * DinvDown(l,j);
                for (l=0; l<Np2; l++) DinvDownNew(l,j) = DinvDown(l,j) - Sj*DinvDown(l,i-Np2) / detratio;
            }
        }
        for (l=0; l<Np2; l++) DinvDownNew(l,i-Np2) = DinvDown(l,i-Np2) / detratio;
        (not QF) ? DinvDown = DinvDownNew : DinvDownQF = DinvDownNew;
    }
}

void WaveFunctions::UpdateInverseQF(bool Move){
    // If move is accepted. DinvUp and DinvDown must be updated. Since we have already computed them, but stored them in the intermidiate values
    // DinvUpQF and DinvDownQF, we can obtain the values from them.
    (Move) ? DinvUp = DinvUpQF : DinvUpQF = DinvUp;
    (Move) ? DinvDown = DinvDownQF : DinvDownQF = DinvDown;
}

void WaveFunctions::QuantumForce(mat r, mat &Qforce){
    bool Analytical = false;

    if (Analytical){
        mat SDup(Np2,Np2), SDdown(Np2,Np2);
        FactorizeDeterminant(SDup,SDdown,r,1);

        rowvec detpart, Jastrow, rij;
        double a, Rij; int i,j;

        for (i=0; i<Nparticles; i++){
            detpart = zeros<rowvec>(Ndimensions); Jastrow = zeros<rowvec>(Ndimensions);

            if (i < Np2) for (j=0; j<Np2; j++) detpart += SDup(j,i) * DinvUpQF(j,i) * r.row(i) / length(r.row(i));
            else         for (j=0; j<Np2; j++) detpart += SDdown(j,i-Np2) * DinvDownQF(j,i-Np2) * r.row(i) / length(r.row(i));

            for (int j=0; j<Nparticles; j++){
                if (i==j) continue;
                rij = r.row(i)-r.row(j); Rij = length(rij);
                a = (particles[i].ms == particles[j].ms) ? 0.25 : 0.5;
                Jastrow += rij / Rij * a / pow(1+beta*Rij,2);
            }
            Qforce.row(i) = 2*(detpart + Jastrow);
        }
    }
    else{
        mat Rplus = zeros<mat>(Nparticles,Ndimensions);
        mat Rminus = zeros<mat>(Nparticles,Ndimensions);
        double Psiplus = 0; double Psiminus = 0;

        Rplus = Rminus = r;
        double Psi = WaveFunction(r,false);

        for (int i=0; i<Nparticles; i++){
            for (int j=0; j<Ndimensions; j++){
                Rplus(i,j)  += h;
                Rminus(i,j) -= h;
                Psiplus = WaveFunction(Rplus,false);
                Psiminus = WaveFunction(Rminus,false);

                Qforce(i,j) = Psiplus - Psiminus;
                Rplus(i,j) = r(i,j); Rminus(i,j) = r(i,j);
            }
        }
        Qforce = Qforce * h / Psi;
    }
}

void WaveFunctions::FactorizeDeterminant(mat &SDup, mat &SDdown, mat R, int d){
    int e, i, y, j;
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 0) continue;
            if (d==0) SDup(e,j) = particles[i].psi(length(R.row(j)));
            if (d==1) SDup(e,j) = particles[i].dpsi(length(R.row(j)));
            if (d==2) SDup(e,j) = particles[i].d2psi(length(R.row(j)));
            if (j == Np2 - 1) e++; // increasing index in Dup
        }
    }
    for (e=0, i=0; i<Nparticles; i++){
        for (j=0; j<Np2; j++){
            if (particles[i].ms != 1) continue;
            y = j + Np2;
            if (d==0) SDdown(e,j) = particles[i].psi(length(R.row(y)));
            if (d==1) SDdown(e,j) = particles[i].dpsi(length(R.row(y)));
            if (d==2) SDdown(e,j) = particles[i].d2psi(length(R.row(y)));
            if (j == Np2 -1) e++;
        }
    }
}

double WaveFunctions::KineticEnergy(mat R){
    double T(0); int i, j;
    double a, Rij, Rki, Rkj; rowvec detpart, Jastrow, rij, rki, rkj;
    mat D2up(Np2,Np2), D2down(Np2,Np2), D1up(Np2,Np2), D1down(Np2,Np2);

    FactorizeDeterminant(D2up,D2down,R,2);
    FactorizeDeterminant(D1up, D1down, R, 1);

    for (i=0; i<Nparticles; i++){
        if (i < Np2) for (j=0; j<Np2; j++) T += D2up(j,i) * DinvUp(j,i);
        else         for (j=0; j<Np2; j++) T += D2down(j,i-Np2) * DinvDown(j,i-Np2);
    }

    for (i=0; i<Nparticles; i++){
        detpart = zeros<rowvec>(Ndimensions); Jastrow = zeros<rowvec>(Ndimensions);

        if (i < Np2) for (j=0; j<Np2; j++) detpart += D1up(j,i) * DinvUp(j,i) * R.row(i) / length(R.row(i));
        else         for (j=0; j<Np2; j++) detpart += D1down(j,i-Np2) * DinvDown(j,i-Np2) * R.row(i) / length(R.row(i));

        for (int j=0; j<Nparticles; j++){
            if (i==j) continue;
            rij = R.row(j)-R.row(i); Rij = length(rij);
            a = (particles[i].ms == particles[j].ms) ? 0.25 : 0.5;
            Jastrow += rij / Rij * a / pow(1+beta*Rij,2);
        }
        T -= 2*dot(detpart,Jastrow);
    }

    for (int k=0; k<Nparticles; k++){
        for (int i=0; i<Nparticles; i++){
            for (int j=0; j<Nparticles; j++){
                if (i==k || j==k) continue;
                rki = R.row(k) - R.row(i); rkj = R.row(k) - R.row(j); Rki = length(rki); Rkj = length(rkj);
                a = (particles[i].ms != particles[j].ms) ? 0.25: 0.5;
                T += dot(rki, rkj) / (Rki*Rkj) * a / pow(1 + beta*Rki,2) * a / pow(1 + beta*Rkj,2);
            }
        if (i == k) continue;
        T += 2*a / (Rki * pow(1+beta*Rki,2)) - 2*a*beta / pow(1+beta*Rki,3);
        }
    }
    return -0.5*T;
}


















