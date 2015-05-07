#include <VMCSolver.h>


void VMCSolver::MCintegration(int N){
    // A VMC method integral to compute QM expectation values

    // Initializing
    mat Rnew = zeros<mat>(Nparticles,Ndimensions);
    mat Rold = zeros<mat>(Nparticles,Ndimensions);
    double PsiNew = 0;
    double PsiOld = 0;
    double Elocal = 0;
    double Etotal = 0;
    double Esquared = 0;
    double AcceptStep = 0;

    // Set initial conditions
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rold(i,j) = step * ( ran2(&idum) - 0.5);
        }
    }

    Rnew = Rold;

    // Total loop for all the cycles
    for (int n=0; n<N; n++){
        PsiOld = atom.WaveFunction(Rold);

        for (int i=0; i<Nparticles; i++){

            for (int j=0; j<Ndimensions; j++){
                Rnew(i,j) = Rold(i,j) + step*( ran2(&idum) - 0.5);
            }
            PsiNew = atom.WaveFunction(Rnew);
            if (ran2(&idum) <= PsiNew*PsiNew / (PsiOld*PsiOld) ){
                for (int j=0; j<Ndimensions; j++){
                    Rold(i,j) = Rnew(i,j);
                    PsiOld = PsiNew;
                }
                AcceptStep ++;
            }
            else{
                for (int j=0; j<Ndimensions; j++){
                    Rnew(i,j) = Rold(i,j);
                }
            }
            (AnalyticalEnergy == true) ? Elocal = LocalEnergyAnalytical(Rnew) : Elocal = LocalEnergy(Rnew);
            Etotal += Elocal;
            Esquared += Elocal * Elocal;
        }
    }
    Energy = Etotal / (N * Nparticles);
    Esquared = Esquared / (N * Nparticles);
    Variance = - Energy*Energy + Esquared;
    AcceptRate = AcceptStep / (N * Nparticles);
}


double VMCSolver::LocalEnergy(mat r){
    mat Rplus = zeros<mat>(Nparticles, Ndimensions);
    mat Rminus = zeros<mat>(Nparticles,Ndimensions);

    Rplus = Rminus = r;
    double Psi = atom.WaveFunction(r);

    // Calculating Kinectic energy
    double T = 0;
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rplus(i,j) = Rplus(i,j) + h;
            Rminus(i,j) = Rminus(i,j) - h;

            double PsiMinus = atom.WaveFunction(Rminus);
            double PsiPlus = atom.WaveFunction(Rplus);
            T -= (PsiMinus + PsiPlus - 2*Psi);
            Rplus(i,j) = r(i,j); Rminus(i,j) = r(i,j);
        }
    }
    if (Psi != 0) T = T * 0.5 / (h*h) / Psi;
    else cout << "Error when computing Slater Determinant. WaveFunction return 0. " << endl;

    // Calculate potential energy for single particle
    double V = 0;
    for (int i=0; i<Nparticles; i++){
        double Singleparticle = 0;
        for (int j=0; j<Ndimensions; j++){
            Singleparticle += r(i,j) * r(i,j);
        }
        V -= charge/sqrt(Singleparticle);
    }
    // Calculate potential energy between two particles.
    for (int i=0; i<Nparticles; i++){
        for (int k=i+1; k<Nparticles; k++){
            double TwoParticle = 0;
            for (int j=0; j<Ndimensions; j++){
                double distance = r(k,j) - r(i,j);
                TwoParticle += distance * distance;
            }
            V += 1.0 / sqrt(TwoParticle);
        }
    }

    return V + T;
}

double VMCSolver::LocalEnergyAnalytical(mat R){
    // Calculating Kinetic energy.
    double T(0), V(0), TwoParticle(0), Singleparticle(0);

    T += atom.diff2_Slater(R);
    //T += 2.0 * atom.diff_Slater_diff_Corrolation(R);
    //T += atom.diff2_Corrolation(R);

    T = T * -0.5;

    // Calculate potential energy for single particle
    for (int i=0; i<Nparticles; i++){
        Singleparticle = 0;
        for (int j=0; j<Ndimensions; j++) Singleparticle += R(i,j) * R(i,j);

        V -= charge/sqrt(Singleparticle);
    }
    // Calculate potential energy between two particles.
    for (int i=0; i<Nparticles; i++){
        for (int k=i+1; k<Nparticles; k++){
            TwoParticle = 0;
            for (int j=0; j<Ndimensions; j++) TwoParticle += (R(k,j) - R(i,j)) * (R(k,j) - R(i,j));

            //V += 1.0 / sqrt(TwoParticle);
        }
    }
    return T + V;
}

/*
double VMCSolver::LocalEnergyAnalytical(mat R){
    // return dell^2 psi_D / psi_D + dell^2 psi_C / psi_c + 2*dell Psi_D / Psi_D * dell Psi_C / Psi_C

    // Calculating Kinetic energy
    double E = 0;

    // Calculating the derivatives in the kinetic energy
    double SlaterTerm = 0; double CorrolationTerm = 0; double CrossTerm = 0;
    for (int i=0; i<Nparticles; i++){

        // Double derivative of Slater term, Psi_D
        double ri = 0;
        for (int j=0; j<Ndimensions; j++){
            ri += R(i,j)*R(i,j);
        }
        ri = sqrt(ri);
        SlaterTerm = SlaterTerm - 2*alpha/ri + alpha**2;

        // Single derivatives giving the Cross terms, Psi_C
        double rik = 0;
        for (int k=0; k<Nparticles; k++){
            for (int j=0; j<Ndimensions; j++){
                rik += (R(i,j)-R(k,j))*(R(i,j)-R(k,j));
                rirk += R(i,j)*R(k,j);
            }
            CrossTerm = CrossTerm - 2 * alpha * a / (rik * (1+beta*rik)) * (rirk / ri - ri);
        }
    }
    // Double derivative of corrolation term
    for (int k=0; k<Nparticles; k++){
        for (int j=0; j<Nparticles; j++){
            for (int i=0; i<Nparticles; i++){
                for (int d=0; d<Ndimensions; d++);
                    rk = ;
                    ri =

                CorrolationTerm += (rk - ri)*(rk - rj) / (rki * rkj) * a/(1+beta*rki)/(1+beta*rki) * a/(1+beta*rkj)/(1+beta*rkj);
            }
            CorrolationTerm += 2*a / (rkj*(1+beta*rkj)*(1+beta*rkj)) - 2*alpha*beta/ ( (1+beta*rkj)*(1+beta*rkj) );
        }
    }



    return E;
}

*/
