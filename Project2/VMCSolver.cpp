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
            Elocal = LocalEnergy(Rnew);
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
    T = T * 0.5 / (h*h) / Psi;

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

