#include <VMCSolver.h>


void VMCSolver::FindStepLength(){
    step = 1.2;
    Nstep = 10000;
    MCintegration();
    while (AcceptRate > 0.51 or AcceptRate < 0.49){
        step += 0.05*ran2(&idum);
        MCintegration();
    }
    Nstep = 1000000;
}

void VMCSolver::MCintegration(){
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
            Rold(i,j) = step * ( ran2(&idum) - 0.5 );
        }
    }

    Rnew = Rold;

    // Total loop for all the cycles
    for (int n=0; n<Nstep; n++){
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
            Elocal = atom.LocalEnergy(Rnew);
            Etotal += Elocal;
            Esquared += Elocal * Elocal;
        }
    }
    Energy = Etotal / (Nstep * Nparticles);
    Esquared = Esquared / (Nstep * Nparticles);
    Variance = - Energy*Energy + Esquared;
    AcceptRate = AcceptStep / (Nstep * Nparticles);
}

void VMCSolver::ImportanceSampling(){
    mat Rold = zeros<mat>(Nparticles,Ndimensions);
    mat Rnew = zeros<mat>(Nparticles,Ndimensions);
    mat QForceOld = zeros<mat>(Nparticles,Ndimensions);
    mat QForceNew = zeros<mat>(Nparticles,Ndimensions);
    double D = 0.5;
    double PsiOld = 0;
    double PsiNew = 0;
    double GreenFunction;
    double Elocal; double Etotal = 0; double Esquared = 0;
    double dt = 0.05;


    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rold(i,j) = GaussianDeviate(&idum) * dt;
        }
    }

    PsiOld = atom.WaveFunction(Rold);
    PsiNew = PsiOld;

    for (int n=0; n<Nstep; n++){
        PsiOld = atom.WaveFunction(Rold);
        atom.QuantumForce(Rold, QForceOld);

        // Test a movement
        for (int i=0; i<Nparticles; i++){
            GreenFunction = 0.0;
            for (int j=0; j<Ndimensions; j++){
                Rnew(i,j) = Rold(i,j) + GaussianDeviate(&idum)*sqrt(dt) + D * dt * QForceOld(i,j);
            }
            for (int k=0; k<Nparticles; k++){
                if (k != i){
                    for (int j=0; j<Ndimensions; j++){
                        Rnew(k,j) = Rold(k,j);
                    }
                }
            }
            GreenFunction = 0.0;
            PsiNew = atom.WaveFunction(Rnew);
            atom.QuantumForce(Rnew, QForceNew);

            for (int j=0; j<Ndimensions; j++){
                GreenFunction += 0.5*(QForceOld(i,j) + QForceNew(i,j)) * (0.5*D*dt*(QForceOld(i,j)-QForceNew(i,j)) - Rnew(i,j) + Rold(i,j));
                GreenFunction = exp(GreenFunction);
            }

            // Check for move
            if (ran2(&idum) <= GreenFunction * (PsiNew*PsiNew / (PsiOld*PsiOld))){
                PsiOld = PsiNew;
                for (int j=0; j<Ndimensions; j++){
                    Rold(i,j) = Rnew(i,j);
                    QForceOld(i,j) = QForceNew(i,j);
                }
            }
            else {
                PsiNew = PsiOld;
                for (int j=0; j<Ndimensions; j++){
                    Rnew(i,j) = Rold(i,j);
                    QForceNew(i,j) = QForceOld(i,j);
                }
            }
            Elocal = atom.LocalEnergy(Rnew);
            Etotal += Elocal;
            Esquared += Elocal * Elocal;
        }
    }
    Energy = Etotal / (Nstep * Nparticles);
    Esquared = Esquared / (Nstep * Nparticles);
    Variance = - Energy*Energy + Esquared;
}




