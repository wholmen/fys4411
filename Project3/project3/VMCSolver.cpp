#include <VMCSolver.h>

void VMCSolver::FindStepLength(){
    step = 0.2;
    MCintegration(1000);
    while (AcceptRate > 0.51 or AcceptRate < 0.49){
        step += 0.05;
        MCintegration(1000);
        if (step > 2.5){break;}
    }
    cout << step;
}

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
    if (not HFBasis) atom.Initialize_System(Rold);
    Rnew = Rold;

    // Total loop for all the cycles
    for (int n=0; n<N; n++){
        PsiOld = (HFBasis) ? atom2.WaveFunction(Rold) : atom.WaveFunction(Rold);

        for (int i=0; i<Nparticles; i++){

            for (int j=0; j<Ndimensions; j++){
                Rnew(i,j) = Rold(i,j) + step*( ran2(&idum) - 0.5);
            }
            PsiNew = (HFBasis) ? atom2.WaveFunction(Rnew) : atom.WaveFunction(Rnew);
            if (ran2(&idum) <= PsiNew*PsiNew / (PsiOld*PsiOld) ){
                if (not HFBasis) atom.UpdateInverseSD(i, false);
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

void VMCSolver::ImportanceSampling(double DT, int Nstep){
    mat Rold = zeros<mat>(Nparticles,Ndimensions);
    mat Rnew = zeros<mat>(Nparticles,Ndimensions);
    mat QForceOld = zeros<mat>(Nparticles,Ndimensions);
    mat QForceNew = zeros<mat>(Nparticles,Ndimensions);
    double D(0.5), PsiOld(0), PsiNew(0), GreenFunction(0), Elocal, Etotal(0), Esquared(0), dt(DT);

    for (int i=0; i<Nparticles; i++) for (int j=0; j<Ndimensions; j++) Rold(i,j) = GaussianDeviate(&idum) * dt;

    if (not HFBasis) atom.Initialize_System(Rold);
    PsiOld = atom.WaveFunction(Rold);
    PsiNew = PsiOld;

    for (int n=0; n<Nstep; n++){
        PsiOld = atom.WaveFunction(Rold);
        atom.QuantumForce(Rold, QForceOld);

        // Test a movement
        for (int i=0; i<Nparticles; i++){
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
            PsiNew = atom.WaveFunction(Rnew); // Dup and Ddown are updated.
            if (not HFBasis) atom.UpdateInverseSD(i,true); // Make D-1 to compute Quantum force.
            atom.QuantumForce(Rnew, QForceNew);


            for (int j=0; j<Ndimensions; j++){
                GreenFunction += 0.5*(QForceOld(i,j) + QForceNew(i,j)) * (0.5*D*dt*(QForceOld(i,j)-QForceNew(i,j)) - Rnew(i,j) + Rold(i,j));
            }
            GreenFunction = exp(GreenFunction);

            // Check for move
            if (ran2(&idum) <= GreenFunction * (PsiNew*PsiNew / (PsiOld*PsiOld))){
                PsiOld = PsiNew;
                //if (not HFBasis) atom.UpdateInverseSD(i); // Particle moved. Inverse D-1 must be updated
                if (not HFBasis) atom.UpdateInverseQF(true); // true. Move is accepted. The D-1 created is now stored
                for (int j=0; j<Ndimensions; j++){
                    Rold(i,j) = Rnew(i,j);
                    QForceOld(i,j) = QForceNew(i,j);
                }
            }
            else {
                if (not HFBasis) atom.UpdateInverseQF(false); // false. Move is declined. The D-1 created is tossed.
                for (int j=0; j<Ndimensions; j++){
                    Rnew(i,j) = Rold(i,j);
                    QForceNew(i,j) = QForceOld(i,j);
                }
            }
            Elocal = (AnalyticalEnergy == true) ? LocalEnergyAnalytical(Rnew) : LocalEnergy(Rnew);
            Etotal += Elocal;
            Esquared += Elocal * Elocal;
        }
    }
    Energy = Etotal / (Nstep * Nparticles);
    Esquared = Esquared / (Nstep * Nparticles);
    Variance = - Energy*Energy + Esquared;
}

double VMCSolver::LocalEnergy(mat r){
    mat Rplus = zeros<mat>(Nparticles, Ndimensions);
    mat Rminus = zeros<mat>(Nparticles,Ndimensions);

    Rplus = Rminus = r;
    double Psi = (HFBasis) ? atom2.WaveFunction(r) : atom.WaveFunction(r);

    // Calculating Kinectic energy
    double T = 0;
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rplus(i,j) = Rplus(i,j) + h;
            Rminus(i,j) = Rminus(i,j) - h;

            double PsiMinus = (HFBasis) ? atom2.WaveFunction(Rminus) : atom.WaveFunction(Rminus);
            double PsiPlus = (HFBasis) ? atom2.WaveFunction(Rplus) : atom.WaveFunction(Rplus);
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

    T = (HFBasis) ? atom2.KineticEnergy(R) : atom.KineticEnergy(R);

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

            V += 1.0 / sqrt(TwoParticle);
        }
    }
    return T + V;
}
