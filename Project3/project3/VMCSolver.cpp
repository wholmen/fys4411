#include <VMCSolver.h>

void VMCSolver::FindStepLength(double b, double a){
    step = 0.2;
    MCintegration_FindVariables(1000, b, a);
    while (AcceptRate > 0.51 or AcceptRate < 0.49){
        step += 0.05;
        MCintegration_FindVariables(1000, b, a);
        if (step > 2.5){break;}
    }
}
void VMCSolver::MCintegration_FindVariables(int Niterations, double b, double a){
    // A function meant to be used when doing steepest descent and when finding steplength

    // Initializing
    mat Rnew = zeros<mat>(Nparticles,Ndimensions);
    mat Rold = zeros<mat>(Nparticles,Ndimensions);

    atom.beta = b; // Changing beta according to the Steepest descent method.
    atom.alpha = a;

    double PsiNew(0), PsiOld(0), Elocal(0), Etotal(0), AcceptStep(0), ExpectationValue1(0), ExpectationValue2(0), PsiFactor(0);

    // Set initial conditions
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rold(i,j) = step * ( ran2(&idum) - 0.5);
        }
    }
    if (not HFBasis) atom.Initialize_System(Rold);
    Rnew = Rold;

    for (int n=0; n<Niterations; n++){
        PsiOld = (HFBasis) ? atom2.WaveFunction(Rold) : atom.WaveFunction(Rold, false);

        for (int i=0; i<Nparticles; i++){

            for (int j=0; j<Ndimensions; j++){
                Rnew(i,j) = Rold(i,j) + step*( ran2(&idum) - 0.5);
            }
            PsiNew = (HFBasis) ? atom2.WaveFunction(Rnew) : atom.WaveFunction(Rnew, false);
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
            PsiFactor = atom.WaveFunction(Rnew,true) / atom.WaveFunction(Rnew,false);
            ExpectationValue1 += PsiFactor * Elocal;
            ExpectationValue2 += PsiFactor;
            Etotal += Elocal;
        }
    }
    AcceptRate = AcceptStep / Nparticles / Niterations;
    ExpectationValue1 = ExpectationValue1 / Niterations / Nparticles;
    ExpectationValue2 = ExpectationValue2 / Niterations / Nparticles;
    Etotal = Etotal / Niterations / Nparticles;
    DiffEbeta = 2 * (ExpectationValue1 - ExpectationValue2 * Etotal);
}


void VMCSolver::MCintegration(vec &Cycle_Energy, vec &Cycle_Energy2, int ncycle, vec &AllEnergies, mat &AllPositions, int Ncycles, int Niterations, double b, double a){
    // A VMC method integral to compute QM expectation values
    atom.beta = b; atom.alpha = a;
    // Initializing
    mat Rnew = zeros<mat>(Nparticles,Ndimensions);
    mat Rold = zeros<mat>(Nparticles,Ndimensions);

    //for (int ncycle = 0; ncycle < Ncycles; ncycle ++){
    double PsiNew(0), PsiOld(0), Elocal(0), Etotal(0), Esquared(0);
    int e(0);

    // Set initial conditions
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rold(i,j) = step * ( ran2(&idum) - 0.5);
        }
    }
    if (not HFBasis) atom.Initialize_System(Rold);
    Rnew = Rold;

    // Total loop for all the cycles
    for (int n=0; n<Niterations; n++){
        PsiOld = (HFBasis) ? atom2.WaveFunction(Rold) : atom.WaveFunction(Rold, false);

        for (int i=0; i<Nparticles; i++){

            for (int j=0; j<Ndimensions; j++){
                Rnew(i,j) = Rold(i,j) + step*( ran2(&idum) - 0.5);
            }
            PsiNew = (HFBasis) ? atom2.WaveFunction(Rnew) : atom.WaveFunction(Rnew, false);
            if (ran2(&idum) <= PsiNew*PsiNew / (PsiOld*PsiOld) ){
                if (not HFBasis) atom.UpdateInverseSD(i, false);
                for (int j=0; j<Ndimensions; j++){
                    Rold(i,j) = Rnew(i,j);
                    PsiOld = PsiNew;
                }
            }
            else{
                for (int j=0; j<Ndimensions; j++){
                    Rnew(i,j) = Rold(i,j);
                }
            }
            Elocal = (AnalyticalEnergy == true) ? LocalEnergyAnalytical(Rnew) : LocalEnergy(Rnew);
            Etotal += Elocal;
            Esquared += Elocal * Elocal;
            if (ncycle == Ncycles-1){ AllEnergies(e) = Elocal; for (int p=0; p<Nparticles; p++) {AllPositions(e,p) = Rnew(p,0);} e++; }
        }
    }
    Energy = Etotal / (Niterations * Nparticles);
    Esquared = Esquared / (Niterations * Nparticles);

    Cycle_Energy(ncycle) = Energy; Cycle_Energy2(ncycle) = Esquared;
}

void VMCSolver::ImportanceSampling(vec &Cycle_Energy, vec &Cycle_Energy2, int ncycle, vec &AllEnergies, mat &AllPositions, int Ncycles, int Niterations, double b, double a, double DT){
    mat Rold = zeros<mat>(Nparticles,Ndimensions);
    mat Rnew = zeros<mat>(Nparticles,Ndimensions);
    mat QForceOld = zeros<mat>(Nparticles,Ndimensions);
    mat QForceNew = zeros<mat>(Nparticles,Ndimensions);
    double D(0.5), PsiOld(0), PsiNew(0), GreenFunction(0), Elocal, Etotal(0), Esquared(0), dt(DT);

    atom.beta = b; // Changing beta according to the Steepest descent method.
    atom.alpha = a;
    int e(0);
    //for (int ncycle=0; ncycle<Ncycles; ncycle++){
    for (int i=0; i<Nparticles; i++) for (int j=0; j<Ndimensions; j++) Rold(i,j) = GaussianDeviate(&idum) * dt;

    if (not HFBasis) atom.Initialize_System(Rold);
    PsiOld = atom.WaveFunction(Rold,false);
    PsiNew = PsiOld;

    for (int n=0; n<Niterations; n++){
        PsiOld = atom.WaveFunction(Rold,false);
        atom.QuantumForce(Rold, QForceOld);

        // Test a movement
        for (int i=0; i<Nparticles; i++){
            for (int j=0; j<Ndimensions; j++) Rnew(i,j) = Rold(i,j) + GaussianDeviate(&idum)*sqrt(dt) + D * dt * QForceOld(i,j);

            for (int k=0; k<Nparticles; k++){
                if (k == i) continue;
                for (int j=0; j<Ndimensions; j++) Rnew(k,j) = Rold(k,j);
            }

            GreenFunction = 0.0;
            PsiNew = atom.WaveFunction(Rnew,false); // Dup and Ddown are updated.
            if (not HFBasis) atom.UpdateInverseSD(i,true); // Make D-1 to compute Quantum force.
            atom.QuantumForce(Rnew, QForceNew);

            for (int j=0; j<Ndimensions; j++){
                GreenFunction += 0.5*(QForceOld(i,j) + QForceNew(i,j)) * (0.5*D*dt*(QForceOld(i,j)-QForceNew(i,j)) - Rnew(i,j) + Rold(i,j));
            }
            GreenFunction = exp(GreenFunction);

            // Check for move
            if (ran2(&idum) <= GreenFunction * (PsiNew*PsiNew / (PsiOld*PsiOld))){
                PsiOld = PsiNew;
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
            if (ncycle == Ncycles-1) {AllEnergies(e) = Elocal; for (int p=0; p<Nparticles; p++) {AllPositions(e,p) = Rnew(p,0);} e++;}
        }
    }
    Energy = Etotal / (Niterations * Nparticles);
    Esquared = Esquared / (Niterations * Nparticles);
    Cycle_Energy(ncycle) = Energy; Cycle_Energy2(ncycle) = Esquared;
}

double VMCSolver::LocalEnergy(mat r){
    mat Rplus = zeros<mat>(Nparticles, Ndimensions);
    mat Rminus = zeros<mat>(Nparticles,Ndimensions);

    Rplus = Rminus = r;
    double Psi = (HFBasis) ? atom2.WaveFunction(r) : atom.WaveFunction(r,false);

    // Calculating Kinectic energy
    double T = 0;
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rplus(i,j) = Rplus(i,j) + h;
            Rminus(i,j) = Rminus(i,j) - h;

            double PsiMinus = (HFBasis) ? atom2.WaveFunction(Rminus) : atom.WaveFunction(Rminus,false);
            double PsiPlus = (HFBasis) ? atom2.WaveFunction(Rplus) : atom.WaveFunction(Rplus,false);
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
