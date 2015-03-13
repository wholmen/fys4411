#include <VMCSolver.h>

VMCSolver::VMCSolver() :
        Nstep(100000),
        Nparticles(2),
        Ndimensions(3),
        charge(2),
        h(0.001),
        sqh(sqrt(h)),
        idum(-1),
        alpha(0.5*charge),
        step(2.55),
        beta(0.5*charge),
        AcceptRate(0),
        WFnumber(1)
    {
    }

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
        PsiOld = WaveFunction(Rold);

        for (int i=0; i<Nparticles; i++){

            for (int j=0; j<Ndimensions; j++){
                Rnew(i,j) = Rold(i,j) + step*( ran2(&idum) - 0.5);
            }
            PsiNew = WaveFunction(Rnew);
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

    PsiOld = WaveFunction(Rold);
    PsiNew = PsiOld;

    for (int n=0; n<Nstep; n++){
        PsiOld = WaveFunction(Rold);
        QuantumForce(Rold, QForceOld);

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
            PsiNew = WaveFunction(Rnew);
            QuantumForce(Rnew, QForceNew);

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
            Elocal = LocalEnergy(Rnew);
            Etotal += Elocal;
            Esquared += Elocal * Elocal;
        }
    }
    Energy = Etotal / (Nstep * Nparticles);
    Esquared = Esquared / (Nstep * Nparticles);
    Variance = - Energy*Energy + Esquared;
}

void VMCSolver::QuantumForce(mat r, mat &Qforce){
    mat Rplus = zeros<mat>(Nparticles,Ndimensions);
    mat Rminus = zeros<mat>(Nparticles,Ndimensions);
    double Psiplus = 0; double Psiminus = 0;

    Rplus = Rminus = r;
    double Psi = WaveFunction(r);

    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rplus(i,j)  += h;
            Rminus(i,j) -= h;
            Psiplus = WaveFunction(Rplus);
            Psiminus = WaveFunction(Rminus);

            Qforce(i,j) = Psiplus - Psiminus;
            Rplus(i,j) = r(i,j); Rminus(i,j) = r(i,j);
        }
    }
    Qforce = Qforce / h / Psi;
}

double VMCSolver::LocalEnergy(mat r){
    mat Rplus = zeros<mat>(Nparticles, Ndimensions);
    mat Rminus = zeros<mat>(Nparticles,Ndimensions);

    Rplus = Rminus = r;
    double Psi = WaveFunction(r);

    // Calculating Kinectic energy
    double T = 0;
    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            Rplus(i,j) = Rplus(i,j) + h;
            Rminus(i,j) = Rminus(i,j) - h;

            double PsiMinus = WaveFunction(Rminus);
            double PsiPlus = WaveFunction(Rplus);
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

double VMCSolver::WaveFunction(mat R){
    if (WFnumber == 1){
        double argument = 0; double r1;
        for (int i=0; i<Nparticles; i++){
            r1 = 0;
            for (int j=0; j<Ndimensions; j++){
                r1 += R(i,j) * R(i,j);
            }
            r1 = sqrt(r1);
            argument += r1;
        }
        return exp(-alpha*argument);
    }
    if (WFnumber == 2){
        double argument = 0; double dargument = 0; double r1; double r12 = 0;
        for (int i=0; i<Nparticles; i++){
            r1 = 0; r12 = 0;

            for (int j=0; j<Ndimensions; j++){
                r1 += R(i,j) * R(i,j);
            }
            argument += sqrt(r1);

            for (int k=i+1; k<Nparticles; k++){
                for (int j=0; j<Ndimensions; j++){
                    r12 += (R(i,j) - R(k,j) )*(R(i,j) - R(k,j) );
                }
                dargument += sqrt(r12);
            }
        }
        return exp(-alpha*argument) * exp(dargument/(2*(1+beta*dargument)));
    }
    else {return false;}
}

