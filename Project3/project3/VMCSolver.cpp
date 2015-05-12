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
    atom.Initialize_System(Rold);
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
                atom.UpdateInverseSD(Rnew,i);
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
    T += 2.0 * atom.diff_Slater_diff_Corrolation(R);
    T += atom.diff2_Corrolation(R);
    T = T * -0.5;
    /*
    double r1(0), r2(0), r12(0), r1r2(0), T1, T2, T3;
    for (int i=0;i<Ndimensions;i++){
        r1 += R(0,i)*R(0,i);
        r2 += R(1,i)*R(1,i);
        r12 += pow(R(0,i)-R(1,i),2);
        r1r2 += R(0,i)*R(1,i);
    }
    r1 = sqrt(r1); r2 = sqrt(r2); r12 = sqrt(r12);
    double alpha = 1.7; double beta = 0.3;

    T1 = alpha/r1 + alpha/r2 - alpha*alpha;
    T2 = 1.0/(2*pow(1+beta*r12,2)) * alpha*(r1+r2) / r12 * (1-r1r2/r1*r2);
    T3 = 1.0/(2*pow(1+beta*r12,2)) * ( -1.0/(2*pow(1+beta*r12,2)) - 2/r12 + 2*beta/(1+beta*r12));
    cout << T << " " << T1+T2+T3 << endl;
    T = T1 + T2 + T3;
    */
    //cout << "T1 num: " << atom.diff2_Slater(R)*-0.5 << " Analytical: " << T1 << endl;
    //cout << "T2 num: " << -atom.diff_Slater_diff_Corrolation(R) << " Analytical: " << T2 << endl;
    //cout << "T3 num: " << atom.diff2_Corrolation(R)*-0.5 << " Analytical: " << T3 << endl;

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
