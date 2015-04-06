#include <Helium.h>




void Helium::QuantumForce(mat r, mat &Qforce){
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

double Helium::LocalEnergy(mat r){
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

double Helium::WaveFunction(mat R){
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
