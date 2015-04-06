#include <Beryllium.h>

double Beryllium::phi1s(double r){
    return exp(-alpha*r);
}

double Beryllium::phi2s(double r){
    return (1-alpha*r/2.0)*exp(-alpha*r/2.0);
}

double Beryllium::corrolation(mat R){
    double corr = 0.0; double dist; double a;
    for (int i=0; i<Nparticles; i++){
        for (int k=0; k<i; k++){
            dist = 0.0;
            for (int j=0; j<Ndimensions; j++){
                dist += (R(i,j) - R(k,j))*(R(i,j) - R(k,j));
            }
            dist = sqrt(dist);

            // Testing if spin is unequal or equal. Using fact that i > k.
            if (i == 2 or i == 3){
                if (k == 0 or k == 1){
                    a = 0.5;
                }
                else {
                    a = 0.25;
                }
            }
            else {
                a = 0.25;
            }

            corr *= exp(a*dist / (2*(1+beta*dist)));
        }
    }
    return corr;
}

double Beryllium::WaveFunction(mat R){
    vec arguments = zeros<vec>(Nparticles);

    for (int i=0; i<Nparticles; i++){
        for (int j=0; j<Ndimensions; j++){
            arguments(i) += R(i,j)*R(i,j);
        }
        arguments(i) = sqrt(arguments(i));
    }
    double wavefunction = ( phi1s(arguments(0)) * phi2s(arguments(1)) - phi1s(arguments(1)) * phi2s(arguments(0)) )
                         *( phi1s(arguments(2)) * phi2s(arguments(3)) - phi1s(arguments(3)) * phi2s(arguments(2)) );

    wavefunction = wavefunction * corrolation(R);

    return wavefunction;
}


double Beryllium::LocalEnergy(mat R){
    double r1, r2, r3, r4, r12, r13, r14, r23, r24, r34;

    for (int j=0; j<Ndimensions; j++){
        r1 += R(0,j)*R(0,j);
        r2 += R(1,j)*R(1,j);
        r3 += R(2,j)*R(2,j);
        r4 += R(3,j)*R(3,j);

        r12 += (R(1,j)-R(0,j)) * (R(1,j)-R(0,j));
        r13 += (R(2,j)-R(0,j)) * (R(2,j)-R(0,j));
        r14 += (R(3,j)-R(0,j)) * (R(3,j)-R(0,j));
        r23 += (R(2,j)-R(1,j)) * (R(2,j)-R(1,j));
        r24 += (R(3,j)-R(1,j)) * (R(3,j)-R(1,j));
        r34 += (R(3,j)-R(2,j)) * (R(3,j)-R(2,j));

    }
    r1 = sqrt(r1); r2 = sqrt(r2); r3 = sqrt(r3); r4 = sqrt(r4);

    r12 = sqrt(r12);
    r13 = sqrt(r13);
    r14 = sqrt(r14);
    r23 = sqrt(r23);
    r24 = sqrt(r24);
    r34 = sqrt(r34);

    return (-Z*pow(r1, 2)*r12*r13*r14*r2*r23*r24*r3*exp((9.0L/2.0)*alpha*(r1 + r2 + r3 + r4)) -
            Z*pow(r1, 2)*r12*r13*r14*r2*r23*r24*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) -
            Z*pow(r1, 2)*r12*r13*r14*r23*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) -
            Z*r1*r12*r13*r14*r2*r23*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) +
            alpha*pow(r1, 2)*r12*r13*r14*r2*r23*r24*r3*((alpha*r1 - 2)*exp((1.0/2.0)*alpha*(2*r1 + r2)) +
            (-alpha*r2 + 2)*exp((1.0/2.0)*alpha*(r1 + 2*r2)))*(-0.03125*alpha*r4*(4*(alpha*r3 - 2)*exp((1.0/2.0)*alpha*(2*r3 + r4)) +
            (-alpha*r4 + 6)*exp((1.0/2.0)*alpha*(r3 + 2*r4))) + 0.25*(alpha*r3 - 2)*exp((1.0/2.0)*alpha*(2*r3 + r4)) +
            0.125*(-alpha*r4 + 4)*exp((1.0/2.0)*alpha*(r3 + 2*r4)))*exp(3*alpha*(r1 + r2 + r3 + r4) +
            (1.0/2.0)*r12*(beta*r12 + 1) + (1.0/2.0)*r13*(beta*r13 + 1) + (1.0/2.0)*r14*(beta*r14 + 1) + r23*(beta*r23 + 1) +
            (1.0/2.0)*r24*(beta*r24 + 1)) + alpha*pow(r1, 2)*r12*r13*r14*r2*r23*r24*r4*((alpha*r1 - 2)*exp((1.0/2.0)*alpha*(2*r1 + r2)) +
            (-alpha*r2 + 2)*exp((1.0/2.0)*alpha*(r1 + 2*r2)))*(0.03125*alpha*r3*((-alpha*r3 + 6)*exp((1.0/2.0)*alpha*(2*r3 + r4)) +
            4*(alpha*r4 - 2)*exp((1.0/2.0)*alpha*(r3 + 2*r4))) - 0.125*(-alpha*r3 + 4)*exp((1.0/2.0)*alpha*(2*r3 + r4)) -
            0.25*(alpha*r4 - 2)*exp((1.0/2.0)*alpha*(r3 + 2*r4)))*exp(3*alpha*(r1 + r2 + r3 + r4) + (1.0/2.0)*r12*(beta*r12 + 1) +
            (1.0/2.0)*r13*(beta*r13 + 1) + (1.0/2.0)*r14*(beta*r14 + 1) + r23*(beta*r23 + 1) + (1.0/2.0)*r24*(beta*r24 + 1)) +
            alpha*r12*r13*r14*r2*r23*r24*r3*r4*((alpha*r3 - 2)*exp((1.0/2.0)*alpha*(2*r3 + r4)) + (-alpha*r4 + 2)*exp((1.0/2.0)*alpha*(r3 + 2*r4)))*
            (0.03125*alpha*pow(r1, 2)*((-alpha*r1 + 6)*exp((1.0/2.0)*alpha*(2*r1 + r2)) +
            4*(alpha*r2 - 2)*exp((1.0/2.0)*alpha*(r1 + 2*r2))) - 0.03125*alpha*pow(r2, 2)*(4*(alpha*r1 - 2)*exp((1.0/2.0)*alpha*(2*r1 + r2)) + (-alpha*r2 + 6)*exp((1.0/2.0)*alpha*(r1 + 2*r2))) -
            0.125*r1*((-alpha*r1 + 4)*exp((1.0/2.0)*alpha*(2*r1 + r2)) + 2*(alpha*r2 - 2)*exp((1.0/2.0)*alpha*(r1 + 2*r2))) + 0.125*r2*(2*(alpha*r1 - 2)*exp((1.0/2.0)*alpha*(2*r1 + r2)) +
            (-alpha*r2 + 4)*exp((1.0/2.0)*alpha*(r1 + 2*r2))))*exp(3*alpha*(r1 + r2 + r3 + r4) + (1.0/2.0)*r12*(beta*r12 + 1) + (1.0/2.0)*r13*(beta*r13 + 1) +
            (1.0/2.0)*r14*(beta*r14 + 1) + r23*(beta*r23 + 1) + (1.0/2.0)*r24*(beta*r24 + 1)) + pow(r1, 2)*r12*r13*r14*r2*pow(r23, 2)*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) +
            pow(r1, 2)*r12*r13*r14*r2*r23*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) + pow(r1, 2)*r12*r13*r14*r2*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) +
            pow(r1, 2)*r12*r13*r2*r23*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) + pow(r1, 2)*r12*r14*r2*r23*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)) +
            pow(r1, 2)*r13*r14*r2*r23*r24*r3*r4*exp((9.0/2.0)*alpha*(r1 + r2 + r3 + r4)))*exp(-9.0/2.0*alpha*(r1 + r2 + r3 + r4))/(pow(r1, 2)*r12*r13*r14*r2*r23*r24*r3*r4);
}


