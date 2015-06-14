#include <GaussianOrbitals.h>

double GaussianOrbitals::factorial(int a){
    double fact(1.0);
    for (int i=a; i>0; i--) fact *= i;
    return fact;
}

double GaussianOrbitals::N(double a, int i, int j, int k){
    return pow(2*a/pi,0.75) * sqrt( pow(4.0*a,i+j+k) / ( factorial(factorial(2*i-1)) * factorial(factorial(2*j-1)) * factorial(factorial(2*k-1)) ));
}

double GaussianOrbitals::Xi(rowvec r, int i, int j, int k, double a, double c){
    double r2 = dot(r,r);
    double xi = c * N(a,i,j,k) * pow(r(0),i) * pow(r(1),j) * pow(r(2),k) * exp(-a*r2);
    return xi;
}

rowvec GaussianOrbitals::XiD(rowvec r, int i, int j, int k, double a, double c){
    rowvec xi(3);
    double r2 = dot(r,r); double x = r(0); double y = r(1); double z = r(2);

    xi(0) = N(a,i,j,k)*c*pow(y,j)*pow(z,k)*(-2*a*pow(x,i+2) + i*pow(x, i))*exp(-a*r2) / x;
    xi(1) = N(a,i,j,k)*c*pow(x,i)*pow(z,k)*(-2*a*pow(y,j+2) + j*pow(y, j))*exp(-a*r2) / y;
    xi(2) = N(a,i,j,k)*c*pow(x,i)*pow(y,j)*(-2*a*pow(z,k+2) + k*pow(z, k))*exp(-a*r2) / z;

    return xi;
}

double GaussianOrbitals::Xi2D(rowvec r, int i, int j, int k, double a, double c){
    double r2 = dot(r,r); double x = r(0); double y = r(1); double z = r(2);
    double xi = N(a,i,j,k)*c*pow(x,i-2)*pow(y,j-2)*pow(z,k-2) * (4*pow(a,2)*pow(x,4)*pow(y,2)*pow(z,2)
         + 4*pow(a,2)*pow(x,2)*pow(y,4)*pow(z,2) + 4*pow(a,2)*pow(x,2)*pow(y,2)*pow(z,4)
         - 4*a*i*pow(x,2)*pow(y,2)*pow(z,2) - 4*a*j*pow(x,2)*pow(y,2)*pow(z,2)
         - 4*a*k*pow(x,2)*pow(y,2)*pow(z,2) - 6*a*pow(x,2)*pow(y,2)*pow(z,2) + pow(i,2)*pow(y,2)*pow(z,2)
         - i*pow(y,2)*pow(z,2) + pow(j,2)*pow(x,2)*pow(z,2) - j*pow(x,2)*pow(z,2)
         + pow(k,2)*pow(x,2)*pow(y,2) - k*pow(x,2)*pow(y,2)) * exp(-a*r2);
    return xi;
}

double GaussianOrbitals::psi(rowvec r, int row){
    double wavefunc(0.0);
    if (Nparticles == 2){
        double k1 = 0.4579; double k2 = 0.6573; double a1 = 13.6267; double a2 = 1.99935; double a3 = 0.382993; double c1 = 0.17523; double c2 = 0.893483; double c3 = 1.0;
        wavefunc = k1*Xi(r,0,0,0,a1,c1) + k2*(Xi(r,0,0,0,a2,c2) + Xi(r,0,0,0,a3,c3) );
    }
    else if (Nparticles == 4){
        double k1, k2, k3, k4, k5, k6, k7, k8, k9;
        if (row == 0) {k1 = -9.9281e-01; k2=-7.6425e-02; k3=2.8727e-02; k4=1.2898e-16; k5=-2.3257e-19; k6=5.6097e-19; k7=1.2016e-16; k8=-4.6874e-19; k9=1.1319e-18; }
        if (row == 1) {k1 = -2.1571e-01; k2 = 2.2934e-01; k3 = 8.2235e-01; k4=5.1721e-16; k5=4.5670e-18; k6=-1.1040e-17; k7=8.5306e-16; k8=7.0721e-18; k9=-1.7060e-17; }
        if (row != 0 && row != 1) cout << "wrong row";
        double a1(71.8876), a2(10.7289), a3(2.22205), a4(1.29548), a5(0.268881), a6(0.07735), a7(1.29548), a8(0.268881), a9(0.07735);
        double c1(0.0644263), c2(0.366096), c3(0.6959340), c4(-0.4210640), c5(1.2240700), c6(1.0), c7(0.2051320), c8(0.8825280), c9(1.0);

        wavefunc = k1*(Xi(r,0,0,0,a1,c1) + Xi(r,0,0,0,a2,c2)  + Xi(r,0,0,0,a3,c3) )
                  +k2*(Xi(r,0,0,0,a4,c4) + Xi(r,0,0,0,a5,c5)) + k3*Xi(r,0,0,0,a6,c6)
                  +k4*(Xi(r,1,0,0,a7,c7) + Xi(r,1,0,0,a8,c8)) + k5*Xi(r,1,0,0,a9,c9)
                  +k6*(Xi(r,0,1,0,a7,c7) + Xi(r,0,1,0,a8,c8)) + k7*Xi(r,0,1,0,a9,c9)
                  +k8*(Xi(r,0,0,1,a7,c7) + Xi(r,0,0,1,a8,c8)) + k9*Xi(r,0,0,1,a9,c9);
    }
    else if (Nparticles == 10){
        double k1, k2, k3, k4, k5, k6, k7, k8, k9;
        if (row == 0) {k1=-9.8077e-01; k2=-9.3714e-02; k3= 2.2863e-02; k4=-9.9519e-19; k5=-1.2125e-18; k6=-4.1800e-19; k7=-1.6696e-19; k8= 1.2125e-18; k9= 3.8779e-19; }
        if (row == 1) {k1=-2.6062e-01; k2= 2.5858e-01; k3= 8.1619e-01; k4=-5.6186e-18; k5=-2.8615e-16; k6= 4.6199e-17; k7=-4.2405e-18; k8=-2.9426e-16; k9= 5.0519e-17; }
        if (row == 2) {k1= 1.1596e-16; k2=-2.0106e-16; k3=-3.2361e-16; k4= 2.7155e-02; k5=-5.6207e-01; k6= 9.1139e-03; k7= 2.8890e-02; k8=-5.9797e-01; k9= 9.6959e-03; }
        if (row == 3) {k1=-8.3716e-18; k2=-9.7173e-17; k3= 1.3237e-16; k4=-4.0320e-01; k5=-2.5833e-02; k6=-3.9180e-01; k7=-4.2895e-01; k8=-2.7482e-02; k9=-4.1683e-01; }
        if (row == 4) {k1=-1.9554e-17; k2=-7.3738e-17; k3= 1.5789e-16; k4= 3.9171e-01; k5= 1.2375e-02; k6=-4.0392e-01; k7= 4.1673e-01; k8= 1.3166e-02; k9=-4.2972e-01; }

        double a1(515.724), a2(77.6538), a3(16.8136), a4(12.483), a5(2.66451), a6(0.60625), a7(12.483), a8(2.66451), a9(0.60625);
        double c1(0.058143), c2(0.347951), c3(0.710714), c4(-0.409922), c5(1.22431), c6(1.0), c7( 0.24746), c8(0.851743), c9(1.0);

        wavefunc = k1*(Xi(r,0,0,0,a1,c1) + Xi(r,0,0,0,a2,c2)  + Xi(r,0,0,0,a3,c3) )
                  +k2*(Xi(r,0,0,0,a4,c4) + Xi(r,0,0,0,a5,c5)) + k3*Xi(r,0,0,0,a6,c6)
                  +k4*(Xi(r,1,0,0,a7,c7) + Xi(r,1,0,0,a8,c8)) + k5*Xi(r,1,0,0,a9,c9)
                  +k6*(Xi(r,0,1,0,a7,c7) + Xi(r,0,1,0,a8,c8)) + k7*Xi(r,0,1,0,a9,c9)
                  +k8*(Xi(r,0,0,1,a7,c7) + Xi(r,0,0,1,a8,c8)) + k9*Xi(r,0,0,1,a9,c9);
    }
    return wavefunc;
}

rowvec GaussianOrbitals::dpsi(rowvec r, int row){
    rowvec wavefunc(3);
    if (Nparticles == 2){
        double k1 = 0.4579; double k2 = 0.6573; double a1 = 13.6267; double a2 = 1.99935; double a3 = 0.382993; double c1 = 0.17523; double c2 = 0.893483; double c3 = 1.0;
        wavefunc = k1*XiD(r,0,0,0,a1,c1) + k2*(XiD(r,0,0,0,a2,c2) + XiD(r,0,0,0,a3,c3) );
    }
    else if (Nparticles == 4){
        double k1, k2, k3, k4, k5, k6, k7, k8, k9;
        if (row == 0) {k1 = -9.9281e-01; k2=-7.6425e-02; k3=2.8727e-02; k4=1.2898e-16; k5=-2.3257e-19; k6=5.6097e-19; k7=1.2016e-16; k8=-4.6874e-19; k9=1.1319e-18; }
        if (row == 1) {k1 = -2.1571e-01; k2 = 2.2934e-01; k3 = 8.2235e-01; k4=5.1721e-16; k5=4.5670e-18; k6=-1.1040e-17; k7=8.5306e-16; k8=7.0721e-18; k9=-1.7060e-17; }

        double a1(71.8876), a2(10.7289), a3(2.22205), a4(1.29548), a5(0.268881), a6(0.07735), a7(1.29548), a8(0.268881), a9(0.07735);
        double c1(0.0644263), c2(0.366096), c3(0.6959340), c4(-0.4210640), c5(1.2240700), c6(1.0), c7(0.2051320), c8(0.8825280), c9(1.0);

        wavefunc = k1*(XiD(r,0,0,0,a1,c1) + XiD(r,0,0,0,a2,c2)  + XiD(r,0,0,0,a3,c3) )
                  +k2*(XiD(r,0,0,0,a4,c4) + XiD(r,0,0,0,a5,c5)) + k3*XiD(r,0,0,0,a6,c6)
                  +k4*(XiD(r,1,0,0,a7,c7) + XiD(r,1,0,0,a8,c8)) + k5*XiD(r,1,0,0,a9,c9)
                  +k6*(XiD(r,0,1,0,a7,c7) + XiD(r,0,1,0,a8,c8)) + k7*XiD(r,0,1,0,a9,c9)
                  +k8*(XiD(r,0,0,1,a7,c7) + XiD(r,0,0,1,a8,c8)) + k9*XiD(r,0,0,1,a9,c9);
    }
    else if (Nparticles == 10){
            double k1, k2, k3, k4, k5, k6, k7, k8, k9;
            if (row == 0) {k1=-9.8077e-01; k2=-9.3714e-02; k3= 2.2863e-02; k4=-9.9519e-19; k5=-1.2125e-18; k6=-4.1800e-19; k7=-1.6696e-19; k8= 1.2125e-18; k9= 3.8779e-19; }
            if (row == 1) {k1=-2.6062e-01; k2= 2.5858e-01; k3= 8.1619e-01; k4=-5.6186e-18; k5=-2.8615e-16; k6= 4.6199e-17; k7=-4.2405e-18; k8=-2.9426e-16; k9= 5.0519e-17; }
            if (row == 2) {k1= 1.1596e-16; k2=-2.0106e-16; k3=-3.2361e-16; k4= 2.7155e-02; k5=-5.6207e-01; k6= 9.1139e-03; k7= 2.8890e-02; k8=-5.9797e-01; k9= 9.6959e-03; }
            if (row == 3) {k1=-8.3716e-18; k2=-9.7173e-17; k3= 1.3237e-16; k4=-4.0320e-01; k5=-2.5833e-02; k6=-3.9180e-01; k7=-4.2895e-01; k8=-2.7482e-02; k9=-4.1683e-01; }
            if (row == 4) {k1=-1.9554e-17; k2=-7.3738e-17; k3= 1.5789e-16; k4= 3.9171e-01; k5= 1.2375e-02; k6=-4.0392e-01; k7= 4.1673e-01; k8= 1.3166e-02; k9=-4.2972e-01; }

            double a1(515.724), a2(77.6538), a3(16.8136), a4(12.483), a5(2.66451), a6(0.60625), a7(12.483), a8(2.66451), a9(0.60625);
            double c1(0.058143), c2(0.347951), c3(0.710714), c4(-0.409922), c5(1.22431), c6(1.0), c7( 0.24746), c8(0.851743), c9(1.0);

            wavefunc = k1*(XiD(r,0,0,0,a1,c1) + XiD(r,0,0,0,a2,c2)  + XiD(r,0,0,0,a3,c3) )
                      +k2*(XiD(r,0,0,0,a4,c4) + XiD(r,0,0,0,a5,c5)) + k3*XiD(r,0,0,0,a6,c6)
                      +k4*(XiD(r,1,0,0,a7,c7) + XiD(r,1,0,0,a8,c8)) + k5*XiD(r,1,0,0,a9,c9)
                      +k6*(XiD(r,0,1,0,a7,c7) + XiD(r,0,1,0,a8,c8)) + k7*XiD(r,0,1,0,a9,c9)
                      +k8*(XiD(r,0,0,1,a7,c7) + XiD(r,0,0,1,a8,c8)) + k9*XiD(r,0,0,1,a9,c9);
    }
    return wavefunc;
}

double GaussianOrbitals::d2psi(rowvec r, int row){
    double wavefunc(0.0);
    if (Nparticles == 2){
        double k1 = 0.4579; double k2 = 0.6573; double a1 = 13.6267; double a2 = 1.99935; double a3 = 0.382993; double c1 = 0.17523; double c2 = 0.893483; double c3 = 1.0;
        wavefunc = k1*Xi2D(r,0,0,0,a1,c1) + k2*(Xi2D(r,0,0,0,a2,c2) + Xi2D(r,0,0,0,a3,c3) );
    }
    else if (Nparticles == 4){
        double k1, k2, k3, k4, k5, k6, k7, k8, k9;
        if (row == 0) {k1=-9.9281e-01; k2=-7.6425e-02; k3=2.8727e-02; k4=1.2898e-16; k5=-2.3257e-19; k6= 5.6097e-19; k7=1.2016e-16; k8=-4.6874e-19; k9= 1.1319e-18;}
        if (row == 1) {k1=-2.1571e-01; k2= 2.2934e-01; k3=8.2235e-01; k4=5.1721e-16; k5= 4.5670e-18; k6=-1.1040e-17; k7=8.5306e-16; k8= 7.0721e-18; k9=-1.7060e-17;}

        double a1(71.8876), a2(10.7289), a3(2.22205), a4(1.29548), a5(0.268881), a6(0.07735), a7(1.29548), a8(0.268881), a9(0.07735);
        double c1(0.0644263), c2(0.366096), c3(0.6959340), c4(-0.4210640), c5(1.2240700), c6(1.0), c7(0.2051320), c8(0.8825280), c9(1.0);

        wavefunc = k1*(Xi2D(r,0,0,0,a1,c1) + Xi2D(r,0,0,0,a2,c2)  + Xi2D(r,0,0,0,a3,c3) )
                  +k2*(Xi2D(r,0,0,0,a4,c4) + Xi2D(r,0,0,0,a5,c5)) + k3*Xi2D(r,0,0,0,a6,c6)
                  +k4*(Xi2D(r,1,0,0,a7,c7) + Xi2D(r,1,0,0,a8,c8)) + k5*Xi2D(r,1,0,0,a9,c9)
                  +k6*(Xi2D(r,0,1,0,a7,c7) + Xi2D(r,0,1,0,a8,c8)) + k7*Xi2D(r,0,1,0,a9,c9)
                  +k8*(Xi2D(r,0,0,1,a7,c7) + Xi2D(r,0,0,1,a8,c8)) + k9*Xi2D(r,0,0,1,a9,c9);
    }
    else if (Nparticles == 10){
            double k1, k2, k3, k4, k5, k6, k7, k8, k9;
            if (row == 0) {k1=-9.8077e-01; k2=-9.3714e-02; k3= 2.2863e-02; k4=-9.9519e-19; k5=-1.2125e-18; k6=-4.1800e-19; k7=-1.6696e-19; k8= 1.2125e-18; k9= 3.8779e-19; }
            if (row == 1) {k1=-2.6062e-01; k2= 2.5858e-01; k3= 8.1619e-01; k4=-5.6186e-18; k5=-2.8615e-16; k6= 4.6199e-17; k7=-4.2405e-18; k8=-2.9426e-16; k9= 5.0519e-17; }
            if (row == 2) {k1= 1.1596e-16; k2=-2.0106e-16; k3=-3.2361e-16; k4= 2.7155e-02; k5=-5.6207e-01; k6= 9.1139e-03; k7= 2.8890e-02; k8=-5.9797e-01; k9= 9.6959e-03; }
            if (row == 3) {k1=-8.3716e-18; k2=-9.7173e-17; k3= 1.3237e-16; k4=-4.0320e-01; k5=-2.5833e-02; k6=-3.9180e-01; k7=-4.2895e-01; k8=-2.7482e-02; k9=-4.1683e-01; }
            if (row == 4) {k1=-1.9554e-17; k2=-7.3738e-17; k3= 1.5789e-16; k4= 3.9171e-01; k5= 1.2375e-02; k6=-4.0392e-01; k7= 4.1673e-01; k8= 1.3166e-02; k9=-4.2972e-01; }

            double a1(515.724), a2(77.6538), a3(16.8136), a4(12.483), a5(2.66451), a6(0.60625), a7(12.483), a8(2.66451), a9(0.60625);
            double c1(0.058143), c2(0.347951), c3(0.710714), c4(-0.409922), c5(1.22431), c6(1.0), c7( 0.24746), c8(0.851743), c9(1.0);

            wavefunc = k1*(Xi2D(r,0,0,0,a1,c1) + Xi2D(r,0,0,0,a2,c2)  + Xi2D(r,0,0,0,a3,c3) )
                      +k2*(Xi2D(r,0,0,0,a4,c4) + Xi2D(r,0,0,0,a5,c5)) + k3*Xi2D(r,0,0,0,a6,c6)
                      +k4*(Xi2D(r,1,0,0,a7,c7) + Xi2D(r,1,0,0,a8,c8)) + k5*Xi2D(r,1,0,0,a9,c9)
                      +k6*(Xi2D(r,0,1,0,a7,c7) + Xi2D(r,0,1,0,a8,c8)) + k7*Xi2D(r,0,1,0,a9,c9)
                      +k8*(Xi2D(r,0,0,1,a7,c7) + Xi2D(r,0,0,1,a8,c8)) + k9*Xi2D(r,0,0,1,a9,c9);
    }
    return wavefunc;
}













