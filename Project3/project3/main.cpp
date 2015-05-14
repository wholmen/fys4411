#include <iostream>
#include <armadillo>
#include <WaveFunctions.h>
#include <VMCSolver.h>

using namespace std;
using namespace arma;

vec ConjugateGradient(mat A, vec b, vec x0);

class WaveFunction;

class WaveFunction{
private:
    double GTO(double r2, double ksi);
    int i;
    double dist;
public:
    double Helium(mat R);

};

double WaveFunction::GTO(double r2, double ksi){
    return 2*ksi/pow(pi,0.75) * exp(-ksi*r2);
}

double WaveFunction::Helium(mat R){
    double ksi1(13.6267), ksi2(1.99935), ksi3(0.382993);
    double d1(0.1752300), d2(0.8934830), d3(1.00000000);

    double phi1s = 0;


    for (i=0; i<2; i++){
        dist = 0;
        for (int j=0; j<3; j++) dist += R(i,j)*R(i,j);
        phi1s += d1*GTO(dist, ksi1) + d2*GTO(dist, ksi2) + d3 * GTO(dist,ksi3);
    }
    return phi1s;
}


int main()
{

    WaveFunctions helium = WaveFunctions(1.68,0.3,2);
    VMCSolver solve = VMCSolver(helium);
    solve.AnalyticalEnergy = true;
    solve.MCintegration(100000);
    cout << solve.Energy << endl;

    WaveFunctions helium1 = WaveFunctions(1.68,0.3,2);
    VMCSolver solve1 = VMCSolver(helium1);
    solve1.AnalyticalEnergy = true;
    solve1.ImportanceSampling(0.005,100000);
    cout << solve1.Energy;

    /*
    GaussianOrbitals helium2 = GaussianOrbitals(2);
    VMCSolver solve2 = VMCSolver(helium2);
    solve2.AnalyticalEnergy = false;
    solve2.FindStepLength();
    solve2.MCintegration(10000);
    cout << solve2.Energy;
    */
    //solve.AnalyticalEnergy = false;
    //solve.MCintegration(10000);
    //cout << solve.Energy;

    return 0;
}





vec ConjugateGradient(mat A, vec b, vec x0){
    int dim = x0.n_elem;
    vec v(dim), r(dim), z(dim), x(dim);
    double c, t, d, tolerance(1e-4);

    x = x0;
    r = b - A*x;
    v = r;
    c = dot(r,r);

    int i=0; int max_iter = dim;
    while (i <= max_iter){
        z = A*v;
        t = c/dot(v,z);
        x = x + t*v;
        r = r-t*z;
        d = dot(r,r);
        if (sqrt(d) < tolerance) {break;}
        v = r + (d/c)*v;
        c = d;
        i++;
    }
    return x;
}

