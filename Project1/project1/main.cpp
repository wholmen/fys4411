#include <iostream>
#include <VMCSolver.h>
#include <armadillo>
#include <lib.h>
#include <fstream>
#include <Beryllium.h>
#include <Atom.h>

using namespace arma;
using namespace std;

void Helium_FindAlphaZoom();
void Helium_FindAlpha();
void Helium_FindBeta();
void ImportanceSampling_FindDt();
void FinalHelium();

int main()
{
    //Helium_FindAlpha();
    //Helium_FindAlphaZoom();
    //Helium_FindBeta();
    //ImportanceSampling_FindDt();

    FinalHelium();

    return 0;
}


void FinalHelium(){
    double alpha = 1.75;
    double beta = 0.3;

    Helium h = Helium(alpha,beta,2);
    Atom a = Atom(h);
    VMCSolver solve = VMCSolver(a,false);

    solve.FindStepLength();
    solve.MCintegration();
    cout << solve.Energy << " " << solve.step;
}

void ImportanceSampling_FindDt(){
    Helium h = Helium(1.68, 1, 2);
    Atom a = Atom(h);
    VMCSolver solve = VMCSolver(a,false);
    //solve.MCintegration();

    vec delta_t = zeros<vec>(6);
    delta_t(0) = 1; delta_t(1) = 0.1; delta_t(2) = 0.01; delta_t(3) = 0.001; delta_t(4) = 0.0001; delta_t(5) = 0.00001;
    ofstream myfile;
    myfile.open("Compare_dt_Importance.txt");
    for (int i=0; i<6; i++){
        solve.ImportanceSampling(delta_t(i));
        myfile << delta_t(i) << " " << solve.Energy << endl;
    }
}

void Helium_FindBeta(){
    int Nbeta = 8;
    vec beta = zeros<vec>(Nbeta);
    for (int i=0; i<Nbeta; i++){
        beta(i) = 0.3 + i*0.2;
    }

    int Nalpha = 8;
    vec alpha = zeros<vec>(Nalpha);
    for (int i=0; i<Nalpha; i++){
        alpha(i) = 1.6 + i*0.025;
    }

    ofstream myfile;
    myfile.open("Helium_atom_wf2.txt");
    for (int i=0; i<Nalpha; i++){
        for (int j=0; j<Nbeta; j+=1){
            Helium helium = Helium(alpha(i), beta(j), 2);
            Atom atom = Atom(helium);
            VMCSolver solver = VMCSolver(atom,false);

            solver.FindStepLength();
            solver.MCintegration();
            myfile << alpha(i) << " " << beta(j) << " " << solver.step << " " << solver.Energy << " " << solver.Variance << " " << endl;
        }
    }
    myfile.close();
}

void Helium_FindAlphaZoom(){
    int Nalpha = 11;
    vec alpha = zeros<vec>(Nalpha);
    for (int i=0; i<Nalpha; i++){
        alpha(i) = 1.6 + i/50.0;
    }

    ofstream myfile; myfile.open("Helium_atom_zoom.txt");
    for (int i=0; i<Nalpha; i++){
        Helium helium = Helium(alpha(i),1,1);
        Atom atom = Atom(helium);
        VMCSolver solver = VMCSolver(atom,false);

        solver.FindStepLength();
        solver.MCintegration();
        myfile << alpha(i) << " " << solver.step << " " << solver.Energy << " " << solver.Variance << " " << endl;

    }
    myfile.close();
}

void Helium_FindAlpha(){
    int Nalpha = 31;
    vec alpha = zeros<vec>(Nalpha);
    for (int i=0; i<Nalpha; i++){
        alpha(i) = 0.0 + i/10.0;
    }

    ofstream myfile; myfile.open("Helium_atom.txt");
    for (int i=0; i<Nalpha; i++){
        Helium helium = Helium(alpha(i),1,1);
        Atom atom = Atom(helium);
        VMCSolver solver = VMCSolver(atom,false);

        solver.FindStepLength();
        solver.MCintegration();
        myfile << alpha(i) << " " << solver.step << " " << solver.Energy << " " << solver.Variance << " " << endl;
    }
    myfile.close();
}

