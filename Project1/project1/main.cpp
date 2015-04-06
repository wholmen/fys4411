#include <iostream>
#include <VMCSolver.h>
#include <armadillo>
#include <complex>
#include <lib.h>
#include <fstream>
#include <Beryllium.h>
#include <Atom.h>

using namespace arma;
using namespace std;

void HeliumAlphaZoom();
void HeliumAlpha();
void HeliumAlphaBeta();

int main()
{

    Helium h = Helium(1.7, 2, 1);
    Atom a = Atom(h);
    VMCSolver solve = VMCSolver(a);
    solve.MCintegration();

    Beryllium b = Beryllium(1.7,2);
    Atom c = Atom(b);
    VMCSolver solve2 = VMCSolver(c);
    solve2.MCintegration();
    cout << solve2.Energy;

    return 0;
}

/*

void HeliumAlphaBeta(){
    int Nalpha = 20;
    vec alpha = zeros<vec>(Nalpha);
    for (int i=0; i<Nalpha; i++){
        alpha(i) = 1.65 + i/80.0;
    }

    int Nbeta = 1;
    vec beta = zeros<vec>(Nbeta);
    for (int i=0; i<Nbeta; i++){
        beta(i) = 1 + i;
    }

    ofstream myfile;
    myfile.open("Helium_atom_wf2.txt");
    for (int i=0; i<Nalpha; i+=1){
        for (int j=0; j<Nbeta; j+=1){
            VMCSolver helium = VMCSolver();
            helium.alpha = alpha(i);
            helium.beta = beta(j);
            helium.WFnumber = 2;
            helium.FindStepLength();

            helium.MCintegration();
            myfile << alpha(i) << " " << beta(j) << " " << helium.step << " " << helium.Energy << " " << helium.Variance << " " << helium.AcceptRate << endl;
        }
    }
    myfile.close();
}

void HeliumAlphaZoom(){
    int Nalpha = 11;
    vec alpha = zeros<vec>(Nalpha);
    for (int i=0; i<Nalpha; i++){
        alpha(i) = 1.6 + i/50.0;
    }

    ofstream myfile; myfile.open("Helium_atom_zoom.txt");
    for (int i=0; i<Nalpha; i++){
        VMCSolver helium = VMCSolver();
        helium.alpha = alpha(i);
        helium.FindStepLength();

        helium.MCintegration();
        myfile << alpha(i) << " " << helium.step << " " << helium.Energy << " " << helium.Variance << " " << helium.AcceptRate << endl;

    }
    myfile.close();
}

void HeliumAlpha(){
    int Nalpha = 21;
    vec alpha = zeros<vec>(Nalpha);
    for (int i=0; i<Nalpha; i++){
        alpha(i) = 1.3 + i/30.0;
    }

    ofstream myfile; myfile.open("Helium_atom.txt");
    for (int i=0; i<Nalpha; i++){
        VMCSolver helium = VMCSolver();
        helium.alpha = alpha(i);
        helium.FindStepLength();

        helium.MCintegration();
        myfile << alpha(i) << " " << helium.step << " " << helium.Energy << " " << helium.Variance << " " << helium.AcceptRate << endl;

    }
    myfile.close();
}


*/
