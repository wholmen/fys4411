#include <iostream>
#include <VMCSolver.h>
#include <armadillo>
#include <complex>
#include <lib.h>
#include <fstream>

using namespace arma;
using namespace std;

void HeliumAlphaZoom();
void HeliumAlpha();

int main()
{
    HeliumAlpha();
    HeliumAlphaZoom();

    /*
    myfile.open("Helium_atom_wf2.txt");
    for (int i=0; i<Nalpha; i+=4){
        for (int j=0; j<Nalpha; j+=4){
            VMCSolver helium = VMCSolver();
            helium.alpha = alpha(i);
            helium.beta = alpha(j);
            helium.WFnumber = 2;
            helium.FindStepLength();

            helium.MCintegration();
            myfile << alpha(i) << " " << alpha(j) << " " << helium.step << " " << helium.Energy << " " << helium.Variance << " " << helium.AcceptRate << endl;
        }
    }
    myfile.close();
    */

    /*
    VMCSolver helium = VMCSolver();
    helium.alpha = 1.0;
    helium.FindStepLength();
    helium.MCintegration();
    cout << helium.step << endl << helium.Energy;

    */
    /*
    VMCSolver helium = VMCSolver();
    helium.alpha = 1.0;
    helium.ImportanceSampling();
    cout << helium.Energy;
    */

    
    return 0;
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



