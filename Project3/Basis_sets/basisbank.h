//This file is maintained by an external python script and should not be edited manually.
#ifndef BASISBANK_H
#define BASISBANK_H
#include <armadillo>
#include <basis.h>
#include <primitive.h>
using namespace std;
using namespace arma;
 
class basisbank{
public:
    basisbank(basis BS);
    basisbank();
    basis bs;
    string basistype;};
#endif // BASISBANK_H