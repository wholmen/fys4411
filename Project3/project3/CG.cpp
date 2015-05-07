#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma;

double E_function(vec &x);

void dE_function(vec &x, vec &g);

void lnsearch(int n, vec &xold, double fold, vec &g, vec &p, vec &x, double *f, double stepmax, int *check, double (*func)(vec &p));

void dfpmin( vec &p, int n, double gtol, int *iter, double *fret, double (*func)(vec &p), void (*dfunc)(vec &p, vec &g));

#define SQUARE(a) ( a == 0.0 ? 0.0 : a*a ) // Return a^2 if a != 0

#define FMAX(a,b) ( a > b ? a : b) // Return the largest number

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STEPMAX 100.0

double E_function(vec &x) {return x(0)*x(0)*0.5 + 1.0/(8*x(0)*x(0));} // Expectation value of energy

void dE_function(vec &x, vec &g) {g(0) = x(0) - 1.0/(4*x(0)*x(0)*x(0));} // Derivative of <E_L> wrt alpha



void dfpmin(vec &p, int n, double gtol, int *iter, double *fret, double (*func)(vec &), void (*dfunc)(vec &, vec &)) {
    // Declaring different variables.
    int check, i, its, j;
    double den, fac, fad, fae, fp, stepmax, sum(0.0), sumdg, sumxi, temp, test;
    vec dg(n), g(n), hdg(n), pnew(n), xi(n);
    mat hessian(n,n);

    fp = (*func)(p); // Evaluate E_function for p and assigning result to fp
    (*dfunc)(p,g); // Evaluate dE_function for p and assigning result to g

    for (i=0; i<n; i++){
        for (j=0; j<n; j++) hessian(i,j)=0.0; // Zeroing out off-diagonal elements
        hessian(i,i)=1.0; // All diagonal elements in matrix A is one. A is a unit matrix
        xi(i) = -g(i); // Vector xi = -g
        sum += p(i)*p(i); // Summing up all elements in vector, p.
    }
    stepmax = STEPMAX * FMAX( sqrt(sum), (double)n );

    for (its=1; its<=ITMAX; its++){ // Starting main loop over all iterations in the iterative process
        *iter = its; // Just storing the amount of iterations used. iter is known outside dfpmin

        lnsearch(n,p,fp,g,xi,pnew,fret,stepmax,&check,func); // func is E_function. Function changes from xold to x
        fp = *fret; // Energy given by new function
        for (i=0; i<n; i++){
            xi(i) = pnew(i) - p(i); // xi is now the difference between old alpha and new alpha.
            p(i) = pnew(i);  // p represents alpha in this example. beta in project. Updating: pold = pnew.
        }
        test = 0.0;
        for (i=0; i<n; i++){
            temp = fabs(xi(i)) / FMAX( fabs(p(i)), 1.0); // (anew - aold) / anew . Or divided by 1 if anew < 1
            if (temp > test) test=temp; // Picks out the i giving the largest temp and saves it as test
        }
        if (test < TOLX) return; // Function is finished. Test was passed and the change in alpha is small enough.

        for (i=0; i<n; i++) dg(i) = g(i); // Storing old g
        (*dfunc)(p,g); // Calculating dE with new p. Stored in g
        test = 0.0;
        den = FMAX(*fret, 1.0);  // Den is the largest value of E(p_new) and 1.0
        for (i=0; i<n; i++){
            temp = fabs(g(i)) * FMAX(fabs(p(i)), 1.0) / den;
            if (temp > test) test = temp;
        }
        if (test < gtol) return; // Function is finished. The derivative is very small and test is passed.

        for (i=0; i<n; i++) dg(i) = g(i) - dg(i); // dg is now difference between new and old dE.
        for (i=0; i<n; i++) {
            hdg(i) = 0.0;
            for (j=0; j<n; j++) hdg(i) += hessian(i,j)*dg(j);
        }
        fac = fae = sumdg = sumxi = 0.0;
        for (i=0; i<n; i++){
            fac += dg(i)*xi(i); // delta dE times delta p.
            fae += dg(i)*hdg(i);
            sumdg += SQUARE(dg(i));
            sumxi += SQUARE(xi(i));
        }
        if (fac*fac > EPS){
            fac = 1.0/fac;
            fad = 1.0/fae;
            for (i=0; i<n; i++) dg(i) = fac*xi(i) - fad*hdg(i);
            for (i=0; i<n; i++) {
                for (j=0; j<n; j++){
                    hessian(i,j) += fac*xi(i)*xi(j) - fad*hdg(i)*hdg(j) + fae*dg(i)*dg(j);
                }
            }
        }
        for (i=0; i<n; i++) {
            xi(i) = 0.0;
            for (j=0; j<n; j++) xi(i) -= hessian(i,j)*g(j);
        }
    }
    cout << "too many iterations in dfpmin" << endl;
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STEPMAX

#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsearch(int n, vec &xold, double fold, vec &g, vec &p, vec &x, double *f, double stepmax, int *check, double (*func)(vec &)){
    int i;
    double a, alam(1.0), alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope(0.0), sum(0.0), temp, test(0.0), tmplam;

    *check = 0;
    for (i=0; i<n; i++) sum += p(i)*p(i); // Summing up all elements in vector p;
    sum = sqrt(sum);

    if (sum>stepmax) for (i=0; i<n; i++) p(i) *= stepmax/sum; //stepmax/sum = 100 if sum>n, else 100*n/sum

    for (i=0; i<n; i++) slope += g(i)*p(i);

    for (i=0; i<n; i++){
        temp = fabs(p(i)) / FMAX( fabs(xold(i)), 1.0);
        if (temp > test) test=temp; // Picks out the i giving the largest temp and saves it as test
    }

    alamin = TOLX / test;

    for (;;) { // Infinite loop
        for (i=0; i<n; i++) x(i) = xold(i) + alam*p(i);
        *f = (*func)(x); // Evaluating E_function (xnew). Storing it in *f. Aka fret in dfpmin.
        if (alam < alamin){
            for (i=0; i<n; i++) x(i) = xold(i);
            *check = 1;
            return; // Leaving function. No change done to x.
        }
        else if (*f <= fold+ALF*alam*slope) return; // Leaving function. Satisfied
        else{
            if (alam == 1.0) tmplam = -slope / (2.0*(*f-fold-slope));
            else{
                rhs1 = *f-fold-alam*slope;
                rhs2 = f2-fold2-alam2*slope;
                a = (rhs1/(alam*alam) - rhs2/(alam*alam2)) / (alam-alam2);
                b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2)) / (alam-alam2);
                if (a == 0) tmplam = -slope/(2.0*b);
                else {
                    disc = b*b - 3.0*a*slope;
                    if ( disc < 0.0 ) cout << "Roundoff problem in lnsearch." << endl;
                    else tmplam = (-b+sqrt(disc)) / (3.0*a);
                }
                if (tmplam > 0.5*alam) tmplam = 0.5*alam;
            }
        }
        alam2 = alam;
        f2 = *f;
        fold2 = fold;
        alam = FMAX(tmplam, 0.1*alam);
    }
}
#undef ALF
#undef TOLX



