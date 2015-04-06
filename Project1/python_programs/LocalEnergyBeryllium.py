from sympy import *

alpha, beta, r1, r2, r3, r4, r12, r13, r14, r23, r24, r34, Z = symbols('alpha beta r1 r2 r3 r4 r12 r13 r14 r23 r24 r23 Z')

function = (exp(-alpha*r1) * (1-alpha*r2/2)*exp(-alpha*r2/2) - exp(-alpha*r2) * (1-alpha*r1/2)*exp(-alpha*r1/2) ) \
* (exp(-alpha*r3) * (1-alpha*r4/2)*exp(-alpha*r4/2) - exp(-alpha*r4) * (1-alpha*r3/2)*exp(-alpha*r3/2) ) \
* exp(r12 / 2*(1+beta*r12)) * exp(r13 / 2*(1+beta*r13)) * exp(r14 / 2*(1+beta*r14)) * exp(r23 / 2*(1+beta*r23)) * exp(r24 / 2*(1+beta*r24)) * exp(r34 / 2*(1+beta*r34))

LocalEnergy = -0.5/r1**2 * diff(r1**2 * diff(function,r1), r1) - 0.5/r1**2 * diff(r2**2 * diff(function,r2), r2) - 0.5/r3**2 * diff(r3**2 * diff(function,r3), r3) \
-0.5/r4**2 * diff(r4**2 * diff(function,r4), r4) - Z / r1 - Z / r2 - Z / r3 - Z / r4 + 1/r12 + 1/r13 + 1/r14 + 1/r23 + 1/r24 + r34

LocalEnergy = simplify(LocalEnergy)

print ccode(LocalEnergy)