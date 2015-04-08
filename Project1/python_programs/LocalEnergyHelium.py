from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing

alpha, Z, beta = symbols('alpha Z beta')
x1, y1, z1 = symbols('x1 y1 z1')
x2, y2, z2 = symbols('x2 y2 z2')

r1 = x1 + y1 + z1
r2 = x2 + y2 + z2

R1 = sqrt(x1*x1 + y1*y1 + z1*z1)
R2 = sqrt(x2*x2 + y2*y2 + z2*z2)

r12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) )
r21 = sqrt( (-x1+x2)**2 + (-y1+y2)**2 + (-z1+z2)**2 )
r1r2 = x1*x2 + y1*y2 + z1*z2

r12s = Symbol('r12'); r1s = Symbol('r1'); r2s = Symbol('r2'); R1s = Symbol('R1'); R2s = Symbol('R2'); r1r2s = Symbol('r1r2'); r21s = Symbol('r21')

phi = exp(-alpha*(R1 + R2))

ddr1 = (diff(diff(phi,x1),x1) + diff(diff(phi,y1),y1) + diff(diff(phi,z1),z1)).simplify().subs(R1,R1s).subs(R2,R2s)
ddr2 = (diff(diff(phi,x2),x2) + diff(diff(phi,y2),y2) + diff(diff(phi,z2),z2)).simplify().subs(R1,R1s).subs(R2,R2s)

EL1 = (0.5*ddr1 / phi + 0.5*ddr2 / phi - Z/R1 - Z/R2 + 1/r12).simplify().subs(R1,R1s).subs(R2,R2s).subs(r12,r12s).collect(Z)

EL1s = Symbol('EL1')
phi = phi * exp( r12 / (2*(1+beta* r12)))

ddr12 = (diff(r12,x1) + diff(r12,y1) + diff(r12,z1)).subs(r12,r12s).factor().subs(r1,r1s).subs(r2,r2s)
ddr21 = (diff(r12,x2) + diff(r12,y2) + diff(r12,z2)).subs(r12,r12s).factor().subs(r1,r1s).subs(r2,r2s)



#EL2 = -0.5*(-alpha**2 + 2*alpha*ddr1JF - d2dr1JF) - 0.5*(-alpha**2 + 2*alpha*ddr2JF - d2dr2JF) - 2/R1 - 2/R2 + 1/r12



ddr1 = (diff(diff(phi,x1),x1) + diff(diff(phi,y1),y1) + diff(diff(phi,z1),z1)).subs(R1,R1s).subs(R2,R2s).subs(r12,r12s).simplify()

ddr2 = (diff(diff(phi,x2),x2) + diff(diff(phi,y2),y2) + diff(diff(phi,z2),z2)).subs(R1,R1s).subs(R2,R2s).subs(r12,r12s).simplify()

T = (ddr1 + ddr2).simplify().subs(R1,R1s).subs(R2,R2s).subs(r12,r12s).subs(R1**2,R1s**2).subs(R2**2,R2s**2).subs(r12**2,r12s**2).subs(r21,r12s).subs(r21**2,r12s**2)

print T.simplify()


#EL2 = (0.5*ddr1 / phi + 0.5*ddr2 / phi - Z/R1 - Z/R2 + 1/r12).subs(R1,R1s).subs(R2,R2s).subs(r12,R12).subs(r1r2,r1r2s).simplify().subs(EL1,EL1s)

#print EL2





"""
from sympy import *

alpha, beta = symbols('alpha beta')
r1 = symbols('r1', real=True)
r2 = symbols('r2', real=True)
r12 = Abs(r1 - r2)
R12 = Symbol('r12')


function1 = exp(-alpha * (r1 + r2))
function2 = exp(-alpha * (r1 + r2)) * exp( r12 / (2*(1+beta* r12)))

print diff(r12, r1)


dx = diff(r12,r1)
print dx

#print dx.subs(r12,R12)


def calculate(function):
	LocalEnergy = -0.5 /r1**2 * diff(r1**2 * diff(function,r1), r1) / function - 0.5/r2**2 * diff(r2**2 * diff(function,r2), r2) / function - 2 / r1 - 2 / r2 + 1/r12
	#LocalEnergy = expand(LocalEnergy)
	#LocalEnergy = factor(LocalEnergy)
	LocalEnergy = simplify(LocalEnergy)
	return LocalEnergy


LocalEnergy = calculate(function2)

print printing.ccode(LocalEnergy)



#double sign(double a) {
#	return (a>0) ? 1 : -1;
#}

"""