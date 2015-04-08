from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing

alpha, Z, beta = symbols('alpha Z beta')
x1, y1, z1 = symbols('x1 y1 z1')
x2, y2, z2 = symbols('x2 y2 z3')

r1 = x1 + y1 + z1
r2 = x2 + y2 + z2

R1 = sqrt(x1*x1 + y1*y1 + z1*z1)
R2 = sqrt(x2*x2 + y2*y2 + z2*z2)

r12 = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) )

R12 = Symbol('r12'); r1s = Symbol('r1'); r2s = Symbol('r2'); R1s = Symbol('R1'); R2s = Symbol('R2')

phi = exp(-alpha*(R1 + R2))

ddr1 = (diff(diff(phi,x1),x1) + diff(diff(phi,y1),y1) + diff(diff(phi,z1),z1)).simplify().subs(R1,R1s).subs(R2,R2s)
ddr2 = (diff(diff(phi,x2),x2) + diff(diff(phi,y2),y2) + diff(diff(phi,z2),z2)).simplify().subs(R1,R1s).subs(R2,R2s)

EL1 = (0.5*ddr1 / phi + 0.5*ddr2 / phi - Z/R1 - Z/R2 + 1/r12).simplify().subs(R1,R1s).subs(R2,R2s).subs(r12,R12).collect(Z)

EL1s = Symbol('EL1')
phi = phi * exp( r12 / (2*(1+beta* r12)))

ddr1 = (diff(diff(phi,x1),x1) + diff(diff(phi,y1),y1) + diff(diff(phi,z1),z1))
ddr2 = (diff(diff(phi,x2),x2) + diff(diff(phi,y2),y2) + diff(diff(phi,z2),z2))

EL2 = (0.5*ddr1 / phi + 0.5*ddr2 / phi - Z/R1 - Z/R2 + 1/r12).subs(R1,R1s).subs(R2,R2s).subs(r12,R12).subs(EL1,EL1s).simplify().factor().subs(r1,r1s).subs(r2,r2s).simplify()

print EL2