from sympy import *

r, alpha = symbols('r alpha')
f = (1-alpha*r/2.0)*exp(-alpha*r/2.0)

d1f = diff(f,r)
d2f = 1 / r**2 * diff( r**2 * diff(f,r), r)

f2 = alpha*r*exp(-alpha*r/2.0)
d1f2 = diff(f2,r)
d2f2 = 1 / r**2 * diff( r**2 * diff(f2,r), r)

print printing.ccode(d1f.simplify())
print printing.ccode(d2f.simplify())
print printing.ccode(d1f2.simplify())
print printing.ccode(d2f2.simplify())
