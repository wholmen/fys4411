from pylab import * 
from sympy import *

x, y, z, i, j, k, a, c, N = symbols('x y z i j k a c N')
r = sqrt(x**2 + y**2 + z**2)
r2 = x**2 + y**2 + z**2
rs = Symbol('r'); rs2 = Symbol('r2')

xi = c*N*x**i*y**j*z**k*exp(-a * r**2)


print printing.ccode(diff(xi,z).simplify().subs(r2,rs2))

dxi = (diff(xi,x) + diff(xi,y) + diff(xi,z) ).subs(r2,rs2).simplify().subs(r,rs);
d2xi = (diff(diff(xi,x),x) + diff(diff(xi,y),y) + diff(diff(xi,z),z)).subs(r2,rs2).simplify()

print printing.ccode(d2xi)
print printing.ccode(dxi)