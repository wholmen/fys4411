from pylab import *
import os

alpha = []; step = []; energy = []; variance = []
infile = open("../build-project1-Desktop-Debug/Helium_atom.txt",'r')
for line in infile:
	splitline = line.split()
	alpha.append(float(splitline[0]))
	step.append(float(splitline[1]))
	energy.append(float(splitline[2]))
	variance.append(float(splitline[3]))

plot(array(alpha),-array(energy),label="energy")
plot(alpha,variance,'-',label="variance")
legend(loc="lower left")
title('Absolute value of Energy and variance as a function of alpha')
xlabel(r"$\alpha$", fontsize=20)
grid('on')
savefig('EnergyVariance_helium1.png')
show()


alpha = []; step = []; energy = []; variance = []
infile = open("../build-project1-Desktop-Debug/Helium_atom_zoom.txt",'r')
for line in infile:
	splitline = line.split()
	alpha.append(float(splitline[0]))
	step.append(float(splitline[1]))
	energy.append(float(splitline[2]))
	variance.append(float(splitline[3]))

plot(array(alpha),array(energy),label="energy")
legend(loc="lower left")
title('Energy as a function of alpha')
xlabel(r"$\alpha$", fontsize=20)
grid('on')
savefig('EnergyVariance_helium2.png')
show()

plot(alpha,variance,'-',label="variance")

legend(loc="lower left")
title('Variance as a function of alpha')
xlabel(r"$\alpha$", fontsize=20)
grid('on')
savefig('EnergyVariance_helium3.png')
show()








alpha = []; beta = []; step = []; energy = []; variance = []
infile = open("../build-project1-Desktop-Debug/Helium_atom_wf2.txt",'r')
for line in infile:
	splitline = line.split()
	alpha.append(float(splitline[0]))
	beta.append(float(splitline[1]))
	step.append(float(splitline[2]))
	energy.append(float(splitline[3]))
	variance.append(float(splitline[4]))

plot(alpha,energy)
hold('on')
plot(beta,energy,'r')
show()