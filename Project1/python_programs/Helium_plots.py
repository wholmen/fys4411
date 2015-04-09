from pylab import *
import os

def Helium_FindAlpha():
	# Energy vs alpha for large alpha. 
	alpha = []; step = []; energy = []; variance = []
	infile = open("../build-project1-Desktop-Debug/Helium_atom.txt",'r')
	for line in infile:
		splitline = line.split()
		alpha.append(float(splitline[0]))
		step.append(float(splitline[1]))
		energy.append(float(splitline[2]))
		variance.append(float(splitline[3]))

	plot(array(alpha),array(energy),label="energy")
	plot(alpha,variance,'-',label="variance")
	legend(loc="lower left")
	title('Energy and variance as a function of alpha')
	xlabel(r"$\alpha$", fontsize=20)
	grid('on')
	ylim([-3,3])
	savefig('EnergyVariance_helium1.png')
	show()

def Helium_FindAlphaZoom():
	# Zoom in on small region of alpha
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


def Helium_FindBeta():
	alpha = []; beta = []; step = []; energy = []; variance = []
	infile = open("../build-project1-Desktop-Debug/Helium_atom_wf2.txt",'r')
	for line in infile:
		splitline = line.split()
		alpha.append(float(splitline[0]))
		beta.append(float(splitline[1]))
		step.append(float(splitline[2]))
		energy.append(float(splitline[3]))
		variance.append(float(splitline[4]))


	#plot(alpha,energy,label="Energy")
	hold('on')
	plot(alpha, energy,'o',label="Energy")
	title("Energy as a function of alpha and beta. Different beta's are plotted vertically.")
	xlabel(r"$\alpha$", fontsize=20)
	grid('on')
	legend()
	savefig("EnergyVariance_helium5")
	show()

	i = argmin(energy)
	print min(energy), alpha[i], beta[i]

def Helium_FindDT():
	dt = []; energy = []
	infile = open('../build-project1-Desktop-Debug/Compare_dt_Importance.txt','r')
	for line in infile:
		splitline = line.split()
		dt.append(float(splitline[0]))
		energy.append(float(splitline[1]))

	dt = array(dt); energy = array(energy)
	plot(dt, energy)
	xscale('log', nonposy='clip')
	title('Energy as dt varies using Importance sampling')
	xlabel('dt')
	ylabel('energy')
	ylim([-4,-2.5])
	savefig('ImportanceSampling_Helium_dt.png')
	show()


#Helium_FindAlpha()
#Helium_FindAlphaZoom()
#Helium_FindBeta()