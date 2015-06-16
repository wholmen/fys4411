from pylab import *

def density(x):
	x = sorted(x)
	h = 0.1
	N = -(x[0] - x[-1])/h
	Nhalf = int(N)/2+1
	hist = linspace(-Nhalf*h, Nhalf*h, int(N)+3)
	psi = zeros(len(hist))

	for i in range(len(x)):
		for n in range(len(hist)-1):
			if x[i] > hist[n] and x[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	psi = psi/len(x)*2
	plot(hist,psi)
	xlabel('x')
	ylabel(r'$\rho(x)$')
	title('Probability density as a function of x')
	savefig('ProbabilityDensityBeryllium.png')
	show()

def readfile(filename,energy,x0,x1,x2=0,x3=0):

	infile0 = open(filename,'r')
	n0,nall = infile0.readline().split();
	for line in infile0:
		splitline = line.split()
		energy.append(float(splitline[0]))
		x0.append(float(splitline[1]))
		x1.append(float(splitline[2]))
		if x2!= 0: x2.append(float(splitline[3]))
		if x3!= 0: x3.append(float(splitline[4]))
	return energy, x0,x1,x2,x3

def blocking(e,e2):
	N = len(e)	# Total amount of samples
	#n = array([float(1+i*100) for i in range(100)]) # different amount of blocks
	n = linspace(1,N/100,100)
	nb = zeros(len(n)) # Block sizes
	for i in range(len(n)):
		nb[i] = int(N/n[i]) # Size of blocks corresponding to amount of blocks. indexes of n and nb are connected.

	variance = zeros(len(n)) # Variance for this block size

	for i in range(len(n)): # Looping through every amount of blocks. More and more blocks for every i
		MeanEnergy = zeros(n[i])  # An array of energy mean for ever block
		MeanEnergy2 = zeros(n[i]) # An array of energy squared mean for every block

		for j in range(len(MeanEnergy)): # Looping through every index in energy mean to fill it up with values

			for k in range(int(j*nb[i]), int((j+1)*nb[i])): # Looping through all samples in block j. 
				MeanEnergy[j] += e[k]  # Adding the sample
				MeanEnergy2[j] += e2[k] # Adding the sample

		MeanEnergy = MeanEnergy/nb[i] # Dividing by amount of samples in each block
		MeanEnergy2 = MeanEnergy2/nb[i]

		# Now I have a list of means for e and e^2 if each block. 

		m = sum(MeanEnergy)/len(MeanEnergy); m2 = sum(MeanEnergy2)/len(MeanEnergy2)
		variance[i] = m2 - m*m

	plot(n,variance,'o')
	xlabel('number of blocks')
	ylabel(r"$\sigma$")
	title('Error in the energy as a function of the number of blocks')
	savefig('BerylliumBlocking.png')
	show()

"""
energy = []; x0 = []; x1 = []; x2 = []; x3 = []
energy, x0,x1,x2,x3 = readfile('../build-project3-Desktop-Debug/Beryllium_my_rank_0.txt',energy,x0,x1,x2,x3)
energy, x0,x1,x2,x3 = readfile('../build-project3-Desktop-Debug/Beryllium_my_rank_1.txt',energy,x0,x1,x2,x3)
energy, x0,x1,x2,x3 = readfile('../build-project3-Desktop-Debug/Beryllium_my_rank_2.txt',energy,x0,x1,x2,x3)
energy, x0,x1,x2,x3 = readfile('../build-project3-Desktop-Debug/Beryllium_my_rank_3.txt',energy,x0,x1,x2,x3)


x0 = array(x0); x1 = array(x1); x2 = array(x2); x3 = array(x3)
density(x0)
#density(x1)
#density(x2)
#density(x3)


energy = array(energy)
blocking(energy,energy*energy)
"""
energy = []; r0 = []; r1 = []; 
energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_0.txt',energy,r0,r1)
energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_1.txt',energy,r0,r1)
energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_2.txt',energy,r0,r1)
energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_3.txt',energy,r0,r1)

x0 = array(x0); x1 = array(x1); x2 = array(x2); x3 = array(x3)
density(x0)

energy = array(energy)
blocking(energy,energy*energy)

# Finding optimal alpha for Beryllium.
"""
infile = open('../build-project3-Desktop-Debug/Beryllium_CompareAlpha.txt','r')
E = []; alpha = []; Error = []
for line in infile:
	splitline = line.split()
	if len(splitline) != 0:
		if splitline[0] == 'Alpha:':
			alpha.append(float(splitline[1])); 
		if splitline[0] == 'Energy:':
			E.append(float(splitline[1]))
		if splitline[0] == 'Error:':
			Error.append(float(splitline[1]))

alpha = array(alpha); E = array(E)
plot(alpha,E)
title('Energy of Beryllium plotted against alpha')
xlabel('alpha')
ylabel('Energy')
savefig('FindOptimalAlphaBeryllium.png')
show()
"""

"""
# FInding optimal alpha for Neon
infile = open('../build-project3-Desktop-Debug/Neon_findalpha.txt','r')
E = []; alpha = []; Error = []
for line in infile:
	splitline = line.split()
	if len(splitline) != 0:
		if splitline[0] == 'Alpha:':
			alpha.append(float(splitline[1])); 
		if splitline[0] == 'Energy:':
			E.append(float(splitline[1]))
		if splitline[0] == 'Error:':
			Error.append(float(splitline[1]))

alpha = array(alpha); E = array(E)
plot(alpha,E)
title('Energy of Neon plotted against alpha')
xlabel('alpha')
ylabel('Energy')
savefig('FindOptimalAlphaNeon.png')
show()
"""