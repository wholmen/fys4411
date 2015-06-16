from pylab import *

def density(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9):
	x0 = sorted(x0); x1 = sorted(x1); x2 = sorted(x2); x3 = sorted(x3)
	x4 = sorted(x4); x5 = sorted(x5); x6 = sorted(x6)
	x7 = sorted(x7); x8 = sorted(x8); x9 = sorted(x9)
	h = 0.01

	N = -(x0[0] - x0[-1])/h

	hist = linspace(0, N*h, int(N)+3)
	psi = zeros(len(hist))

	for i in range(len(x0)):
		for n in range(len(hist)-1):
			if x0[i] > hist[n] and x0[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x1)):
		for n in range(len(hist)-1):
			if x1[i] > hist[n] and x1[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x2)):
		for n in range(len(hist)-1):
			if x2[i] > hist[n] and x2[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x3)):
		for n in range(len(hist)-1):
			if x3[i] > hist[n] and x3[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x4)):
		for n in range(len(hist)-1):
			if x4[i] > hist[n] and x4[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x5)):
		for n in range(len(hist)-1):
			if x5[i] > hist[n] and x5[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x6)):
		for n in range(len(hist)-1):
			if x6[i] > hist[n] and x6[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x7)):
		for n in range(len(hist)-1):
			if x7[i] > hist[n] and x7[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x8)):
		for n in range(len(hist)-1):
			if x8[i] > hist[n] and x8[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	for i in range(len(x9)):
		for n in range(len(hist)-1):
			if x9[i] > hist[n] and x9[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	psi = psi/len(x0)*2
	plot(hist,psi)
	xlabel('r')
	ylabel(r'$\rho(r)$')
	title('Probability density as a function of r with Jastrow factor')
	savefig('ProbabilityDensityNeon_Jastrow.png')
	show()

def readfile(filename,energy,x0,x1,x2=0,x3=0,x4=0,x5=0,x6=0,x7=0,x8=0,x9=0):

	infile0 = open(filename,'r')
	#n0,nall = infile0.readline().split();
	for line in infile0:
		splitline = line.split()
		energy.append(float(splitline[0]))
		x0.append(float(splitline[1]))
		x1.append(float(splitline[2]))
		if x2!= 0: x2.append(float(splitline[3]))
		if x3!= 0: x3.append(float(splitline[4]))
		if x2!= 0: x4.append(float(splitline[5]))
		if x3!= 0: x5.append(float(splitline[6]))
		if x2!= 0: x6.append(float(splitline[7]))
		if x3!= 0: x7.append(float(splitline[8]))
		if x2!= 0: x8.append(float(splitline[9]))
		if x3!= 0: x9.append(float(splitline[10]))
	if x2==0: return energy, x0,x1
	#elif x2!=0x4==0: return energy, x0,x1,x2,x3
	else: return energy, x0,x1,x2,x3,x4,x5,x6,x7,x8,x9


def blocking(e,e2):
	N = len(e)	# Total amount of samples
	
	n = linspace(10,10000,11); # Amount of samples in a block
	nb = N / n # Amount of blocks
	variance = zeros(len(n)) # Variance for this block size

	for i in range(len(n)): # Looping through the different block sizes
		localvar = zeros(nb[i])
		locale = zeros(nb[i]); locale2 = zeros(nb[i])
		for j in range( int(nb[i])): # Looping through all the blocks 
			
			E = 0; E2 = 0; elem=0;
			for k in range( j*int(n[i]), (j+1)*int(n[i]) ):
				E += e[k]; E2 += e[k]*e[k]
			
			Esum = E / n[i]; E2sum = E2 / n[i]
			localvar[j] = (E2sum - Esum**2)
			#locale[j] = Esum; locale2[j] = Esum**2
		
		variance[i] = sum(localvar) / nb[i]
		#variance[i] = 1.0/nb[i] * (1.0/nb[i]*sum(locale2) - (1.0/nb[i]*sum(locale))**2)

	plot(n,variance,'-')
	xlabel('number of samples in blocks')
	ylabel(r"$\sigma$")
	title('Error in the energy as a function of samples in block')
	#savefig('HeliumBlocking.png')
	show()



# Blocking and radial density Helium
"""
energy = []; r0 = []; r1 = []; 
energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_0.txt',energy,r0,r1)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_1.txt',energy,r0,r1)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_2.txt',energy,r0,r1)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_3.txt',energy,r0,r1)

r0 = array(r0); r1 = array(r1)#; x2 = array(x2); x3 = array(x3)
density(r0,r1)

energy = array(energy)
blocking(energy,energy*energy)
"""

# Blocking and radial density Beryllium
"""
energy = []; r0 = []; r1 = []; r2 = []; r3 = [] 
energy, r0,r1,r2,r3 = readfile('../build-project3_gaussian-Desktop-Debug/Beryllium_my_rank_0.txt',energy,r0,r1,r2,r3)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_1.txt',energy,r0,r1)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_2.txt',energy,r0,r1)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_3.txt',energy,r0,r1)

r0 = sqrt(array(r0)); r1 = sqrt(array(r1)); r2 = sqrt(array(r2)); r3 = sqrt(array(r3))
print "lol"
density(r0,r1,r2,r3)

energy = array(energy)
blocking(energy,energy*energy)
"""

# Blocking and radial density Neon
energy = []; r0 = []; r1 = []; r2 = []; r3 = []; r4 = []; r5 = []; r6 = []; r7 = []; r8 = []; r9 = []
energy, r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 = readfile('../build-project3_gaussian-Desktop-Debug/Neon_my_rank_0.txt',energy,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9)
energy, r0a,r1a,r2a,r3a,r4a,r5a,r6a,r7a,r8a,r9a = readfile('../build-project3_gaussian-Desktop-Debug/Neon_my_rank_1.txt',energy,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9)
r0 = r0+r0a; r1 = r1+r1a; r2 = r2+r2a; r3 = r3+r3a; r4 = r4+r4a; r5 = r5+r5a; r6 = r6+r6a; r7 = r7+r7a; r8 = r8+r8a; r9 = r9+r9a;
energy, r0a,r1a,r2a,r3a,r4a,r5a,r6a,r7a,r8a,r9a = readfile('../build-project3_gaussian-Desktop-Debug/Neon_my_rank_2.txt',energy,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9)
r0 = r0+r0a; r1 = r1+r1a; r2 = r2+r2a; r3 = r3+r3a; r4 = r4+r4a; r5 = r5+r5a; r6 = r6+r6a; r7 = r7+r7a; r8 = r8+r8a; r9 = r9+r9a;
energy, r0a,r1a,r2a,r3a,r4a,r5a,r6a,r7a,r8a,r9a = readfile('../build-project3_gaussian-Desktop-Debug/Neon_my_rank_3.txt',energy,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9)
r0 = r0+r0a; r1 = r1+r1a; r2 = r2+r2a; r3 = r3+r3a; r4 = r4+r4a; r5 = r5+r5a; r6 = r6+r6a; r7 = r7+r7a; r8 = r8+r8a; r9 = r9+r9a;
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_2.txt',energy,r0,r1)
#energy, r0,r1 = readfile('../build-project3_gaussian-Desktop-Debug/Helium_my_rank_3.txt',energy,r0,r1)

r0 = sqrt(array(r0)); r1 = sqrt(array(r1)); r2 = sqrt(array(r2)); r3 = sqrt(array(r3)); r4 = sqrt(array(r4)); r5 = sqrt(array(r5)); r6 = sqrt(array(r6))
r7 = sqrt(array(r7)); r8 = sqrt(array(r8)); r9 = sqrt(array(r9))
print "lol"
density(r0,r1,r2,r3,r4,r5,r6,r7,r8,r9)

energy = array(energy)

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