from pylab import *

def density(x):
	x = sorted(x)
	h = 0.1
	N = -(x[0] - x[-1])/h
	Nhalf = int(N)/2+1
	hist = linspace(-Nhalf*h, Nhalf*h, int(N)+3)
	psi = zeros(len(hist))

	for i in range(len(x)):
		for n in range(len(hist)):
			if x[i] > hist[n] and x[i] < hist[n+1]:
				#print x[i], hist[n], hist[n+1]
				psi[n] += 1

	psi = psi/len(x)*2
	plot(hist,psi)
	xlabel('x')
	ylabel(r'$\rho(x)$')
	title('Probability density as a function of x')
	savefig('ProbabilityDensity.png')
	show()

def blocking():
	N = len(e)	# Total amount of samples
	n = array([float(1+i*100) for i in range(100)]) # different amount of blocks
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
	savefig('HeliumBlocking.png')
	show()

e = []
e2 = []
x = []
y = []
z = []

infile = open('../build-project1-Desktop-Debug/Blocking_file.txt','r')
for line in infile:
	spl = line.split()
	e.append(float(spl[0]))
	e2.append(float(spl[1]))
	x.append(float(spl[2]))
	#y.append(float(spl[3]))
	#z.append(float(spl[4]))

e = array(e); e2 = array(e2); x = array(x); y = array(y); z = array(z);

density(x)
blocking()
