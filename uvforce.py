import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import random, copy, sys

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)



alpha =1.0 # self inhibition rate
beta = 1.0 # inhibition decay
gamma = 1.0 # bleaching receptivity
delta = .05 # self resensitization
r=.005
Nmax = 10
iters = 20000.0
dt = .01
duty =.5
maxintensity = 1.0
R_0 = 1.0

# alpha = float(sys.argv[1])
# beta = float(sys.argv[2])
# r = float(sys.argv[3])
# iters = int(sys.argv[4])
# maxintensity = float(sys.argv[5])
# duty = float(sys.argv[6])


length = iters
period = length/5.0

Ms = np.zeros(iters)
Mdots =np.zeros_like(Ms)
Resists = np.zeros_like(Ms)
Resistdots = np.zeros_like(Ms)

Resists[0] = 0
#Ms[0] = alpha*maxintensity/beta*.01
force = np.zeros_like(Ms)


# construct intensity profile. duty in percentage, length must be divisible by period
def heaviside(x):
    return 1 * (x>=0)

def construct_intensity(duty, maximum, period, length):
    intensity = np.zeros(int(length))
    for i in range(0, int(length/period)):
        framecount = 0
        while framecount < period:
            frame = i*int(period)+framecount
            if framecount/period < duty:
                intensity[frame] = maximum
            else:
                intensity[frame] = 0
                
            framecount+=1
    return intensity

intensity = construct_intensity(duty, maxintensity, period, length)

#main simulation    


for i in range(1,len(Ms)):
    # IFFL case
    dResist = (alpha*intensity[i-1] - beta*Resists[i-1])*dt
    # NFL case
    #dResist = (intensity[i-1] - beta*Resists[i-1])*dt
    Resists[i] = Resists[i-1] + dResist
    dM = (gamma*intensity[i-1]*heaviside(R_0 - Resists[i-1]) - delta*Ms[i-1])*dt
    Ms[i] = Ms[i-1] + dM
    Mdots[i] = dM/dt 
    Resistdots[i] = dResist/dt
    #dM = r*(-beta*Ms[i-1]**2 + alpha*intensity[i-1]*Ms[i-1])*dt
    #Mdots[i] = dM/dt
    #Ms[i] = Ms[i-1] + dM
    # if Ns[i]<0:
    #     Ns[i]=0
    #force[i] = max(gamma*Mdots[i],0)

for j in range(0,len(Mdots)):
    force[j] = max([Mdots[j],0])

#print "Ns is: ", Ns, "\n\n", "Force is: ", force
plt.figure('Ms')
plt.plot(np.arange(0,iters*dt,dt), Mdots)
plt.title('$F$ vs $t$ for %i iterations ' % (iters,))
plt.xlabel('$t$')
plt.ylabel('$M$')
plt.tight_layout(pad=0.2)
matplotlib.pyplot.show()

# plt.savefig("ms-iters%i-duty%1.1f-alpha%02.2f-beta%02.2f-r%02.3f.pdf" % (iters,duty,alpha, beta, r))

# plt.figure('Mdots')
# plt.plot(np.arange(0,iters*dt,dt),Mdots)
# plt.title('$\dot M$ vs $t$ for %i iterations with $\\alpha=%g$, $\\beta=%g$, $r=%g$, duty$=%g$' % (iters, alpha, beta, r,duty))
# plt.xlabel('$t$')
# plt.ylabel('$\dot M$')
# plt.tight_layout(pad=0.2)
# plt.savefig("mdots-iters%i-duty%1.1f-alpha%02.2f-beta%02.2f-r%02.3f.pdf" % (iters,duty,alpha, beta, r))

