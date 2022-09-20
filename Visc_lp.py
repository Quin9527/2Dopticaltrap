import numpy as np
import math
import matplotlib.pylab as plt
import Visc_lp

g = 9.8
Boltzmann = 1.38064852e-23
M = 28.9 * 1.6606e-27 #moleculare mass of air
miu = 1.827e-5
T = 300

class vis_(object):

    def __init__(self,radius,pressure):

        #self.p = presure
        self.r = radius
        #self.temp = T
        self.p= pressure
        #self.v = velocity

    def vis_(self):

        l = (miu/self.p)*np.sqrt(math.pi*Boltzmann*T/(2*M))
        vis_ = 6*math.pi*miu*self.r*(1+1.63*l/self.r)**(-1)

        return vis_


if __name__ == '__main__':
    x = np.linspace(10e-3,101300,10000)
    y = vis_(5e-6,x).vis_()
    plt.plot(1e9*y)
    plt.xlabel('Pressure (Pa)',fontsize = 15)
    plt.ylabel('Viscous drag (nNs/m)  ',fontsize = 15)
    plt.show()
