import numpy as np
import math
import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
g = 9.8
Boltzmann = 1.38064852e-23
M = 28.9 * 1.6606e-27 #moleculare mass of air
T = 300

class buo(object):

    def __init__(self,radius,pressure):

        #self.p = presure
        self.r = radius
        #self.temp = T
        self.den = M*pressure/(Boltzmann*T)

    def buo(self):

        buo = self.den*(4/3*math.pi*self.r**3)*g
        return buo

if __name__ == '__main__':

    x = np.linspace(0,101300,1000)
    y = buo(5e-6,x).buo()
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    # plt.axes(xscale='log')
    plt.plot(x,1e15*y)
    plt.xlabel('Pressure (Pa)',fontsize=15)
    plt.ylabel('buoyancy (fN)',fontsize=15)
    plt.show()






