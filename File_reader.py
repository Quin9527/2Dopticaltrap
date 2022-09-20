import numpy as np
import matplotlib.pylab as plt
import scipy.integrate as inte
import math
import sympy
import os
from textwrap import fill

g = 9.8
mass = 970*4/3*math.pi*(5e-6)**3
if __name__ == '__main__':

    def eng(y):

        Energy = []
        pos = y['pos']
        vels = y['vels']
        for i in range(len(y['pos'])):
            Eg = (1 / 2) * mass * vels[i] ** 2
            Energy.append(Eg)
        return Energy

    def eng_2d(y):

        Energy = []
        pos = y['pos_z']
        vels = y['vels_z']
        for i in range(len(y['pos_z'])):
            Eg = (1 / 2) * mass * vels[i] ** 2
            Energy.append(Eg)
        return Energy


    def plot(y,filenumber,legend):

        plt.figure(1)
        i = 0
        j = 0
        k = 0
        plt.xlabel('Time(s)')
        plt.ylabel('Position (mm)')
        while i < filenumber:
            plt.plot(y[i]['t'], 1e3*y[i]['pos'],label=legend[i])
            i += 1
        plt.legend(loc='best')
        plt.grid()

        plt.figure(2)
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (mm/s)')
        while j < filenumber:
            plt.plot(y[j]['t'], 1e3*y[j]['vels'],label=legend[j])
            j += 1
        plt.legend(loc='best')
        plt.grid()

        plt.figure(3)
        plt.axes(yscale="log")
        plt.plot(y[k]['pres'][0,:],y[k]['pres'][1,:])
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (Pa)')

        plt.figure(4)
        plt.xlabel('Time(s)')
        plt.ylabel('Kinetic Energy (J)')
        while k < filenumber:
            Eng = eng(y[k])
            plt.plot(y[k]['t'], Eng, label=legend[k])
            k += 1
        plt.legend(loc='best')

        plt.grid()
        plt.show()

    def plot_2d1d(y,filenumber,legend):

        plt.figure(1)

        plt.xlabel('Time(s)')
        plt.ylabel('Position (mm)')

        plt.plot(y['t'], 1e3*y['pos_z'],label=legend)
        plt.legend(loc='best')
        plt.grid()

        plt.figure(2)
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (mm/s)')
        plt.plot(y['t'], 1e3*y['vels_z'],label=legend)
        plt.legend(loc='best')
        plt.grid()

        plt.figure(3)
        plt.xlabel('Time(s)')
        plt.ylabel('Kinetic Energy (J)')
        Eng = eng_2d(y)
        plt.plot(y['t'], Eng, label=legend)
        plt.legend(loc='best')

        plt.grid()
        plt.show()

    os.chdir('E:\Imperial\Project\Simulation results')
    filenumber = int(input("number of files:"))
    y = []
    i = 0
    # Model v4 reader
    # name = (['Allforce.npz','WithoutPhforce.npz','Withoutnoise.npz'])
    # legend = (['Allforce','No photophoretic force','No brownian motion'])
    # while i < filenumber:
    #    # name = str(input("file name:"))
    #     tempfile = np.load(name[i],allow_pickle=True)
    #     y.append(tempfile)
    #     i+=1
    #
    # plot(y,filenumber,legend)

    # 2D model reader
    name = (['2Dmodel.npz','2Dmodel_1dtest.npz'])
    legend = ("Force included:\nOptical force,\nBrownian motion,\nviscos drag,\nPhotophoretic force ")
    while i < filenumber:
       # name = str(input("file name:"))
        tempfile = np.load(name[i],allow_pickle=True)
        y.append(tempfile)
        i+=1

    plot_2d1d(y[1],filenumber,legend)


