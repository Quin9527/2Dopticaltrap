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


    def plot(y, Twodim,legend,marker):

        if Twodim:
            # plt.figure(1)
            # plt.plot(y['t'], 1e3 * y['pos_z'])
            # plt.xlabel('Time(s)')
            # plt.ylabel('Position_z (mm)')
            # plt.grid()
            # plt.figure(2)
            # plt.plot(y['t'], 1e3 * y['pos_r'])
            # plt.xlabel('Time (s)')
            # plt.ylabel('Position_r (μm)')
            # plt.grid()
            # plt.figure(3)
            # plt.plot(y['t'], 1e3 * y['vels_r'])
            # plt.xlabel('Time(s)')
            # plt.ylabel('Velocity_r (mm/s)')
            # plt.grid()
            # plt.figure(4)
            # plt.plot(y['t'], 1e3 * y['vels_z'])
            # plt.xlabel('Time (s)')
            # plt.ylabel('velocity_z (mm/s)')
            # plt.grid()
            # plt.figure(5)
            plt.plot(1e6 * y['pos_r'], 1e3 * y['pos_z'],label = legend)
            plt.xlabel('Position, a (μm)',fontsize = 15)
            plt.ylabel('Position, z (mm)',fontsize = 15)
            # plt.grid()
            # plt.show()
        else:

            i = len(y['t'])/20
            plt.plot(y['t'], 1e3 * y['pos_z'])#,label=legend,marker=marker,markersize=5,markevery=int(i)
            plt.xlabel('Time(s)',fontsize=14)
            plt.ylabel('Position (mm)',fontsize=14)
            plt.grid()
            # plt.figure(2)
            # plt.plot(y['t'], 1e3 * y['vels_z'])#,label=legend
            # plt.xlabel('Time (s)',fontsize=14)
            # plt.ylabel('velocity (mm/s)',fontsize=14)
            # plt.grid()
            # plt.show()


    os.chdir('E:\Imperial\Project\Simulation results\Frequency')
    # filename = str(input("file:"))
    i = 0
    # Twodim = False
    # #
    # #
    # filenumber = 1
    # y=[]

    # name = (['phm_1D.npz','Final_nophf_1D.npz','RK23_1D.npz'])
    # legend = (['low_1D.npz','No photophoretic force','RK23'])
    # marker = (['D', 's', 'x', 'o'])
    # while i < filenumber:
    #     # name = str(input("file name:"))
    #     tempfile = np.load(name[i],allow_pickle=True)
    #     y.append(tempfile)
    #     plot(y[i], Twodim,legend[i],marker[i])
    #     i+=1
    Twodim = True
    os.chdir('E:\Imperial\Project\Simulation results\TWODIM')
    filenumber = 2
    y=[]
    plt.figure(1)
    # Model v4 reader
    name = (['15tilt_hp_1_2D.npz','30tilt_hp_1_2D.npz','15_atm_2D.npz','B_Euler_τ_10s_1D.npz'])
    legend = (['5 degree','10 degree','100 μm','Euler τ'])
    marker = (['D','s','x','o'])
    while i < filenumber:
        # name = str(input("file name:"))
        tempfile = np.load(name[i],allow_pickle=True)
        y=tempfile
        plot(y, Twodim,legend[i],marker[i])
        i+=1
    plt.grid()
    plt.legend(loc='best')
    # 2D model reader
    plt.show()
    # tempfile = np.load(filename,allow_pickle=True)
    # plot(tempfile,Twodim=False)


