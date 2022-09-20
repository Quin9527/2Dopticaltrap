import matplotlib.pyplot
import scipy.integrate
import sympy
import math
import numpy as np
import pylab as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


m0 = 4 * math.pi * 1e-7  # magnetic permeability of vacuum [H/m]
e0 = 8.8541878176e-12  # electric permittivity of vacuum [F/m]
c0 = math.sqrt(1 / (m0 * e0))  # speed of light in vacuum [m/s]
ni = 1
nt = 1.405  # silica particle refractive index

class optf(object):

    def __init__(self,beamwaist,beampower,dist_r,dist_z,radius,rayleigh,fs_switch,fg_switch,onedim):

        self.bw = beamwaist
        self.pi = beampower
        self.z = dist_z
        self.x = dist_r
        self.rayleigh = rayleigh
        self.r = radius
        self.fs_bool = fs_switch
        self.fg_bool = fg_switch
        self.onedim = onedim
        if self.fs_bool==False and self.fg_bool==False:
            return 0

    def theta_2(self,theta_1):

        theta_2 = math.asin((ni * math.sin(theta_1)) / nt)
        return theta_2

    def R(self,theta_1):

        theta_2 = self.theta_2(theta_1)
        R = 1/2*((np.sin(theta_1-theta_2)/np.sin(theta_1+theta_2))**2+
                 (np.tan(theta_1-theta_2)/np.tan(theta_1+theta_2))**2)
        return R

    def Intens_(self,ra):

        waist_z = self.cal_waist(self.z)
        P0 = 2 * self.pi / (math.pi * waist_z ** 2)
        I = P0* np.exp(-2 * ra ** 2 / waist_z**2)

        return I

    def cal_waist(self,z):

        waist_z = self.bw * math.sqrt(1+(z/self.rayleigh)**2)
        return waist_z

    def rho(self,theta_1,phi):

        rho = np.sqrt(self.x**2+self.r**2*np.sin(theta_1)**2-2*self.x*self.r*
                      np.sin(theta_1)*np.cos(phi))
        return rho

    def fs_func(self,theta_1,phi):

        theta_2 = self.theta_2(theta_1)
        R = self.R(theta_1)
        rho = self.rho(theta_1,phi)
        T = 1-R
        return ni/c0*self.Intens_(rho)*(1+
                   R*np.cos(2*theta_1)-T**2*(np.cos(2*theta_1-
                             2*theta_2)+R*np.cos(2*theta_1))/(1+R**2+2*R*np.cos(2*theta_2)))*self.r**2*\
                                    np.sin(2*theta_1)

    def fg_func(self,theta_1,phi):

        theta_2 = self.theta_2(theta_1)
        R = self.R(theta_1)
        rho = self.rho(theta_1,phi)
        T = 1-R
        return -ni/c0*self.Intens_(rho)*(
                   R*np.sin(2*theta_1)-T**2*(np.sin(2*theta_1-
                             2*theta_2)+R*np.sin(2*theta_1))/(1+R**2+2*R*np.cos(2*theta_2)))*self.r**2*\
                                    np.sin(2*theta_1)*np.cos(phi)

    def cal_f(self):

        Fs = 0
        Fg = 0
        if self.fs_bool:
            fs = self.fs_func
            # if self.onedim:
            #     theta_step = 60
            #     phi_step = 1
            #     dtheta_1 = math.pi/(2*theta_step)
            #     dphi = 2*math.pi/phi_step
            #     for i in np.arange(1e-3,math.pi/2,math.pi/(2*theta_step)):
            #             Fs += fs(i,2*math.pi)*dphi*dtheta_1
            #     return Fs
            # else:
            Fs = scipy.integrate.dblquad(fs,0,2*math.pi,
                                         lambda theta_1 : 0,lambda theta_1 : math.pi/2)

            return Fs[0]

        if self.fg_bool:

            fg = self.fs_func
            # dtheta_1 = math.pi/(2*theta_step)
            # dphi = 2*math.pi/phi_step
            # for i in np.arange(-math.pi/2+1e-3,math.pi/2,math.pi/(2*theta_step)):
            #     for j in np.arange(1e-3,2*math.pi,2*math.pi/phi_step):
            #         Fg += fg(i,j)*dphi*dtheta_1

            fg = self.fg_func
            Fg = scipy.integrate.dblquad(fg,0,2*math.pi,
                                         lambda theta_1 : 0,lambda theta_1 : math.pi/2)
            return Fg[0]


if __name__ == '__main__':
    z = np.linspace(-500e-6,500e-6,1000)
    wavelength = 532e-9
    M2=1.168
    rayleigh = math.pi * 2e-6 ** 2 / (wavelength*M2)

    # plt.figure(1)
    # plt.title('Optical force along z-direction')
    i = 0
    x = np.linspace(-10e-6,10e-6,200)
    #
    # while i <= 3:
    #     Fs_z = []
    #     for index in z:
    #
    #         Fs_z.append(optf(2e-6,0.4,0,index, 1e-6*2**i, rayleigh,True,False,False).cal_f())
    #     ndarray = np.array(Fs_z)
    #     char = str(2**i)
    #     plt.plot(1000*z,1e9*ndarray,label='Radius:'+char+'μm')
    #     i+=1
    #
    #
    # plt.legend(loc='best')
    # plt.xlabel('Distance from focal spot (mm)')
    # plt.ylabel('Optical scatter force (nN)')

    # plt.figure(2)
    # plt.title('Optical force along transverse direction')
    #
    # i = 0
    # while i <= 3:
    #     Fg_r = []
    #     for index in x:
    #         Fg_r .append(-optf(2e-6,0.4,index,0, 1e-6*2**i, rayleigh,False,True,False).cal_f())
    #     char = str(2**i)
    #     ndarray = np.array(Fg_r)
    #     plt.plot(1e6*x,1e9*ndarray,label='Radius:'+char+'μm')
    #     i+=1
    #
    # plt.legend(loc='best')
    # plt.xlabel('Distance from beam center axis (μm)')
    # plt.ylabel('Optical gradient force (nN)')
    # #

    # plt.figure(1)
    # plt.title('Optical force along z-direction\nParticle radius: 1.5 μm')
    # i = 0
    # x = np.linspace(-10e-6,10e-6,200)
    #
    # while i <= 1:
    #     Fs_z = []
    #     for index in z:
    #
    #         Fs_z.append(optf(1.9e-6+0.1e-6*i,0.4,0,index, 5e-6, rayleigh,True,False,False).cal_f())
    #     ndarray = np.array(Fs_z)
    #     char = str(2**i-1)
    #     plt.plot(1000*z,1e9*ndarray,label='Offset:'+char+'μm')
    #     i+=1
    #
    #
    # plt.legend(loc='best')
    # plt.xlabel('Distance from focal spot (mm)')
    # plt.ylabel('Optical scatter force (nN)')
    # plt.show()

    # plt.figure(1)
    # plt.title('Optical force as a function of offset\nParticle radius: 1.5 μm')
    # i = 0
    # x = np.linspace(-10e-6,10e-6,200)
    #
    # while i <= 3:
    #     Fs_z = []
    #     for index in x:
    #
    #         Fs_z.append(-optf(2e-6,0.4,index,(2**i-1)*10*1e-6, 1.5e-6, rayleigh,False,True,False).cal_f())
    #     ndarray = np.array(Fs_z)
    #     char = str(2**i-1)
    #     plt.plot(1e6*x,1e9*ndarray,label='z-direction:'+char+'μm')
    #     i+=1
    #
    #
    # plt.legend(loc='best')
    # plt.xlabel('Distance from beam center axis (μm)')
    # plt.ylabel('Optical gradient force (nN)')


    #
    fig = plt.figure()
    ax3 = plt.axes(projection='3d')
    # plt.title('Optical force in different position\nParticle radius: 1.5 μm')
    i = 0
    x = np.linspace(-1e-6, 1e-6, 40)
    z = np.linspace(-1e-6, 1e-6, 40)
    X,Z = np.meshgrid(x,z)
    Fs_z = np.zeros([40,40])
    for index in x:
        for index_ in z:
            X_index = np.where(x==index)
            Z_index = np.where(z==index_)
            Fs_z[Z_index,X_index]=-optf(2e-6, 0.4, 4*index, 100*index_,
                                       1.5e-6, rayleigh, False, True, False).cal_f()

    ax3.plot_surface(1000*Z*1e3,4*X*1e6,1e9*Fs_z,rstride = 1, cstride=1,cmap='rainbow')
    plt.gca().xaxis.set_major_locator(matplotlib.pyplot.MultipleLocator(0.5))
    ax3.set_xlabel('z, axial distance (mm)',fontsize = 12)
    ax3.set_ylabel('a, transverse offset (μm)',fontsize = 12)
    ax3.set_zlabel('Optical gradient force (nN)',fontsize = 12)


    plt.show()


