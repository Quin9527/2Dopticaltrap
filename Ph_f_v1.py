import numpy as np
import math
import matplotlib.pylab as plt

miu = 1.827e-5 #dynamic viscosity of air
kg = 0.025 #thermal conductivity of air in slip flow regime
T = 300  #ambient temperature
ks = 0.180  #silicone oil thermal conductivity
lam = 532e-9
cs = 1.17
ct = 2.18
cm = 1.14
#J1 = 0.0002 #asymetic parameter negative photophoresis
M = 28.9 * 1.6606e-27 #moleculare mass of air
Boltzmann = 1.38064852e-23
R = 0.287 #gas constant of air
alpha = 0.01 #absorption
n = 1.405 #refractive index of silicone oil
#M = 28.9 * 1.6606e-27

class ph_f(object):

    def __init__(self,power,waist,pressure,z_distance,radius,r_distance):

        self.p = power
        self.bw = waist
        self.press = pressure
        self.z = z_distance
        self.r = radius
        self.x = r_distance
        step = 30
        # self.l = (miu / self.press) * np.sqrt(math.pi * Boltzmann * T / (2 * M))  # mean free path

    def cal_J1(self):

        r = 2 * math.pi*self.r/lam
        k = (4*math.pi**-1)*alpha*lam*10**2
        J1 = 2*n*k*r*(3*(n-1)/8*n**2-(2/5)*n*r*k)
        return J1

    def cal_waist(self):

        rayleigh = math.pi*self.bw**2/lam
        waist_z = self.bw *math.sqrt(1+(self.z/rayleigh)**2)
        return waist_z

    def I(self,r):

        I0 = 2*self.p/(math.pi*self.bw**2)
        I = I0*((self.bw/self.cal_waist())**2)*math.exp(-(2*r**2)/(self.cal_waist()**2))
        return I

    def fp_c(self,mfp,step):

        l = mfp
        kg_ = 15 * miu * R / 4
        dens = self.mfp2dens(mfp)
        J1 = self.cal_J1()
        Qp_c = -4*math.pi*cs*J1*(miu**2)/(dens*ks*T)*(
            (1+3*cm*l/self.r)*(1+2*ct*l/self.r+2*kg/ks))**-1
        Fp_c = 0
        for i in np.arange(-step,step,1):

            Fp_c += Qp_c*self.I(i*self.r/step-self.x)*self.r/step/2

        return Fp_c

    def fp_f(self,mfp,step):

        l = mfp
        kg_ = 15*miu*R/4
        Hm = 4*kg_*self.r/(15*l*ks)
        J1 = self.cal_J1()
        dens = self.mfp2dens(mfp)
        Qp_f = -((math.pi**2)*(miu**2)*J1/(6*dens*ks*T))*(
                ((self.r/l)**2)/(1+Hm))
        Fp_f =0
        for i in np.arange(-step,step,1):

            Fp_f += Qp_f*self.I(i*self.r/step-self.x)*self.r/step/2

        return Fp_f

    def mfp2dens(self,mfp):

        pressure = (miu / mfp) * math.sqrt(math.pi * Boltzmann * T / (2 * M))
        dens = M * pressure / (Boltzmann * T)
        return dens

    def pres2mfp(self,pres):

        mfp = (miu / pres) * np.sqrt(math.pi * Boltzmann * T / (2 * M))
        return mfp

    def total_fp(self,step):

        Kn = self.pres2mfp(self.press)/self.r

        # if Kn < 0.01:
        #     fp_c = self.fp_c(Kn * self.r, step)
        #     return fp_c
        # if Kn > 100000:
        #     fp_f = self.fp_f(Kn * self.r, step)
        #     return fp_f
        # if Kn>0.01 and Kn <= 100000:
        fp_c = self.fp_c(Kn * self.r, step)
        fp_f = self.fp_f(Kn * self.r, step)
        if fp_c + fp_f == 0:
            return 0
        else:
            tran_fp = (fp_f*fp_c/(fp_c+fp_f))
            return tran_fp

    def test_total_fp(self,step):

        Kn = self.pres2mfp(self.press)/self.r
        total_fp=[]
        for index in Kn:
            fp_f = self.fp_f(index * self.r, step)
            fp_c = self.fp_c(index * self.r, step)
            if index < 0.01:
                 total_fp.append(fp_c)
            if index > 100000:
                total_fp.append(fp_f)
            if index>0.01 and index <= 100000:
                total_fp.append(fp_f*fp_c/(fp_c+fp_f))
        return total_fp

if __name__ == '__main__':

    x1 = np.logspace(-5,-3,1000,base=10)
    x2 = np.logspace(-3,5,1000,base=10)
    x = np.append(x1,x2)
    y = []
    # mfp = ph_f(0.4,2.65e-6,x,0,5e-6).pres2mfp(x)
    plt.figure(1)
    plt.axes(xscale="log")
    for i in range(1,4):
     step = 50
     y = np.array((ph_f(0.4, 1.87e-6, x, (i-1)*1e-5, 5e-6, 0).test_total_fp(step)))*1e9

    # y2 = ph_f(0.4,2.65e-6,x,0,5e-6).fp_c(mfp)
    # y3 = ph_f(0.4,2.65e-6,x,0,5e-6).fp_f(mfp)

     plt.plot(x,y,label=str(10*(i-1))+'μm')
    # plt.plot(x,y2)
    # plt.plot(x,y3)
    plt.xlabel('Pressure (Pa)',fontsize=15)
    plt.ylabel('Photophoretic force (nN) ',fontsize=15)
    plt.legend(loc='best')
    plt.grid()
    plt.show()


