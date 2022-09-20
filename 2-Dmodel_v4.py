import numpy as np
import math
import matplotlib.pylab as plt
import scipy.io
import time
import Visc_lp as vi
import Opt_f_v3 as Optf
import Buoyancy_v2 as buo
import Ph_f_v1 as phf
import scipy.integrate as inte
import sympy
from tempfile import TemporaryFile
from matplotlib.animation import FuncAnimation
import os

# Define global physical constants
Boltzmann = 1.38064852e-23
g = 9.8
miu = 1.81e-5    #viscosity for air in 300K
T = 300
pressure = []
time_count = 0

def integrate(pos, vels, forces, mass, dt,radius,pres,noise_switch,temp,onlybrownian):

    vis = vi.vis_(radius,pres).vis_()
    tau = 1/vis
    noise_ = noise(radius,pres,noise_switch,temp,onlybrownian)
    t = tau/mass
    if onlybrownian:

        vels =  noise_ * math.sqrt(tau/(2*mass) * (
                1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt)

    elif noise_switch:

        vels = noise_ * math.sqrt(tau/(2*mass) * (
                1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt) + forces *  dt / mass
    else:

        vels = forces *  dt / mass + vels * math.e**(-(1/(tau*mass))*dt)

    pos += vels * dt

    return [pos,vels]

def pres_cal(Inpres,Endpres,tim,t):

    Inpres_lg = math.log10(Inpres)
    Endpres_lg = math.log10(Endpres)
    pres_ = 10**(Inpres_lg+t*(Endpres_lg-Inpres_lg)/tim)
    return pres_

def integrate_RK4_z(  pos_z, pos_r,tim, mass, vels, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres):

    args = (temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres)
    t = (step-1)*dt

    def v_cal(t , pos_z, pos_r,  vels,mass, *args):

        vis = vi.vis_(r, pres_cal(Inpres,Endpres,tim,t)).vis_()
        tau = 1 / vis
        noise_ = noise(r, pres_cal(Inpres,Endpres,tim,t), switch4, temp, onlynoise)
        if onlynoise:

            vels =  noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt)
            # dv = (noise_ * math.sqrt(tau/(2*mass) * (
            #         1-math.e**(-(2/(tau*mass)*dt))))+vels*(math.e**(-(1/(tau*mass))*dt)-1))/dt
        elif switch4:

            vels = noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt) + computeForce_z(tim, pos_z,
                         pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
                                rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise) *  dt / mass

            # dv = (noise_ * math.sqrt(tau/(2*mass) * (
            #         1-math.e**(-(2/(tau*mass)*dt))))+vels*(math.e**(-(1/(tau*mass))*dt)-1))/dt+ computeForce_z(tim, pos_z,
            #             pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
            #                    rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise) / mass
        else:

            vels = computeForce_z(tim, pos_z,
                        pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
                                rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise) *  dt / mass + vels * math.e**(-(1/(tau*mass))*dt)
            # dv = (vels*(math.e**(-(1/(tau*mass))*dt)-1))/dt+ computeForce_z(tim, pos_z, pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
            #             rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise) / mass
        return vels

        # k1 = dv_cal( tim, pos_z, pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
        #              rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim)
        # k2 = dv_cal( tim, pos_z, pos_r, mass, vels+dt/2*k1, temp, r, pres_cal(Inpres,Endpres,tim,t+dt/2), bw,
        #              rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim)
        # k3 = dv_cal(tim, pos_z, pos_r, mass, vels + dt / 2 * k2, temp, r, pres_cal(Inpres,Endpres,tim,t+dt/2), bw,
        #             rayleigh, power, switch1, switch2, switch3, switch4, onedim, dt, onlynoise, tim)
        # k4 = dv_cal(tim, pos_z, pos_r, mass, vels + dt * k3, temp, r, pres_cal(Inpres,Endpres,tim,t+dt/2), bw,
        #             rayleigh, power, switch1, switch2, switch3, switch4, onedim, dt, onlynoise, tim)
        # vels_1 = vels + dt/6*(k1+2*k2+2*k3+k4)
        # return vels_1

    args = (mass, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim,step, Inpres,Endpres)
    k1 = v_cal(t , pos_z, pos_r, vels,*args)
    k2 = v_cal(t+dt/2 , pos_z+dt/2*k1, pos_r, vels, *args)
    k3 = v_cal(t+dt/2 , pos_z+dt/2*k2, pos_r, vels, *args)
    k4 = v_cal(t +dt, pos_z+dt*k3, pos_r, vels, *args)
    pos_1 = pos_z + dt/6*(k1+2*k2+2*k3+k4)

    return [pos_1,1/6*(k1+2*k2+2*k3+k4)]

def integrate_RK4_r(  pos_z, pos_r,tim, mass, vels, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres):

    args = (temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres)
    t = (step-1)*dt

    def v_cal(t , pos_z, pos_r,  vels,mass, *args):

        vis = vi.vis_(r, pres_cal(Inpres,Endpres,tim,t)).vis_()
        tau = 1 / vis
        noise_ = noise(r, pres_cal(Inpres,Endpres,tim,t), switch4, temp, onlynoise)
        args = ( bw,rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise)
        if onlynoise:

            vels =  noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt)
            # dv = (noise_ * math.sqrt(tau/(2*mass) * (
            #         1-math.e**(-(2/(tau*mass)*dt))))+vels*(math.e**(-(1/(tau*mass))*dt)-1))/dt
        elif switch4:

            vels = noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt) + computeForce_r(tim, pos_z,
                         pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t),*args) *  dt / mass

            # dv = (noise_ * math.sqrt(tau/(2*mass) * (
            #         1-math.e**(-(2/(tau*mass)*dt))))+vels*(math.e**(-(1/(tau*mass))*dt)-1))/dt+ computeForce_z(tim, pos_z,
            #             pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
            #                    rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise) / mass
        else:

            vels = computeForce_r(tim, pos_z,
                        pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t),*args) *  dt / mass + vels * math.e**(-(1/(tau*mass))*dt)
            # dv = (vels*(math.e**(-(1/(tau*mass))*dt)-1))/dt+ computeForce_z(tim, pos_z, pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
            #             rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise) / mass
        return vels

        # k1 = dv_cal( tim, pos_z, pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), bw,
        #              rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim)
        # k2 = dv_cal( tim, pos_z, pos_r, mass, vels+dt/2*k1, temp, r, pres_cal(Inpres,Endpres,tim,t+dt/2), bw,
        #              rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim)
        # k3 = dv_cal(tim, pos_z, pos_r, mass, vels + dt / 2 * k2, temp, r, pres_cal(Inpres,Endpres,tim,t+dt/2), bw,
        #             rayleigh, power, switch1, switch2, switch3, switch4, onedim, dt, onlynoise, tim)
        # k4 = dv_cal(tim, pos_z, pos_r, mass, vels + dt * k3, temp, r, pres_cal(Inpres,Endpres,tim,t+dt/2), bw,
        #             rayleigh, power, switch1, switch2, switch3, switch4, onedim, dt, onlynoise, tim)
        # vels_1 = vels + dt/6*(k1+2*k2+2*k3+k4)
        # return vels_1

    args = (mass, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim,step, Inpres,Endpres)
    k1 = v_cal(t , pos_z, pos_r, vels,*args)
    k2 = v_cal(t+dt/2 , pos_z, pos_r+dt/2*k1, vels, *args)
    k3 = v_cal(t+dt/2 , pos_z, pos_r+dt/2*k2, vels, *args)
    k4 = v_cal(t +dt, pos_z, pos_r+dt*k3, vels, *args)
    pos_1 = pos_r + dt/6*(k1+2*k2+2*k3+k4)

    return [pos_1,1/6*(k1+2*k2+2*k3+k4)]

def integrate_RK2_z(  pos_z, pos_r,tim, mass, vels, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres):

    args = (temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres)
    t = (step-1)*dt

    def v_cal(t , pos_z, pos_r,  vels,mass, *args):

        vis = vi.vis_(r, pres_cal(Inpres,Endpres,tim,t)).vis_()
        tau = 1 / vis
        noise_ = noise(r, pres_cal(Inpres,Endpres,tim,t), switch4, temp, onlynoise)
        args1 = ( bw,rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise)
        if onlynoise:

            vels =  noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt)

        elif switch4:

            vels = noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt) + \
                   computeForce_z(tim, pos_z,
                         pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t),*args1) *  dt / mass

        else:

            vels = computeForce_z(tim, pos_z,
                        pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t), *args1) \
                   *  dt / mass + vels * math.e**(-(1/(tau*mass))*dt)

        return vels

    args = (mass, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim, Inpres,Endpres)
    k1 = v_cal(t , pos_z, pos_r, vels,*args)
    k2 = v_cal(t+dt/2 , pos_z+dt*k1/2, pos_r, vels,*args)
    pos_1 = pos_z + dt/2*(k1+k2)

    return [pos_1,1/2*(k1+k2)]
def integrate_RK2_r(  pos_z, pos_r,tim, mass, vels, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres):

    args = (temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,step, Inpres,Endpres)
    t = (step-1)*dt

    def v_cal(t , pos_z, pos_r,  vels,mass, *args):


        vis = vi.vis_(r, pres_cal(Inpres,Endpres,tim,t)).vis_()
        tau = 1 / vis
        noise_ = noise(r, pres_cal(Inpres,Endpres,tim,t), switch4, temp, onlynoise)
        args1 = ( bw,rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise)
        if onlynoise:

            vels =  noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt)

        elif switch4:

            vels = noise_ * math.sqrt(tau/(2*mass) * (
                    1-math.e**(-(2/(tau*mass)*dt)))) + vels * math.e**(-(1/(tau*mass))*dt) + computeForce_r(tim, pos_z,
                         pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t),*args1) *  dt / mass


        else:

            vels = computeForce_r(tim, pos_z,
                        pos_r, mass, vels, temp, r, pres_cal(Inpres,Endpres,tim,t),
                                  *args1) *  dt / mass + vels * math.e**(-(1/(tau*mass))*dt)

        return vels

    args = (mass, temp, r, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise,tim,step, Inpres,Endpres)
    k1 = v_cal(t , pos_z, pos_r, vels,*args)
    k2 = v_cal(t+dt/2 , pos_z+dt*k1/2, pos_r, vels, *args)
    pos_1 = pos_r + 1/2*(k1+k2)

    return [pos_1,1/2*(k1+k2)]

def noise(r,pres,switch,temp,onlybrwonian):

    vis = vi.vis_(r,pres).vis_()
    damp = vis
    sigma = math.sqrt(2.0 * damp * temp * Boltzmann )
    noise_term = 0
    noise_term = np.random.randn() * sigma
    # print("v:%e" % noise_term + " m/s  ",end=' ')
    if switch or onlybrwonian:
        return noise_term
    else:
        return 0

def computeForce_z( tim, pos_z, pos_r, mass, vels_z, temp, r, pres, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise):

    if onlynoise:
        return 0
    vis = vi.vis_(r,pres).vis_()
    force = 0
    force += buo.buo(r,pres).buo()                                          #Buoyancy
    force += -g * mass                                                      #Gravity
    # if switch1:
    #  force += -vis * vels_z                                                   #Viscos drag
    if switch1:
        force += Optf.optf(bw,power,pos_r,pos_z,r,rayleigh,True,False,onedim).cal_f()   #Opitcal force
    if switch3:
        force += phf.ph_f(power,bw,pres,pos_z,r,pos_r).total_fp(50)                      #photophoretic force

    return force

def computeForce_r( tim, pos_z, pos_r, mass, vels_r, temp, r, pres, bw,
                 rayleigh, power,switch1,switch2,switch3,switch4,onedim,dt,onlynoise):

    # vis = vi.vis_(r,pres).vis_()
    if onlynoise:
        return 0
    force = 0
    # if switch1:
    #  force += -vis * vels_r                                                   #Viscos drag
    if switch1:
     force += -Optf.optf(bw,power,pos_r,pos_z,r,rayleigh,False,True,onedim).cal_f()   #Opitcal force




    return force


def run(**args):

    temp = args['temp']
    radius, wavelength = args['radius'],args['wavelength']
    mass = args['density_of_particle']*4/3*math.pi*radius**3
    bw, M2,power = args['beamwaist'], args['M2'], args['power']
    OF,PF,VD,Noi = args['Opticalforce'], args['Photophoreticforce'] , args['Viscosdrag'], args['Noise']
    Inpos_z, Inv_z, Inpos_r, Inv_r = args['Initialposition_z'], args['Initialvelocity_z'],\
                     args['Initialposition_r'], args['Initialvelocity_r']
    Inpres, Endpres, tim, Twodim = args['Initialpressure'], args['Endpressure']\
        ,args['time'],args['2-dimension']
    dt, onlynoise, method = args['dt'], args['onlybrownian'], args['method']

    nsteps = int(tim/dt)
    global pressure, percent_
    Inpres_lg = math.log10(Inpres)
    Endpres_lg = math.log10(Endpres)
    rayleigh = math.pi * bw ** 2 / (wavelength*M2)

    press = np.logspace(Inpres_lg,Endpres_lg, nsteps,base=10)
    # Dimension
    dim = 2
    pos = np.zeros((dim,nsteps))
    vels = np.zeros((dim,nsteps))
    pos.astype(np.float64)
    vels.astype(np.float64)
    # Initial position
    pos[0,0] = Inpos_z
    # Initial velocity
    vels[0,0] = Inv_z
    pos[1, 0] = Inpos_r
    # Initial velocity
    vels[1, 0] = Inv_r

    step = 0
    output = []

    while step < nsteps:

        step += 1
        args = (temp, radius, bw,
                rayleigh, power, OF, PF, VD, Noi, not Twodim, dt, onlynoise, step, Inpres, Endpres)
        args_eu = (bw, rayleigh, power, VD, OF, PF, Noi, not Twodim, dt, onlynoise)
        if method == 2: # RK4

            # Compute all forces
            if (step == 1):

             tem = integrate_RK4_z(pos[dim - 2, step - 1],pos[dim - 1, step - 1],   tim,  mass, vels[dim - 2, step - 1], *args)
            else:
             tem = integrate_RK4_z(pos[dim - 2, step - 2],pos[dim - 1, step - 2],   tim,  mass, vels[dim - 2, step - 2], *args)
            pos[dim-2,step-1] = tem[0]
            vels[dim-2,step-1] = tem[1]


            if Twodim:

                if (step == 1):

                    tem = integrate_RK4_r(pos[dim - 2, step - 1], pos[dim - 1, step - 1], tim, mass,
                                          vels[dim - 1, step - 1], *args)
                else:
                    tem = integrate_RK4_r(pos[dim - 2, step - 2], pos[dim - 1, step - 2], tim, mass,
                                          vels[dim - 1, step - 2], *args)
                pos[dim - 1, step - 1] = tem[0]
                vels[dim - 1, step - 1] = tem[1]
        if method == 1: # RK2

            # Compute all forces
            if (step == 1):

                tem = integrate_RK2_z(pos[dim - 2, step - 1],pos[dim - 1, step - 1],   tim,  mass, vels[dim - 2, step - 1], *args)
            else:
                tem = integrate_RK2_z(pos[dim - 2, step - 2],pos[dim - 1, step - 2],   tim,  mass, vels[dim - 2, step - 2], *args)
            pos[dim-2,step-1] = tem[0]
            vels[dim-2,step-1] = tem[1]


            if Twodim:

                if (step == 1):

                    tem = integrate_RK2_r(pos[dim - 2, step - 1], pos[dim - 1, step - 1], tim, mass,
                                          vels[dim - 1, step - 1], *args)
                else:
                    tem = integrate_RK2_r(pos[dim - 2, step - 2], pos[dim - 1, step - 2], tim, mass,
                                          vels[dim - 1, step - 2], *args)
                pos[dim - 1, step - 1] = tem[0]
                vels[dim - 1, step - 1] = tem[1]
        if method == 0: #Euler
            # Compute all forces
            if (step == 1):
                forces = computeForce_z(tim, pos[dim - 2, step - 1], pos[dim - 1, step - 1], mass,
                                        vels[dim - 2, step - 1], temp,
                                        radius, press[step - 1], *args_eu)
                tem = integrate(pos[dim - 2, step - 1], vels[dim - 2, step - 1], forces, mass, dt,
                                radius, press[step - 2], Noi, temp, onlynoise)
            else:
                forces = computeForce_z(tim, pos[dim - 2, step - 2], pos[dim - 1, step - 2], mass,
                                        vels[dim - 2, step - 2], temp,
                                        radius, press[step - 2], *args_eu)
                # Move the system in time
                tem = integrate(pos[dim - 2, step - 2], vels[dim - 2, step - 2], forces, mass, dt, radius,
                                press[step - 2], Noi, temp, onlynoise)
            pos[dim - 2, step - 1] = tem[0]
            vels[dim - 2, step - 1] = tem[1]

            # if tem[0] < -1.5*abs(Inpos_z) and tem[1]< 0:
            #     print("particle cannot be trapped")
            #     break

            if Twodim:

                if (step == 1):
                    forces = computeForce_r(tim, pos[dim - 2, step - 1], pos[dim - 1, step - 1], mass,
                                            vels[dim - 1, step - 1], temp, radius,
                                            press[step - 1], *args_eu)
                    tem = integrate(pos[dim - 1, step - 1], vels[dim - 1, step - 1], forces, mass, dt, radius,
                                    press[step - 2], Noi, temp, onlynoise)
                else:
                    forces = computeForce_r(tim, pos[dim - 2, step - 2], pos[dim - 1, step - 2], mass,
                                            vels[dim - 1, step - 2], temp, radius,
                                            press[step - 2], *args_eu)
                    # Move the system in time
                    tem = integrate(pos[dim - 1, step - 2], vels[dim - 1, step - 2], forces, mass, dt, radius,
                                    press[step - 1], Noi, temp, onlynoise)
                pos[dim - 1, step - 1] = tem[0]
                vels[dim - 1, step - 1] = tem[1]
        if  method != 1 and method !=0 and method !=2:
            print("Error: None of the solver is seleted")
            break
        precent = step/nsteps*100
        print("\rr:%e" % pos[dim - 1, step - 1] + " m     "
                                  "z:%e" % pos[dim - 2, step - 1] + " m      V_r:%E"%vels[dim - 1, step -
                                            1]+" m/s V_z:%e"% vels[dim - 2, step - 1]+
                                                   " m/s   pressure:%.3e" % press[step - 1] + "Pa    "
                                                           "progress:%.2f" % precent+"%", end=' ')
    output = np.ndarray((6,nsteps))
    output[0,:] = pos[0,:]
    output[1,:] = vels[0,:]
    output[2,:] = press
    output[3,:] = np.linspace(0,tim,nsteps)
    if Twodim:
        output[4,:] = pos[1,:]
        output[5,:] = vels[1,:]
    global time_count
    time_count = 0

    return output


if __name__ == '__main__':


    def plot(y,Twodim):

        if Twodim:
            plt.figure(1)
            plt.plot(y[3,:], 1e3*y[0,:])
            plt.xlabel('Time(s)')
            plt.ylabel('Position_z (mm)')
            plt.grid()
            plt.figure(2)
            plt.plot(y[3,:], 1e6*y[4,:])
            plt.xlabel('Time (s)')
            plt.ylabel('Position_r (μm)')
            plt.grid()
            plt.figure(3)
            plt.plot(y[3,:], 1e3*y[5, :])
            plt.xlabel('Time(s)')
            plt.ylabel('Velocity_r (mm/s)')
            plt.grid()
            plt.figure(4)
            plt.plot(y[3,:], 1e3*y[1, :])
            plt.xlabel('Time (s)')
            plt.ylabel('velocity_z (mm/s)')
            plt.grid()
            plt.figure(5)
            plt.plot(1e6*y[4, :], 1e3*y[0, :])
            plt.xlabel('Position_r (μm)')
            plt.ylabel('Position_z (mm)')
            plt.grid()
            plt.show()
        else:
            # plt.figure(1)
            plt.plot(y[3,:], 1e3*y[0, :])
            plt.xlabel('Time(s)')
            plt.ylabel('Position_z (mm)')
            # plt.figure(2)
            # plt.plot(y[3,:], 1e6*y[1, :])
            # plt.xlabel('Time (s)')
            # plt.ylabel('velocity (μm/s)')
            # plt.show()

    def savefile(output1,Twodim,filename):

        if Twodim:
            np.savez(filename+'_2D', t= output1[3,:], pos_r=output1[4,:], vels_r=output1[5,:],
                     pos_z=output1[0,:], vels_z=output1[1,:])
        else:
            np.savez(filename+'_1D', t=output1[3,:],
                     pos_z=output1[0,:], vels_z=output1[1,:])


    def setting(optf_switch,phf_switch,noise_switch,onlybrownian,two_dimension_switch,timestep,method):

        params = {
            'temp': 300,
            'density_of_particle': 960,
            'radius': 5e-6,
            'M2' : 1.085,
            'beamwaist' : 1.866e-6,
            'wavelength':532e-9,
            'time':6,
            'power': 0.4,
            'Opticalforce':optf_switch,
            'Photophoreticforce':phf_switch,
            'Viscosdrag':True,
            'Noise':noise_switch,
            'Initialposition_z':1e-3,
            'Initialvelocity_z':0,
            'Initialposition_r': 0,
            'Initialvelocity_r': 0,
            'Initialpressure':101300,
            'Endpressure':1.33e-4,
            '2-dimension':two_dimension_switch,
            'dt':timestep,
            'onlybrownian':False,
            'method': method
        }
        return params

    os.chdir('E:\Imperial\Project\Simulation results\Final results')

    Twodim = False
    set1 = setting(True,True,True,False,Twodim,5e-4,0)
    timestart = time.time()
    output1 = run(**set1)
    timeend = time.time()
    timecost = timeend-timestart
    print("\ntime cost:%.2f"%timecost+"s")
    # savefile(output1,Twodim,'Model4_Eu_τ')
    set2 = setting(True,True,True,False,Twodim,1e-3,1)
    timestart = time.time()
    output2 = run(**set2)
    timeend = time.time()
    timecost = timeend-timestart
    print("\ntime cost:%.2f"%timecost+"s")
    # savefile(output2,Twodim,'Model4_RK2_2τ')
    set3 = setting(True,True,True,False,Twodim,5e-3,2)
    timestart = time.time()
    output3 = run(**set3)
    timeend = time.time()
    timecost = timeend-timestart
    # print("\ntime cost:%.2f"%timecost+"s")
    savefile(output3,Twodim,'Model4_RK4_10τ')
    plot(output1,Twodim)
    plot(output2,Twodim)
    plot(output3,Twodim)
    plt.show()
