import numpy as np
import math
from math import e
from math import sqrt
from math import cos
from math import sin
import matplotlib.pylab as plt

import time
import Visc_lp as vi
import Opt_f_v3 as Optf
import Buoyancy_v2 as buo
import Ph_f_v1 as phf



import os

# Define global physical constants
Boltzmann = 1.38064852e-23
g = 9.8
miu = 1.81e-5  # viscosity for air in 300K
T = 300
pressure = []
time_count = 0


def integrate(mindt,z_direction,t, pos, vels, forces, mass, dt, radius, pres, noise_switch, temp, onlybrownian):
    vis = vi.vis_(radius, pres).vis_()
    tau = 1 / vis
    noise_ = noise(radius, pres, noise_switch, temp, onlybrownian)

    if onlybrownian:

        vels = noise_ * sqrt(tau / (2 * mass) * (
                1 - e ** (-(2 / (tau * mass) * dt)))) + vels * e ** (-(1 / (tau * mass)) * dt)

    elif noise_switch:

        vels = noise_ * sqrt(tau / (2 * mass) * (
                1 - e ** (-(2 / (tau * mass) * dt)))) + vels * e ** (
                           -(1 / (tau * mass)) * dt) + forces * dt / mass
    else:

        vels = forces * dt / mass + vels * e ** (-(1 / (tau * mass)) * dt)

    pos += vels * dt
    if z_direction:
        v = abs(vels)
        if v > 3e-3:
            dt = mindt
        elif v <=3e-3 and v > 1e-4:
            dt = 5e-4
        elif v <=1e-4 and v > 5e-5:
            dt = 1e-3
        else:
            dt = 2e-3

    return [pos, vels, dt]

def pres_cal(Inpres,Endpres,tim,t):

    Inpres_lg = math.log10(Inpres)
    Endpres_lg = math.log10(Endpres)
    pres_ = 10**(Inpres_lg+t*(Endpres_lg-Inpres_lg)/tim)
    return pres_

def noise(r, pres, switch, temp, onlybrwonian):
    vis = vi.vis_(r, pres).vis_()
    damp = vis
    sigma = math.sqrt(2.0 * damp * temp * Boltzmann)
    noise_term = 0
    noise_term = np.random.randn() * sigma
    # print("v:%e" % noise_term + " m/s  ",end=' ')
    if switch or onlybrwonian:
        return noise_term
    else:
        return 0


def computeForce_z(tim, pos_z, pos_r, mass, vels_z, temp, r, pres, bw,
                   rayleigh, power, switch1, switch2, switch3, switch4, onedim, dt, onlynoise,tilt):
    if onlynoise:
        return 0
    force = 0
    force += buo.buo(r, pres).buo()  # Buoyancy
    force += -g * mass   # Gravity

    if switch1:
        force += Optf.optf(bw, power, pos_r, pos_z, r, rayleigh, True, False, onedim).cal_f()* cos(tilt)
        force += Optf.optf(bw, power, pos_r, pos_z, r, rayleigh, False, True, onedim).cal_f() * sin(tilt)# Opitcal force
    if switch3:
        force += phf.ph_f(power, bw, pres, pos_z, r, pos_r).total_fp(50)* cos(tilt)  # photophoretic force

    return force


def computeForce_r(tim, pos_z, pos_r, mass, vels_r, temp, r, pres, bw,
                   rayleigh, power, switch1, switch2, switch3, switch4, onedim, dt, onlynoise,tilt):
    if onlynoise:
        return 0
    force = 0
    # if not switch4:
                                              #Viscos drag
    if switch1:
        force += -Optf.optf(bw, power, pos_r, pos_z, r, rayleigh, False, True, onedim).cal_f() * cos(tilt)  # Opitcal force
        force += Optf.optf(bw, power, pos_r, pos_z, r, rayleigh, True, False, onedim).cal_f()* sin(tilt)
    if switch3:
        force += phf.ph_f(power, bw, pres, pos_z, r, pos_r).total_fp(50)* sin(tilt)  # photophoretic force
    return force


def run(**args):
    temp = args['temp']
    radius, wavelength = args['radius'], args['wavelength']
    mass = args['density_of_particle'] * 4 / 3 * math.pi * radius ** 3
    bw, M2, power = args['beamwaist'], args['M2'], args['power']
    OF, PF, VD, Noi = args['Opticalforce'], args['Photophoreticforce'], args['Viscosdrag'], args['Noise']
    Inpos_z, Inv_z, Inpos_r, Inv_r = args['Initialposition_z'], args['Initialvelocity_z'], \
                                     args['Initialposition_r'], args['Initialvelocity_r']
    Inpres, Endpres, tim, Twodim = args['Initialpressure'], args['Endpressure'] \
        , args['time'], args['2-dimension']
    dt, onlynoise,tilt = args['mindt'], args['onlynoise'],args['tilt']

    if Inpres == 0 or Endpres == 0:
        print("Error: pressure is 0 Pa")
        return

    nsteps = int(tim / dt)
    global pressure, percent_
    Inpres_lg = math.log10(Inpres)
    Endpres_lg = math.log10(Endpres)
    rayleigh = math.pi * bw ** 2 / (wavelength * M2)
    press = np.logspace(Inpres_lg, Endpres_lg, nsteps, base=10)
    # Dimension

    mindt = dt
    # Initial position
    pos_z = []
    pos_r = []
    vels_z = []
    vels_r = []
    timecount = []
    pos_z.append(Inpos_z)
    pos_r.append(Inpos_r)
    # Initial velocity
    vels_z.append(Inv_z)
    vels_r.append(Inv_r)


    rayleigh = math.pi * bw ** 2 / (wavelength*M2)
    t = 0

    output = []
    timecount.append(0)
    tem_r = np.zeros(3)
    tem = np.zeros(3)
    tem[0]=pos_z[0]
    tem_r[0]=pos_r[0]
    tem[1]=vels_z[0]
    tem_r[1]=vels_r[0]

    while True:

        t = t + dt

        args_eu = (bw, rayleigh, power, VD, OF, PF, Noi, not Twodim, dt, onlynoise,tilt)

        forces = computeForce_z(tim, tem[0], tem_r[0], mass,
                                tem[1], temp,
                                radius, pres_cal(Inpres,Endpres,tim,t), *args_eu)
        # Move the system in time
        tem = integrate(mindt,True,t,tem[0], tem[1], forces, mass, dt, radius,
                        pres_cal(Inpres,Endpres,tim,t), Noi, temp, onlynoise)
        pos_z.append(tem[0])
        vels_z.append(tem[1])


        if Twodim:

            forces = computeForce_r(tim, tem[0], tem_r[0], mass,
                                    tem_r[1], temp, radius,
                                    pres_cal(Inpres,Endpres,tim,t), *args_eu)
            # Move the system in time
            tem_r = integrate(mindt,False,t,tem_r[0], tem_r[1], forces, mass, dt, radius,
                            pres_cal(Inpres,Endpres,tim,t), Noi, temp, onlynoise)
            pos_r.append(tem_r[0])
            vels_r.append(tem_r[1])
        else:
            tem_r[0]=0
        dt = tem[2]
        timecount.append(t)
        if t > tim:
            break
        precent = t / tim * 100
        print("\rr:%e" % tem_r[0] + " m     "
                                                  "z:%e" % tem[0] + " m      V_r:%E" % tem_r[1] + " m/s V_z:%e" % tem[1] +
              " m/s   pressure:%.3e" % pres_cal(Inpres,Endpres,tim,t) + "Pa    "
                                                         "progress:%.2f" % precent + "% "   +"  time step:%e" % tem[2] +" t:%e"%t, end=' ')

    output0 = pos_z
    output.append(np.array(output0))
    output1 = vels_z
    output.append(np.array(output1))
    output2 = press
    output.append(np.array(output2))
    output3 = timecount
    output.append(np.array(output3))
    if Twodim:
        output4 = pos_r
        output.append(np.array(output4))
        output5 = vels_r
        output.append(np.array(output5))
    global time_count
    time_count = 0

    return output


if __name__ == '__main__':

    def plot(y, Twodim, str):

        if Twodim:
            plt.figure(1)
            plt.plot(y[3], 1e3 * y[0])
            plt.xlabel('Time(s)')
            plt.ylabel('Position_z (mm)')
            plt.grid()
            plt.figure(2)
            plt.plot(y[3], 1e6 * y[4])
            plt.xlabel('Time (s)')
            plt.ylabel('Position_r (μm)')
            plt.grid()
            plt.figure(3)
            plt.plot(y[3], 1e3 * y[5])
            plt.xlabel('Time(s)')
            plt.ylabel('Velocity_r (mm/s)')
            plt.grid()
            plt.figure(4)
            plt.plot(y[3], 1e3 * y[1])
            plt.xlabel('Time (s)')
            plt.ylabel('velocity_z (mm/s)')
            plt.grid()
            plt.figure(5)
            plt.plot(1e6 * y[4], 1e3 * y[0])
            plt.xlabel('Position_r (μm)')
            plt.ylabel('Position_z (mm)')
            plt.grid()
            plt.show()
        else:
            # plt.figure(1)
            plt.plot(y[3], 1e3 * y[0], label=str)
            plt.xlabel('Time(s)')
            plt.ylabel('Position_z (mm)')
            # plt.figure(2)
            # plt.plot(y[3,:], 1e6*y[1, :])
            # plt.xlabel('Time (s)')
            # plt.ylabel('velocity (μm/s)')
            # plt.show()


    def savefile(output1, Twodim, filename):

        if Twodim:
            np.savez(filename + '_2D', t=output1[3], pos_r=output1[4], vels_r=output1[5],
                     pos_z=output1[0], vels_z=output1[1])
        else:
            np.savez(filename + '_1D', t=output1[3],
                     pos_z=output1[0], vels_z=output1[1])


    def setting(optf_switch, phf_switch, noise_switch, two_dimension_switch, timestep, t, z0,zr,tilt):

        params = {
            'temp': 300,
            'density_of_particle': 960,
            'radius': 5e-6,
            'M2': 1.085,
            'beamwaist': 1.866e-6,
            'wavelength': 532e-9,
            'time': t,
            'power': 0.4,
            'Opticalforce': optf_switch,
            'Photophoreticforce': phf_switch,
            'Viscosdrag': True,
            'Noise': noise_switch,
            'Initialposition_z': z0,
            'Initialvelocity_z': 0,
            'Initialposition_r': zr,
            'Initialvelocity_r': 0,
            'Initialpressure':101300,
            'Endpressure': 101300,
            '2-dimension': two_dimension_switch,
            'mindt': timestep,
            'tilt':tilt,
            'onlynoise': False
            
        }
        return params


    os.chdir('E:\Imperial\Project\Simulation results\TWODIM')

    Twodim = True
    #
    # set1 = setting(True, True,True, Twodim, 1e-4, 1, 2e-3,2e-6)
    # timestart = time.time()
    # output1 = run(**set1)
    # timeend = time.time()
    # timecost = timeend - timestart
    # print("\ntime cost:%.2f" % timecost + "s")
    # savefile(output1, Twodim, '2_medium')

    set3 = setting(True,True,True,Twodim,1e-4,0.2,2e-3,0,math.pi/36)
    timestart = time.time()
    output3 = run(**set3)
    timeend = time.time()
    timecost = timeend-timestart
    print("\ntime cost:%.2f"%timecost+"s")
    savefile(output3,Twodim,'15tilt_hp_1')
    set2 = setting(True,False,True,Twodim,1e-4,0.2,2e-3,0,math.pi/18)
    timestart = time.time()
    output2 = run(**set2)
    timeend = time.time()
    timecost = timeend-timestart
    print("\ntime cost:%.2f"%timecost+"s")
    savefile(output2,Twodim,'30tilt_hp_1')

    # set1 = setting(True,True,True,Twodim,5e-4,600,30e-3)
    # timestart = time.time()
    # output1 = run(**set1)
    # timeend = time.time()
    # timecost = timeend-timestart
    # print("\ntime cost:%.2f"%timecost+"s")
    # savefile(output1,Twodim,'Final_all')
    #
    # set3 = setting(True,False,True,Twodim,5e-4,600,30e-3)
    # timestart = time.time()
    # output3 = run(**set3)
    # timeend = time.time()
    # timecost = timeend-timestart
    # print("\ntime cost:%.2f"%timecost+"s")
    # savefile(output3,Twodim,'Final_nophf')
    # set2 = setting(True,True,False,Twodim,5e-4,600,30e-3)
    # timestart = time.time()
    # output2 = run(**set2)
    # timeend = time.time()
    # timecost = timeend-timestart
    print("\ntime cost:%.2f" % timecost + "s")
    # savefile(output2,Twodim,'Final_nonoi')
    plot(output1,Twodim,'Euler_τ')
    # plot(output2,Twodim,'Euler_τ/10')
    # plot(output3,Twodim,'Euler_5τ')
    # plot(output4,Twodim,'Euler')
    plt.legend(loc='best')
    plt.show()
