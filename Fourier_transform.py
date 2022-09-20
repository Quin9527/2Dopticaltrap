import matplotlib.pyplot as plt
import matplotlib
import math
import numpy as np
import os


if __name__ == '__main__':

    def plot(y):

            plt.figure(1)
            plt.plot(y['t'], 1e3 * y['pos_z'])
            plt.xlabel('Time(s)')
            plt.ylabel('Position_z (mm)')
            plt.grid()
            plt.figure(2)
            plt.plot(y['t'], 1e3 * y['vels_z'])
            plt.xlabel('Time (s)')
            plt.ylabel('velocity (mm/s)')

            plt.grid()
            plt.show()

    def plot_fft(y,str):



        plt.loglog(y,label=str)

        plt.xlabel('Frequency (Hz)',fontsize=14)
        plt.ylabel('Power spectral density  (W/Hz)',fontsize=14)
        plt.grid()
    matplotlib.style.use('fast')
    plt.figure(1)

    def plot_spectrum(data,str):

        x = [ele.real for ele in data]

        # extract imaginary part

        y = [ele.imag for ele in data]
        Aw = np.sqrt(np.power(x,2)+np.power(y,2))
        # plt.scatter(x, y,marker='x')
        #
        # plt.ylabel('Imaginary')
        #
        # plt.xlabel('Real')
        plt.loglog(Aw,label=str)

        plt.ylabel('Frequency spectrum',fontsize=15)

        plt.xlabel('Frequency (Hz)',fontsize=15)

    def plot_psd(data,str,num_fft):

        x = [ele.real for ele in data]

        # extract imaginary part

        y = [ele.imag for ele in data]
        Aw = np.sqrt(np.power(x, 2) + np.power(y, 2))
        Aw = Aw
        Aw = Aw/num_fft
        # plt.scatter(x, y,marker='x')

        # plt.ylabel('Imaginary')

        # plt.xlabel('Real')

        plt.loglog(1e6*Aw,label=str)

        plt.ylabel('Power spectral density',fontsize=15)

        plt.xlabel('Frequency (Hz)',fontsize=15)



    # plt.figure(1)
    num_fft = 100000
    matplotlib.rcParams['agg.path.chunksize'] = 10000
    os.chdir('E:\Imperial\Project\Simulation results\Frequency')


    # filename = ('fs_final_timestep5e-4_1D.npz')
    # tempfile1 = np.load(filename,allow_pickle=True)
    # y1= np.fft.rfft(tempfile1['pos_z'],num_fft)
    # filename = ('fs_final_timestep5e-5_1D.npz')
    # tempfile2 = np.load(filename,allow_pickle=True)
    # filename = ('fs_final_timestep5e-6_1D.npz')
    # tempfile3 = np.load(filename,allow_pickle=True)

    # y2 = np.fft.rfft(tempfile2['pos_z'],num_fft)

    # y3 = np.fft.rfft(tempfile3['pos_z'],num_fft)

   
    # plot_psd(y3,'τ/100',100000)

    # plot_psd(y2,'τ/10',100000)

    # plot_psd(y1,'τ',100000)



    # plt.xlim([1e-1,1e6])
    # plt.legend(loc='best')

    #
    # os.chdir('E:\Imperial\Project\Simulation results\Final results')

    #
    plt.figure(1)
    filename = ('atm_1D.npz')
    tempfile3 = np.load(filename,allow_pickle=True)
    y = np.fft.rfft(tempfile3['pos_z'],num_fft)
    plot_spectrum(y,'(a)' )

    filename = ('lp_1D.npz')
    tempfile4 = np.load(filename,allow_pickle=True)
    y = np.fft.rfft(tempfile4['pos_z'],num_fft)
    plot_spectrum(y, '(b)' )
    filename = ('phm_1D.npz')

    tempfile5 = np.load(filename,allow_pickle=True)
    y = np.fft.rfft(tempfile5['pos_z'],num_fft)
    plot_spectrum(y, '(c)' )
    # filename = ('fs_final_nonoise_ts5e-3_1D.npz')
    # tempfile6 = np.load(filename,allow_pickle=True)
    # plot_fft(np.power(np.fft.rfft(tempfile6['pos_z'],num_fft),2)/num_fft, 'No brownian motion 10τ')
    # filename = ('RK23_timesteps_1D.npz')
    # num_fft = 10000
    # tempfile5 = np.load(filename,allow_pickle=True)
    # y = np.fft.rfft(tempfile5['pos_z'],num_fft)
    # plot_spectrum(y, 'RK23 method, no brownian motion ' )
    # tempfile7 = np.load(filename,allow_pickle=True)
    # plot_fft(np.power(np.fft.rfft(tempfile7['pos_z'],num_fft),2)/num_fft, 'No brownian motion RK23')
    plt.grid()
    plt.legend(loc='best')



    plt.show()
