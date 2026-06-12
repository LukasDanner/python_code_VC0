import numpy as np
import qutip as qt
import sys
import shutil
import os

# parameter class
from Class_Parameters.class_parameters import Parameters

# logging class
from Logging.logging import Logging

# time grid factory
from factory_timegrid import factory_timegrid

# circuit factory
from Circuits.factory_circuits import factory_circuits

from helpers_classical import plot_N_in_1
from helpers_classical import plot_N_in_1_diff_xaxis
from helpers_classical import write_parameters_run
from helpers_classical import plot_phasespace_classical
from helpers_classical import fourier_transform
from helpers_classical import inverse_fourier_transform
from helpers_classical import current_VCO
from helpers_classical import plot_trajectory3d
from helpers_classical import plot_ReIm


from scipy.signal import find_peaks

# config-file (as arrays, if we want to switch between ini-files quickly)
config_filenames = ['configfile_VCO_clSp_RF.ini']

#relative path of root folder
output_path = '../Output/'

# create a parameter object
param_obj = Parameters()

# read input file and set up all parameters
param_obj.setup_all_parameters(config_filenames[0])

# time evolution
from TimeEvolutionClassical.propagator_classical import PropagatorClassical


# create a logging object
# which creates a root folder to store results
logging = Logging(param_obj, output_path=output_path)

shutil.copyfile(os.getcwd() +"/"+config_filenames[0], logging.root + "/configfile.txt")

for N_now in range(param_obj.N_runs):

    # setup new (full and flattened dictionary in parameter object)
    param_obj.setup_current_parameters(N_now)

    print(param_obj.params_now_flat["circuit_name"])
    
    # call the logging object, which will create a new sub-folder
    logging.notify_new_run(N_now, circuit_name=param_obj.params_now_flat["circuit_name"])

    # store the exact run parameters for this run
    write_parameters_run(param_obj.params_now, logging)

    # which parameter is varied + value of it
    vary_name = param_obj.params_now_flat["folder_vary"]
    vary_val = param_obj.params_now_flat[vary_name]

    # --- time grid
    fac_stroboscopic = 2 * np.pi if param_obj.params_now_flat["t_st"] else 1.0
    times = factory_timegrid(**param_obj.params_now_flat)
    
    # steady state times
    t_st = param_obj.params_now_flat["t_st"] * fac_stroboscopic
    i_st = np.where(times >= t_st)[0][0]  

    # zoom times
    t_z0 = param_obj.params_now_flat["t_z0"] * fac_stroboscopic
    t_zf = param_obj.params_now_flat["t_zf"] * fac_stroboscopic
    i_z0 = np.where(times >= t_z0)[0][0]  
    i_zf = np.where(times >= t_zf)[0][0] + 1

    # --- circuit
    circuit = factory_circuits(**param_obj.params_now_flat)

    # --- propagator
    propagator = PropagatorClassical(times, circuit, nsteps=1000,
                                     include_noise=circuit.include_noise,
                                     **param_obj.params_now_flat)

    pathLoad = logging.subdir + "/res.npy"

    # propagate
    if param_obj.params_now_flat["load_sol"] and os.path.exists(pathLoad):
        print("loading solution ")


        propagator.res = np.load(logging.subdir + "/res.npy")

    else:
        print("propagating ")
        propagator.propagate(COMPLEX=True)

        if param_obj.params_now_flat["store_sol"]:

            np.save(logging.subdir + "/res", propagator.res)

    # complex variable y, y' in phase-space 
    alpha = 2 * np.conjugate(propagator.res[:, 0]) + 1.0j * (np.conjugate(propagator.res[:, 1]) + propagator.res[:, 1] * np.exp(2.0j * times * circuit.Qcoil))/circuit.Qcoil

    # definition of complex magnetization
    Mxy = propagator.res[:, 2]


    # --- Fourier transform of steady state
    #res_ft = np.empty_like(propagator.res[i_st:, :])
    res_ft = np.empty([len(times[i_st:]), len(propagator.res[0])], dtype=np.complex_)

    for ll in range(len(propagator.res[0])):
        wvec, res_ft[:, ll] = fourier_transform(times[i_st:] - times[i_st], propagator.res[i_st:, ll]) 

    # frequencies used for zooming
    w_z0 = param_obj.params_now_flat["w_z0"]
    w_zf = param_obj.params_now_flat["w_zf"]
    i_fz0 = np.where(wvec >= w_z0)[0][0]  
    i_fzf = np.where(wvec >= w_zf)[0][0] 

    # index for frequency w=0
    i_f0 = np.where(wvec >= 0.0)[0][0]  

    # index for frequency w=w0
    i_fgamma = np.where(wvec >= 1.0)[0][0]  

    # index of main oscillation frequency of y
    i_osc = i_f0 + np.argmax(np.abs(res_ft[i_f0:, 0])**2)

    # index of main oscillation frequency of Mx
    i_spinosc = i_f0 + np.argmax(np.abs(res_ft[i_f0:, 2])**2)

    # main steady-state oscillation frequency
    w_osc = [wvec[i_osc]]

    print("check", np.abs(propagator.res[-1, 0]), propagator.res[-1, 0] )


    # --- Plots in time space 

    # full result of a
    data = np.zeros([3, len(times)])
    data[0] = np.real(propagator.res[:, 0])
    data[1] = np.imag(propagator.res[:, 0])
    data[2] = np.abs(propagator.res[:, 0])
    plot_N_in_1(times, data,
                xlabel='time ' + r"$\tau =  t \omega_0 / Q$", ylabel='',
                legend=[r"$Re a(\tau)$",r"$Im a(\tau)$", r"$|a(\tau)|$"],
                linestyle=["-", "-", "--"],
                fname=logging.subdir + "/a_t.pdf",
                vline=[t_st])

    # zoomed result of a
    plot_N_in_1(times[i_z0:i_zf], data[:, i_z0:i_zf],
                xlabel='time ' + r"$ t \omega_0$", ylabel='',
                legend=[r"$Re a(\tau)$",r"$Im a(\tau)$", r"$|a(\tau)|$"],
                fname=logging.subdir + "/a_t_zoom.pdf")
    
    # full result of a'
    data = np.zeros([3, len(times)])
    data[0] = np.real(propagator.res[:, 1])
    data[1] = np.imag(propagator.res[:, 1])
    data[2] = np.abs(propagator.res[:, 1])
    plot_N_in_1(times, data,
                xlabel='time ' + r"$\tau =  t \omega_0 / Q$", ylabel='',
                legend=[r"$Re a'(\tau)$",r"$Im a'(\tau)$", r"$|a'(\tau)|$"],
                fname=logging.subdir + "/adot_t.pdf",
                vline=[t_st])

    # zoomed result of a'
    plot_N_in_1(times[i_z0:i_zf], data[:, i_z0:i_zf],
                xlabel='time ' + r"$\tau = t \omega_0 / Q$",
                legend=[r"$Re a'(\tau)$",r"$Im a'(\tau)$", r"$|a'(\tau)|$"],
                fname=logging.subdir + "/adot_t_zoom.pdf")

    # M_{xy}
    plot_ReIm(times, Mxy,
              xlabel='time ' + r"$\tau = t \omega_0 / Q$",
              ylabel=r"$m_{xy}(\tau)$",
              fname=logging.subdir + "/Mxy_ReIm.pdf")    
    
    # zoomed result of M_{xy}
    plot_ReIm(times[i_z0:i_zf], Mxy[i_z0:i_zf],
              xlabel='time ' + r"$\tau = t \omega_0 / Q$",
              ylabel=r"$m_{xy}(\tau)$",
              fname=logging.subdir + "/Mxy_zoom.pdf")
        

    # M_z
    plot_N_in_1(times, [np.real(propagator.res[:, 3])],
              xlabel='time ' + r"$\tau =  t \omega_0/Q$", 
              ylabel=r"$m_{z}(\tau)$",
              fname=logging.subdir + "/Mz_t.pdf")    
    
    # magnetic field Bxy
    data = np.zeros([3, len(times)])
    data[0] = np.real(propagator.res[:, 4])
    data[1] = np.imag(propagator.res[:, 4])
    data[2] = np.abs(propagator.res[:, 4])    
    plot_N_in_1(times, data,
                xlabel='time ' + r"$\tau = t \omega_0 / Q$",
                legend=[r"$Re B_{xy}$",r"$Im B_{xy}$", r"$|B_{xy}|$"],
                fname=logging.subdir + "/Bxy_t.pdf")
    

    # --- plots in Fourier space

    # a(w)
    plot_N_in_1(wvec, [np.log10(np.abs(res_ft[:, 0])**2)],
            xlabel='frequency ' + r"$Q \omega / \omega_0$",
            legend=[r"$log_{10}\left[|a(\omega)|^2\right]$"],
            xlim=[w_z0, w_zf],
            fname=logging.subdir + "/a_FT.pdf",
            vline=[0.0])

    # M_xy(w)
    plot_N_in_1(wvec, [np.log10(np.abs(res_ft[:,2])**2)],
                xlim=[w_z0, w_zf],
                xlabel='frequency ' + r"$Q \omega / \omega_0$", ylabel='',
                legend=[r"$log_{10}\left[|M_x(\omega)|^2\right]$"],
                linestyle = ["--", "-."],
                vline=[0.0], 
                fname=logging.subdir + "/Mxy_FT.pdf")

    # Bxy
    plot_N_in_1(wvec, [np.log10(np.abs(res_ft[:,4])**2)],
                xlabel='frequency ' + r"$Q \omega / \omega_0$", ylabel='',
                legend=[r"$log_{10}\left[|B_{xy}(\omega)/\tilde{B}|^2\right]$"],
                xlim=[-3.2, 6.2],
                vline=[0.0], 
                fname=logging.subdir + "/Bxy_FT.pdf")
    


    # --- phase-space plots

    plot_phasespace_classical([alpha, alpha[i_st:]],
                              legend=["full evolution", "steady state"],
                              fname=logging.subdir + "/res_phsp.pdf")



