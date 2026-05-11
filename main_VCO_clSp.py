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
from helpers_classical import write_parameters_run
from helpers_classical import plot_phasespace_classical
from helpers_classical import fourier_transform
from helpers_classical import inverse_fourier_transform
from helpers_classical import current_VCO
from helpers_classical import plot_trajectory3d

from scipy.signal import find_peaks

# config-file (as arrays, if we want to switch between ini-files quickly)
config_filenames = ['configfile_VCO_clSp.ini']

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
    times = factory_timegrid(**param_obj.params_now_flat)
    
    # steady state times
    t_st = param_obj.params_now_flat["t_st"]
    i_st = np.where(times >= t_st)[0][0]  

    # zoom times
    t_z0 = param_obj.params_now_flat["t_z0"]
    t_zf = param_obj.params_now_flat["t_zf"]
    i_z0 = np.where(times >= t_z0)[0][0]  
    i_zf = np.where(times >= t_zf)[0][0]  

    # --- circuit
    circuit = factory_circuits(**param_obj.params_now_flat)

    # --- propagator
    propagator = PropagatorClassical(times, circuit, nsteps=1000,
                                     include_noise=circuit.include_noise,
                                     **param_obj.params_now_flat)

    # propagate
    if param_obj.params_now_flat["load_sol"]:

        propagator.res = np.load(logging.subdir + "/res.npy")

    else:

        propagator.propagate()

        if param_obj.params_now_flat["store_sol"]:

            np.save(logging.subdir + "/res", propagator.res)

    print(propagator.res)

    peaks, _ = find_peaks(propagator.res[i_st:, 0])
    A_osc = np.mean((propagator.res[i_st:,0])[peaks])

    # --- Fourier transform of steady state
    res_ft = np.empty_like(propagator.res[i_st:, :])
    for ll in range(len(propagator.res[0])):
        wvec, res_ft[:, ll] = fourier_transform(times[i_st:] - times[i_st], propagator.res[i_st:, ll]) 
 
    # frequencies used for zooming
    w_z0 = param_obj.params_now_flat["w_z0"]
    w_zf = param_obj.params_now_flat["w_zf"]
    i_fz0 = np.where(wvec >= w_z0)[0][0]  
    i_fzf = np.where(wvec >= w_zf)[0][0]  

    # index for frequency w= 0
    i_f0 = np.where(wvec >= 0.0)[0][0]  

    # index for frequency w= w0
    i_fw0 = np.where(wvec >= circuit.w0)[0][0]  

    # main steady-state oscillation frequency
    i_osc = i_f0 + np.argmax(np.abs(res_ft[i_f0:, 0])**2)
    w_osc = wvec[i_osc]


    A_osc0 = 4 * np.sqrt(2/3) * np.sqrt(1 - 1/circuit.alpha_od) 
    w_osc0 = 1 - (circuit.alpha_od-1)**2 / 16 / circuit.Qcoil**2
    print("steady-state oscillation amplitude and frequency: ", A_osc, w_osc)
    print("approx. oscillation amplitude and freq. from paper", A_osc0, w_osc0)

    # plots in Fourier space
    plot_N_in_1(wvec, [np.log10(np.abs(res_ft[:, 0])**2)],
            xlabel='frequency ' + r"$\omega / \omega_0$", ylabel='',
            legend=[r"$log_{10}\left[|y(\omega)|^2\right]$"],
            xlim=[w_z0, w_zf],
            fname=logging.subdir + "/res_w.pdf",
            vline=[circuit.w0]) 

    # plots in time space

    # full result
    plot_N_in_1(times, np.transpose(propagator.res[:, :2]),
                xlabel='time ' + r"$\tau =  t \omega_0$", ylabel='',
                legend=[r"$y(\tau)$", r"$\dot{y}(\tau)$"],
                fname=logging.subdir + "/res_t.pdf",
                vline=[t_st])

    # zoomed result
    plot_N_in_1(times[i_z0:i_zf], np.transpose(propagator.res[i_z0:i_zf, :2]),
                xlabel='time ' + r"$ t \omega_0$", ylabel='',
                legend=[r"$y(\tau)$", r"$\dot{y}(\tau)$"],
                fname=logging.subdir + "/res_t_zoom.pdf")
    
    # result as complex variable in phase-space 
    alpha =  propagator.res[:, 0] + 1.0j*propagator.res[:, 1]
    alpha_rot = alpha  * np.exp(1.0j * w_osc * times)

    # phase-space plot
    plot_phasespace_classical([alpha_rot, alpha_rot[i_st:]],
                              legend=["full evolution", "steady state"],
                              fname=logging.subdir + "/res_phsp.pdf")


    # current of the VCO
    I_d = current_VCO(propagator.res[:, 0], circuit.n)

    # FT of the current 
    wvec, I_d_FT = fourier_transform(times[i_st:] - times[i_st], I_d[i_st:]) 

    # full result
    plot_N_in_1(times[i_z0:i_zf], [I_d[i_z0:i_zf]],
                xlabel='time ' + r"$ t \omega_0$", ylabel='',
                legend=[r"$i_d / I_{bias}(\tau)$"],
                fname=logging.subdir + "/current_zoom.pdf")
    
    # plots in Fourier space
    plot_N_in_1(wvec, [np.log10(np.abs(I_d_FT)**2)],
            xlabel='frequency ' + r"$\omega / \omega_0$", ylabel='',
            legend=[r"$log_{10}\left[|i_d / I_{bias}(\omega)|^2\right]$"],
            xlim=[w_z0, w_zf],
            fname=logging.subdir + "/current_w.pdf",
            vline=[circuit.w0]) 
    
    Mxy = propagator.res[:, 2] + 1.0j * propagator.res[:, 3]
    Mxy_rot = Mxy * np.exp(1.0j * times * w_osc)


    # plot spin motion
    plot_trajectory3d(np.real(Mxy_rot), np.imag(Mxy_rot), propagator.res[:, 4],
                      xlabel=r'$M_x/M_0$', ylabel=r'$M_y/M_0$', zlabel=r'$M_z/M_0$',
                      title="in rotating frame",
                      fname=logging.subdir + "/spin_traj_rot.pdf",
                      elev=15, azim=-75,
                      sphere=True,
                      equal_axes=True)
    
    plot_phasespace_classical([Mxy_rot, Mxy_rot[i_st:]],
                            legend=["full evolution", "steady state"],
                            fname=logging.subdir + "/Mxy_rot.pdf")