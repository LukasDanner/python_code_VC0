import sys

import numpy as np
import os
import qutip
from qutip import *
import numbers


import matplotlib as mpl
#mpl.use('pdf')
import matplotlib.pyplot as plt

from scipy.fftpack import fft, fftfreq, fftshift, ifft


# all these values can be interpreted as 'True'
bool_arr = ['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']

# name of temporary figures
fig_tmp_name = 'Output/Temporary/tmp/fig_tmp'
fig_tmp_count = 1

# name of temporary data files
dat_tmp_name = 'Output/Temporary/tmp/dat_tmp'
dat_tmp_count = 1

def current_VCO(voltage, n):

    return - voltage + voltage**3 / (4 * n)**2

def current_L(voltage, dt, alpha_od, Q):

    F = np.zeros_like(voltage)
    F[1:] = np.cumsum((voltage[1:] + voltage[:-1]) / 2) * dt * Q / 2 / alpha_od

    return F


def lookup_in_dict(dict_, key, default=None):
    '''
    searches a key in a dictionary and returns the respective value
    input:
        - 'dict_'       the dictionary in which we look up the value of the key
        - 'key'         the key of which we want the value
        - 'default'     a given default value which should be given if the key is not found in dict_

    returns the value of the key
    '''

    try:
        # try to find the key in dictionary and set the value of the key
        val = dict_[key]

    except KeyError:

        # if key was not found, set val to the given default value
        val = default

    # return the value of the key
    return val

def str_to_array(str_, separator=','):
    '''
    input:
        - 'str_'        a string of characters which are separated by one specific character 'separator'
        - 'separator'   a character which separates the string 'str_'

    returns an array, where we split the string after each separator
    '''

    return str_.replace(' ', '').split(separator)


def change_datatype_array_elements(arr, new_type):
    '''
    changes the data-type of all array-elements

    input:
        - 'arr'             the array
        - 'new_type'        the new data type

    output:
    the same array, but all element are typecasted
    '''

    global datatypes_global
    global bool_arr

    # special case where the new data type is a numpy array itself
    # then: we expect that 'arr' contains single elements which contain opening '[' and closing ']'brakets
    if new_type == 'array':

        # count how many arrays we have
        count =0

        for ll in range(len(arr)):

            # count the number of opening brackets:
            if (arr[ll])[0] == '[':
                count +=1

        # size of a single array
        size = int(len(arr)/count)

        # allocate an array of arrays, where each array has given size
        ret_arr = [np.zeros(size) for kk in range(count)]

        # go through the given array
        for ll in range(len(arr)):

            # set the value in the ret_arr
            # eventually get rid of brackets before typecast
            ret_arr[int(ll // size)][ll%size] = np.double(arr[ll].strip('[]'))

        return ret_arr

    # normal cases: simple type cast sufficient
    # go through all array elements
    for ii in range(len(arr)):

        # type cast
        if new_type == 'float':

            arr[ii] = float(arr[ii])

        if new_type == 'np.double':

            arr[ii] = np.double(arr[ii])

        if new_type == 'int':

            arr[ii] = int(arr[ii])

        if new_type == 'bool':

            if arr[ii].lower() in bool_arr:

                arr[ii] = True

            else:

                arr[ii] = False

    return arr

def max_length_array_in_dict(dict_):
    '''
    input
        - 'dict_'       a dictionary with keys and values
                        the values at the keys are arrays with different lengths

    this function finds and returns the maximum length of the arrays
    '''

    # initialize max. length
    N_max = 0

    # go through all keys
    for key in dict_:

        # adjust the max. length if the value is an array with larger length
        if len(dict_[key]) > N_max:
            N_max = len(dict_[key])

    return N_max

def arr_to_arr_of_given_length(arr_, N):
    '''
    input:
        - 'arr_'        an array with some length
        - 'N'           the wanted length of the array

    output:
        - an array with length N with the property
            - if input array has length N, then just the input array
            - if input array does not have length N, then an array [arr_[0],arr_[0],...,arr_[0]]
    '''

    # if arr_ fits to wanted length, simply return the arr_ again
    if len(arr_) == N:

        return arr_

    # if arr_ has smaller length
    # then just return an array with length N, where all values are set to arr_[0]
    else:

        ret_arr = [arr_[0] for ll in range(N)]
        return ret_arr



def flatten_nested_dicts(nested_dict, arr_index):
    '''
    input:
        - 'nested_dict'
        the nested dictionary, where...
            the keys are section names
            the values at the keys are itself dictionaries, where...
                the keys are parameter names
                the values are arrays with parameter values for all loop-runs
        - 'arr_index'
          the index of the loop-run

    flattens the outer dictionary of a nested dictionary
    for the parameter values, we pick the value in the array at index arr_index

    returns a flattened dictionary
    '''

    # allocate a new, empty dictionary
    dict_ = dict()

    # go through all section names (=keys of the outer dictionary)
    for section_name in nested_dict:

        # go through all keys of the inner dictionary (= parameter names of the ini-file)
        for key in nested_dict[section_name]:

            # add the key to the dict
            # add the value of the array at index arr_index as value for the key
            dict_[key] = nested_dict[section_name][key][arr_index]

    # return the flat dictionary
    return dict_


def search_nested_dict(nested_dict, wanted_key):
    '''
    input:
        - 'nested_dict'
        the nested dictionary, where...
            the keys are section names
            the values at the keys are itself dictionaries, where...
                the keys are parameter names
                the values are arrays with parameter values for all loop-runs
        - 'wanted_key'
          the wanted key

    search in the inner dictionary for the key
    return a dictionary where:
        key   = wanted_key
        value = all values at the key
    '''

    # go through outer-dictionary keys
    for section_name in nested_dict:

        # go through inner-dictionary keys
        for key in nested_dict[section_name]:

            # if this key is the wanted key
            if key == wanted_key:

                # return a dictionary consisting of this key and its values
                return {key: nested_dict[section_name][key]}


def write_parameters_run(dict, logging, fname="config_exact_params.txt"):
    '''


    '''

    # open a new file in the current ouput-subfolder
    with open(logging.subdir+"/"+fname, 'w') as f:

        # loop through outer dict
        for key, value in dict.items():

            # write the key (=section name)
            f.write('[%s]\n' % key)

            # loop through inner dict
            for key2, value2 in value.items():

                # write parameter name (key2) and its value (value2) to the file
                f.write('%s\t\t%s\n' % (key2, value2))

def log_to_file(data,
                fname=None,     # file name
                header='',      # name of header in file
                comments='#',   # for comments (any line with commands starts with this character)
                newline='\n',   # character for new line
                append=True     # if file shall be appended or overwritten
                ):
    '''
    stores data into a (text-)file
    '''

    # catch error if no file name is specified
    if fname is None:

        global dat_tmp_count
        fname = dat_tmp_name + str(dat_tmp_count) + '.txt'
        dat_tmp_count += 1

    # if file shall be appended
    if append and os.path.exists(fname):

        with open(fname, 'a') as file:

            # do not write header twice
            np.savetxt(file, data, newline=newline, comments=comments)

    else:

        np.savetxt(fname, data, header=header, newline=newline, comments=comments)



def plot_ReIm(x, y,
              xlabel='', ylabel='', title=None,
              fname=None, hline=None, xlim=None, ylim=None):
    '''
    quick plot of the real part and imaginary part of a complex function

    input:
        x           array specifying the x-axis
        y           array specifying the y-axis (can be complex)
        xlabel      label for x-axis
        ylabel      label for y-axis
        fname       filename as full path
        hline       array, specifying y-positions of horizontal lines
                    eg: [1, 2] will mean that two horizontal lines are plotted at y=1 and y=2
        xlim        array, specifying the lower and upper bound for x-axis
        ylim        array, specifying the lower and upper bound for y-axis
    '''

    # two subplots:

    fig, axs = plt.subplots(2)

    # upper plot is real part
    axs[0].plot(x, np.real(y))

    # bottom plot is imag part
    axs[1].plot(x, np.imag(y))

    # set axis-labels
    axs[1].set(xlabel=xlabel)
    axs[0].set(ylabel='Re[' + ylabel + ']')
    axs[1].set(ylabel='Im[' + ylabel + ']')

    # print title
    if title is not None:
        axs[0].set(title=title)

    # if given, draw horizontal lines
    if hline is not None:

        # go through all lines given
        for kk in range(len(hline)):

            # plot this horizontal line
            plt.hlines(hline[kk], x[0], x[-1], colors="black", linestyles=":")

    # set limits for x-axis
    if xlim is not None:
        axs[0].set_xlim(xlim)
        axs[1].set_xlim(xlim)

    # set limits for y-axis
    if ylim is not None:
        axs[0].set_ylim(ylim)
        axs[1].set_ylim(ylim)

    # handle case where no file name is given
    if fname is None:

        global fig_tmp_count
        plt.savefig(fig_tmp_name + str(fig_tmp_count) + '.pdf')
        fig_tmp_count += 1

    else:

        # save the figure
        plt.savefig(fname)

    # close the plot
    plt.close()


def plot_N_in_1(x, ymat,
                xlabel='', ylabel='',
                legend=None, xlim=None, ylim=None, linestyle=None, fname=None,
                hline=None, vline=None, title=None, yscale=None):
    '''
    quick plot of the real part and imaginary part of a complex function
    input:
        x           array specifying the x-axis
        ymat        2d-array specifying the functions to plot (must be real)
                    ymat[kk] are the functions to be plotted and
                    shall be an array with the length of x
        xlabel      label for x-axis
        ylabel      label for y-axis
        legend      array, with labels for different functions
        xlim        array, specifying the lower and upper bound for x-axis
        ylim        array, specifying the lower and upper bound for y-axis
        linestyle   array, specifying matplotlib-linestyles for different functions
        fname       filename as full path
        hline       array, specifying y-positions of horizontal lines
                    eg: [1, 2] will mean that two horizontal lines are plotted at y=1 and y=2
        vline       array, specifying x-positions of vertical lines
                    eg: [1, 2] will mean that two vertical lines are plotted at x=1 and x=2
        yscale      scale for y-axis, if we want to plot logarithmically
    '''

    # open the plot
    fig, ax = plt.subplots()

    # print title
    if title is not None:
        plt.title(title)

    # if no legends are given, make None-array
    if legend is None:
        legend = [None for kk in range(len(ymat))]

    # if no linestyles are given, then all lines will be normal
    if linestyle is None:
        linestyle = ['-' for kk in range(len(ymat))]

    # set y-limits
    if ylim is not None:
        ax.set_ylim(ylim)

    # set x-limits
    if xlim is not None:
        ax.set_xlim(xlim)

    # set scale for y-axis
    if yscale is not None:
        ax.set_yscale(yscale)

    # set horizontal lines
    # see explanations in plot_ReIm
    if hline is not None:

        for kk in range(len(hline)):
            plt.hlines(hline[kk], x[0], x[-1], colors="black", linestyles=":")


    # set vertical lines
    # see explanations in plot_ReIm
    if vline is not None:

        for kk in range(len(vline)):
            plt.axvline(vline[kk], color="black", linestyle=":")

    # plot all functions
    for kk in range(len(ymat)):

        # ymat[kk] must be an array with same length as the array x
        plt.plot(x, ymat[kk], label=legend[kk], linestyle=linestyle[kk])

    # print legend at best location available
    plt.legend(loc='best')

    # set xlabel
    ax.set(xlabel=xlabel)

    # set ylabel
    ax.set(ylabel=ylabel)

    # store the result
    if fname is None:

        global fig_tmp_count
        plt.savefig(fig_tmp_name + str(fig_tmp_count) + '.pdf')
        fig_tmp_count += 1

    else:

        plt.savefig(fname)

    # close the plot
    plt.close()


def plot_N_in_1_diff_xaxis(xmat, ymat,
                           xlabel='', ylabel='',
                           legend=None, xlim=None, ylim=None, linestyle=None, fname=None,
                           hline=None, vline=None, title=None, yscale=None):
    '''
    quick plot of the real part and imaginary part of a complex function
    input:
        x           array specifying the x-axis
        ymat        2d-array specifying the functions to plot (must be real)
                    ymat[kk] are the functions to be plotted and
                    shall be an array with the length of x
        xlabel      label for x-axis
        ylabel      label for y-axis
        legend      array, with labels for different functions
        xlim        array, specifying the lower and upper bound for x-axis
        ylim        array, specifying the lower and upper bound for y-axis
        linestyle   array, specifying matplotlib-linestyles for different functions
        fname       filename as full path
        hline       array, specifying y-positions of horizontal lines
                    eg: [1, 2] will mean that two horizontal lines are plotted at y=1 and y=2
        vline       array, specifying x-positions of vertical lines
                    eg: [1, 2] will mean that two vertical lines are plotted at x=1 and x=2
        yscale      scale for y-axis, if we want to plot logarithmically
    '''

    # open the plot
    fig, ax = plt.subplots()

    # print title
    if title is not None:
        plt.title(title)

    # if no legends are given, make None-array
    if legend is None:
        legend = [None for kk in range(len(ymat))]

    # if no linestyles are given, then all lines will be normal
    if linestyle is None:
        linestyle = ['-' for kk in range(len(ymat))]

    # set y-limits
    if ylim is not None:
        ax.set_ylim(ylim)

    # set x-limits
    if xlim is not None:
        ax.set_xlim(xlim)

    # set scale for y-axis
    if yscale is not None:
        ax.set_yscale(yscale)

    # set horizontal lines
    # see explanations in plot_ReIm
    if hline is not None:

        for kk in range(len(hline)):
            plt.hlines(hline[kk], x[0], x[-1], colors="black", linestyles=":")


    # set vertical lines
    # see explanations in plot_ReIm
    if vline is not None:

        for kk in range(len(vline)):
            plt.axvline(vline[kk], color="black", linestyle=":")

    # plot all functions
    for kk in range(len(ymat)):

        # ymat[kk] must be an array with same length as the array x
        plt.plot(xmat[kk], ymat[kk], label=legend[kk], linestyle=linestyle[kk])

    # print legend at best location available
    plt.legend(loc='best')

    # set xlabel
    ax.set(xlabel=xlabel)

    # set ylabel
    ax.set(ylabel=ylabel)

    # store the result
    if fname is None:

        global fig_tmp_count
        plt.savefig(fig_tmp_name + str(fig_tmp_count) + '.pdf')
        fig_tmp_count += 1

    else:

        plt.savefig(fname)

    # close the plot
    plt.close()


def plot_phasespace_classical(alpha, xlim=None, ylim=None, 
                              xlabel="position " + r"$ y$" , ylabel="momentum " + r"$\dot{y}$", 
                              title="",
                              legend=None, 
                              circle=None,
                              fname=None,
                              points=None,
                              ):


    N_traj = len(alpha)

    # open the plot
    fig, ax = plt.subplots()

    if xlim is None:

        xmax = np.max(np.abs(np.real(np.concatenate(alpha))))
        xlim = [-xmax, xmax]

    if ylim is None:

        pmax = np.max(np.abs(np.imag(np.concatenate(alpha))))
        ylim = [-pmax, pmax]

    # if no legends are given, make None-array
    if legend is None:
        legend = [None for kk in range(N_traj)]


    # aspect ratio of plot
    # set aspect ratio of plot correctly
    xlength = xlim[1] - xlim[0]
    plength = ylim[1] - ylim[0]
    ax.set_box_aspect(plength / xlength)

    # set title of the plot
    plt.title(title)

    # set xlabel
    ax.set(xlabel=xlabel)

    # set ylabel
    ax.set(ylabel=ylabel)

    ax.set(xlim=xlim)

    ax.set(ylim=ylim)

    if circle is not None:

        vec = circle * np.exp(1.0j * np.linspace(0, 2 * np.pi, 201))
        plt.plot(np.real(vec), np.imag(vec), '--', linewidth=0.04, color="black")

    if points is not None:

        for ll in range(len(points)):
            plt.plot(np.real(points[ll]), np.imag(points[ll]), 'o')
        plt.gca().set_prop_cycle(None)

    # actual plot
    for ll in range(N_traj):

        plt.plot(np.real(alpha[ll]), np.imag(alpha[ll]), '--',  
                linewidth=0.02, 
                # markersize=0.2, marker=".", 
                label=legend[ll])
        

    plt.legend(loc='best')
    
    # save the figure
    if fname is None:

        global fig_tmp_count
        plt.savefig(fig_tmp_name + str(fig_tmp_count) + '.pdf')
        fig_tmp_count += 1

    else:

        plt.savefig(fname)

    # close the plot
    plt.close()

def plot_trajectory3d(x, y, z,
                      xlabel='Sx', ylabel='Sy', zlabel='Sz',
                      title=None,
                      fname=None,
                      elev=25, azim=45,
                      sphere=True,
                      equal_axes=True):
    '''
    quick 3D plot of a trajectory (eg: spin precession)

    input:
        x           array specifying x-component versus time
        y           array specifying y-component versus time
        z           array specifying z-component versus time

        xlabel      label for x-axis
        ylabel      label for y-axis
        zlabel      label for z-axis

        title       plot title

        fname       filename as full path

        elev        elevation angle for 3D view
        azim        azimuth angle for 3D view

        sphere      bool, if True plot sphere with radius
                    equal to magnitude of initial vector

        equal_axes  bool, if True use equal scaling on all axes
    '''

    import numpy as np
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import Axes3D

    # create figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot trajectory
    ax.plot(x, y, z, linewidth=2)

    # mark initial point
    ax.scatter(x[0], y[0], z[0],
               color='green',
               s=60,
               label='initial')

    # mark final point
    ax.scatter(x[-1], y[-1], z[-1],
               color='red',
               s=60,
               label='final')

    # sphere with radius of initial vector
    if sphere:

        # radius from initial values
        r = np.sqrt(x[0]**2 + y[0]**2 + z[0]**2)

        # spherical coordinates
        u = np.linspace(0, 2*np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # sphere surface
        xs = r * np.outer(np.cos(u), np.sin(v))
        ys = r * np.outer(np.sin(u), np.sin(v))
        zs = r * np.outer(np.ones(np.size(u)), np.cos(v))

        # plot sphere
        ax.plot_surface(xs, ys, zs,
                        color='lightblue',
                        alpha=0.15,
                        linewidth=0)

    # determine symmetric axis limits around origin
    maxval = np.max(np.abs(np.concatenate([x, y, z])))

    ax.set_xlim([-maxval, maxval])
    ax.set_ylim([-maxval, maxval])
    ax.set_zlim([-maxval, maxval])

    # draw coordinate axes through origin
    ax.plot([-maxval, maxval], [0, 0], [0, 0],
            color='black', linewidth=1)

    ax.plot([0, 0], [-maxval, maxval], [0, 0],
            color='black', linewidth=1)

    ax.plot([0, 0], [0, 0], [-maxval, maxval],
            color='black', linewidth=1)
    
    # set labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    # set title
    if title is not None:
        ax.set_title(title)

    # equal axis scaling
    if equal_axes:

        max_range = np.array([
            x.max() - x.min(),
            y.max() - y.min(),
            z.max() - z.min()
        ]).max() / 2.0

        mid_x = (x.max() + x.min()) * 0.5
        mid_y = (y.max() + y.min()) * 0.5
        mid_z = (z.max() + z.min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # set viewing angle
    ax.view_init(elev=elev, azim=azim)

    # legend
    ax.legend()

    # handle case where no file name is given
    if fname is None:

        global fig_tmp_count
        plt.savefig(fig_tmp_name + str(fig_tmp_count) + '.pdf')
        fig_tmp_count += 1

    else:

        # save figure
        plt.savefig(fname)

    # close plot
    plt.close()

def fourier_transform(tvec, func, norm="forward"):

    # time step
    dt = np.abs(tvec[1] - tvec[0])

    # number of discrete time steps
    Ntimes = len(tvec)# -1



    # frequency vector
    w_vec = fftfreq(Ntimes, dt)

    # fourier transform
    res_ft = np.fft.fft(func, norm=norm)  #/ len(tvec)

    # re-order frequencies from most negative to most positive (centre on 0)
    idx = np.array([], dtype='int')
    idx = np.append(idx, np.where(w_vec < 0.0))
    idx = np.append(idx, np.where(w_vec >= 0.0))

    w_vec = 2 * np.pi * w_vec[idx]
    res_ft = res_ft[idx]

    #print(w_vec)

    return w_vec, res_ft


def inverse_fourier_transform(wvec, ft, norm="forward"):

    dt = 2*np.pi/(wvec[-1] - wvec[0])

    tf = 2*np.pi/(wvec[1] - wvec[0])

    #print(tf)

    # this is changed now, I don't know if true
    #tvec = np.linspace(0, tf-dt, len(wvec))
    #tvec = np.linspace(0, tf, len(wvec))

    tvec = np.linspace(-tf/2, tf/2, len(wvec))


    ft = np.fft.ifftshift(ft)
    func = np.fft.ifft(ft, norm=norm) #* np.sqrt(len(tvec))
    func = fftshift(func)

    return tvec, func


def wvec_scale_to_gamma(wvec, w0, Qfac):

    return Qfac * (wvec - w0) / w0

def find_peak_indices(w, f, thres=0.0):

    res = []

    for ll in range(1, len(w)-1):

        if np.abs(f[ll]) > np.abs(f[ll-1]) + thres and np.abs(f[ll]) > np.abs(f[ll+1]) + thres:

            res += [ll]

    # sort the array
    res = np.array(res)
    idx = np.argsort(f[res])[::-1]

    return res[idx]


def find_positive_peak_indices(w, f, thres=0.0):

    res = []

    for ll in range(1, len(w)-1):

        if np.abs(f[ll]) > np.abs(f[ll-1]) + thres and np.abs(f[ll]) > np.abs(f[ll+1]) + thres and w[ll] >= 0.0:

            res += [ll]

    # sort the array
    res = np.array(res)
    idx = np.argsort(f[res])[::-1]

    print(res, f[res], idx)

    return res[idx]

def find_max_index(arr):

    return np.argmax(arr)

def find_avg_oscillation_freq(w, f):

    i0 = np.where(w >= 0.5)[0][0]

    norm = sum(np.abs(f[i0:])**2)
    mean = sum(w[i0:] * np.abs(f[i0:])**2)
    #print(mean, norm)

    return mean/norm


def find_avg_slope(t, f):

    slope, boff = np.polyfit(t, f, 1)

    return slope, boff