#
#  ANALYSIS OF TUNNEL IONISATION SIMULATION
#

simulation_to_analyse = '.'

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cmap_map import cmap_map

# IMPORT SMILEI's POST-PROCESSING TOOL
# ------------------------------------

import happi

import subprocess

### RUN SMILEI SIMULATION
subprocess.run(["../../smilei", "input.py"]) 


# LOADING SIMULATION & VARIABLES FROM NAMELIST
# ------------------------------------------------------

S  = happi.Open(simulation_to_analyse ,verbose=False)

t0  = 2.*np.pi
Lv  = 0.
Lp  = S.namelist.Lsim[0]
dt  = S.namelist.Main.timestep
Z_A = S.namelist.Z

a0  = S.namelist.a0
resx = S.namelist.resx
Tsim = S.namelist.Tsim
Lsim = S.namelist.Lsim

print('- ref. ang. frequency w0 = '+str(S.namelist.Main.reference_angular_frequency_SI))

### PLOT FIELD

def plot_field():
    Data = S.Field(0, "Ey").getData()
    i = int(len(Data)/3*1.913)
    x = np.linspace(0, Lsim, Data[i].size)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x/t0, Data[i], linewidth=2.5)

    plt.xlim(x[0]/t0,x[-1]/t0)
    plt.xlabel(r'$x/\lambda$')
    plt.ylabel('Ey [a0]') 


# SIMULATION ANALYSIS & COMPARISON TO RATE EQUATIONS
# --------------------------------------------------


# get corresponding time-steps
t    = dt * np.array(S.ParticleBinning(0).get()['times'])


# Species diagnostic number -> ionization model name
diag_model_dict = {0: 'Tunneling',
                   1: 'Tong-Lin',
                   2: 'Barrier suppression (BSI)',
                   3: 'Full PPT (m != 0)'}


### PLOT CHARGE STATES FOR ONE SPECIES

def plot_charge_states(diag, compare_to=None):
    '''
        diag = 0, 1, 2, or 3 --- diagnostic number
        compare_to = None, 0, 1, 2, or 3 --- if not None, the function plots another diagnostic for comparison
    '''

    # assign a color to each charge state
    cmap = mpl.colormaps.get_cmap('gist_ncar')
    cmap_dark = cmap_map(lambda x: 0.8*x, cmap)
    charge = np.linspace(0, Z_A+1)

    n    = np.array( S.ParticleBinning(diag).getData() )  
    n00  = n[0,0]
    n    = n/n00

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for Z in range(Z_A+1):
        rgba = cmap(Z/float(Z_A+1))
        ax.plot(t/t0, n[:,Z], color=rgba, linewidth=2, label='Z = %d'%Z)
        title = diag_model_dict[diag]

    if compare_to != None:
        n2   = array( S.ParticleBinning(compare_to).getData() )  
        n200  = n2[0,0]
        n2    = n2/n200
        for Z in range(Z_A+1):
            rgba = cmap_dark(Z/float(Z_A+1))
            ax.plot(t/t0, n2[:,Z], '--', color=rgba, linewidth=1, label='Z = %d'%Z)
        title = 'Solid: ' + title + ';    dashed: ' + diag_model_dict[compare_to]
    ax.legend(loc='upper left', ncol=1 if compare_to==None else 2, handletextpad=0.1, labelspacing=0.1)
    ax.set_xlabel(r'$t/T$')
    ax.set_ylabel(r'$n_i/n_{total}$') 
    ax.set_title(title)


def compare_models(tmin=10.5, tmax=14):
    '''
        tmin and tmax set the plot limits alon the x-axis
    '''
    cmap = mpl.colormaps.get_cmap('gist_ncar')
    charge = np.linspace(0, Z_A+1)

    ### plot the density of each charge state as a function of time
    fig, *axes = plt.subplots(nrows=4, ncols=1)
    for i in range(4):
        ax = axes[0][i]
        n    = np.array( S.ParticleBinning(i).getData() )  
        n00  = n[0,0]
        n    = n/n00
        for Z in range(Z_A+1):
            rgba = cmap(Z/float(Z_A+1))
            ax.plot(t/t0, n[:,Z], color=rgba, linewidth=1.2, label='Z = %d'%Z)
            ax.set_xlim([tmin, tmax])
            ax.grid()
            ax.text(tmin+0.05, 
                    0.85, 
                    diag_model_dict[i], 
                    bbox=dict(boxstyle='round', facecolor='w', edgecolor='k', alpha=0.05))
        ax.set_ylabel(r'$n_i/n_{total}$') 
        if i<3:
            ax.set_xticklabels([])
        if i==0:
            ax.legend(loc='upper left', 
                      ncol=Z/3, 
                      handletextpad=0.1, 
                      labelspacing=0.1, 
                      bbox_to_anchor=(0., 1.6))
    axes[0][-1].set_xlabel('t/T')

    plt.subplots_adjust(hspace=0)

    fig.set_size_inches(7, 7)


if  __name__ == '__main__':
    plot_field()

    plot_charge_states(0)
    plot_charge_states(1)
    plot_charge_states(2)
    plot_charge_states(3)

    compare_models(tmin=10.5, tmax=14)

    plt.show(block=True)





