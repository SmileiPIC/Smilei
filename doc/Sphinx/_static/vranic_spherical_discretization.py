# ______________________________________________________________________________
#
# Geometric study of the Vranic particle merging
#
# ______________________________________________________________________________

import sys
import numpy as np
import math
from matplotlib import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# ______________________________________________________________________________
# RCParams - personalize the figure output

rcParams['figure.facecolor'] = 'w'
rcParams['font.size'] = 15
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['axes.labelsize'] = 15

rcParams['xtick.major.size'] = 10
rcParams['ytick.major.size'] = 10

rcParams['xtick.minor.size'] = 5
rcParams['ytick.minor.size'] = 5

rcParams['axes.linewidth'] = 1.5

rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5

rcParams['axes3d.grid'] = False

# ______________________________________________________________________________
# Functions

def plot_circle_arc(ax,start_point, end_point):

    interval = [0,0,0]
    delta    = [0,0,0]
    angle    = [0,0,0]
    angle_p  = [0,0,0]

    dim_ref = 3
    r = start_point[0]

    interval[1] = abs(end_point[1] - start_point[1])
    interval[2] = abs(end_point[2] - start_point[2])

    delta[1] = interval[1] / dim_ref
    delta[2] = interval[2] / dim_ref

    for i in range(1,dim_ref+1):

        angle_p[1] = (i-1)*delta[1] + start_point[1]
        angle_p[2] = (i-1)*delta[2] + start_point[2]

        angle[1]   = angle_p[1] + delta[1]
        angle[2]   = angle_p[2] + delta[2]

        x_p = r*math.cos(angle_p[1])*math.cos(angle_p[2])
        y_p = r*math.sin(angle_p[1])*math.cos(angle_p[2])
        z_p = r*math.sin(angle_p[2])

        x   = r*math.cos(angle[1])*math.cos(angle[2])
        y   = r*math.sin(angle[1])*math.cos(angle[2])
        z   = r*math.sin(angle[2])

        ax.plot([x, x_p],[y, y_p],[z, z_p],color='k',lw=1)

# ______________________________________________________________________________
# Main

if __name__ == "__main__":

    # __________________________________________________________________________
    # Spherical geometry

    fig = figure(figsize=(8, 8))
    ax = fig.gca(projection='3d')

    r = 1.0

    theta_min = -math.pi * 0.5
    theta_max = math.pi * 0.25
    theta_dim = 20
    theta_delta = abs(theta_max - theta_min) / theta_dim

    phi_min = - math.pi * 0.25
    phi_max = math.pi * 0.5
    phi_dim = 20
    phi_delta = abs(phi_max - phi_min) / phi_dim

    for phi_i in range(theta_dim+1):

        phi_p = (phi_i-1)*phi_delta + phi_min
        phi   = phi_p + phi_delta

        for theta_i in range(1,theta_dim+1):

            theta_p = (theta_i-1)*theta_delta + theta_min
            theta   = theta_p + theta_delta

            x_p = r*math.cos(theta_p)*math.cos(phi)
            x   = r*math.cos(theta)*math.cos(phi)

            y_p = r*math.sin(theta_p)*math.cos(phi)
            y = r*math.sin(theta)*math.cos(phi)

            z_p = r*math.sin(phi)
            z = r*math.sin(phi)

            ax.plot([x, x_p],[y, y_p],[z, z_p],color='k',lw=1)

        if (phi_i > 0):
            for theta_i in range(theta_dim+1):

                theta = theta_i*theta_delta + theta_min

                x_p = r*math.cos(theta)*math.cos(phi_p)
                x   = r*math.cos(theta)*math.cos(phi)

                y_p = r*math.sin(theta)*math.cos(phi_p)
                y = r*math.sin(theta)*math.cos(phi)

                z_p = r*math.sin(phi_p)
                z = r*math.sin(phi)

                #ax.plot([x, x_p],[y, y_p],[z, z_p],color='k',lw=1)
                plot_circle_arc(ax,[r,theta,phi_p], [r,theta,phi])

    ax.set_xlim([-0.8,0.8])
    ax.set_ylim([-0.8,0.8])
    ax.set_zlim([-0.8,0.8])

    ax.set_axis_off()

    fig.tight_layout()

    # __________________________________________________________________________
    # Spherical geometry with solid angle compensation

    fig1 = figure(figsize=(8, 8))
    ax1 = fig1.gca(projection='3d')

    r = 1.0

    theta_min = -math.pi * 0.5
    theta_max = math.pi * 0.25
    theta_interval = abs(theta_max - theta_min)
    theta_dim = 20
    theta_delta_ref = theta_interval / theta_dim

    phi_min = -math.pi * 0.25
    phi_max = math.pi * 0.5
    phi_dim = 20
    phi_delta = abs(phi_max - phi_min) / phi_dim

    for phi_i in range(theta_dim+1):

        phi_p = (phi_i-1)*phi_delta + phi_min
        phi   = phi_p + phi_delta

        theta_delta = theta_delta_ref / math.sin(0.5*math.pi + phi_p + 0.5*phi_delta)

        if (theta_delta > abs(theta_interval)):
            theta_delta = abs(theta_interval)

        theta_dim = int(round(theta_interval / theta_delta))

        theta_delta = theta_interval / theta_dim

        for theta_i in range(1,theta_dim+1):

            theta_p = (theta_i-1)*theta_delta + theta_min
            theta   = theta_p + theta_delta

            x_p = r*math.cos(theta_p)*math.cos(phi)
            x   = r*math.cos(theta)*math.cos(phi)

            y_p = r*math.sin(theta_p)*math.cos(phi)
            y = r*math.sin(theta)*math.cos(phi)

            z_p = r*math.sin(phi)
            z = r*math.sin(phi)

            #ax1.plot([x, x_p],[y, y_p],[z, z_p],color='k',lw=1)
            plot_circle_arc(ax1,[r,theta_p,phi], [r,theta,phi])

        if (phi_i > 0):
            for theta_i in range(theta_dim+1):

                theta = theta_i*theta_delta + theta_min

                x_p = r*math.cos(theta)*math.cos(phi_p)
                x   = r*math.cos(theta)*math.cos(phi)

                y_p = r*math.sin(theta)*math.cos(phi_p)
                y = r*math.sin(theta)*math.cos(phi)

                z_p = r*math.sin(phi_p)
                z = r*math.sin(phi)

                #ax1.plot([x, x_p],[y, y_p],[z, z_p],color='k',lw=1)
                plot_circle_arc(ax1,[r,theta,phi_p], [r,theta,phi])

    ax1.set_xlim([-0.8,0.8])
    ax1.set_ylim([-0.8,0.8])
    ax1.set_zlim([-0.8,0.8])

    ax1.set_axis_off()

    fig1.tight_layout()

    show()
