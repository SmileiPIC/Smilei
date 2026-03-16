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

def plot_3d_rectangle(origin,v1,v2,alpha = 1.,color = 'red',factor=1):

    x = np.array([[origin[0],origin[0]+factor*v1[0]],[origin[0]+factor*v2[0],origin[0]+factor*(v1[0]+v2[0])]])
    y = np.array([[origin[1],origin[1]+factor*v1[1]],[origin[1]+factor*v2[1],origin[1]+factor*(v1[1]+v2[1])]])
    z = np.array([[origin[2],origin[2]+factor*v1[2]],[origin[2]+factor*v2[2],origin[2]+factor*(v1[2]+v2[2])]])

    ax.plot_surface(x, y, z, alpha=alpha,color=color)

def perpendicular_angle(v1,v2,color='k',factor=1):

    # Perpendicular angle
    x = np.array([v1[0],v1[0]+v2[0]]) * factor
    y = np.array([v1[1],v1[1]+v2[1]]) * factor
    z = np.array([v1[2],v1[2]+v2[2]]) * factor
    ax.plot(x,y,z,color=color)
    x = np.array([v1[0]+v2[0],v2[0]]) * factor
    y = np.array([v1[1]+v2[1],v2[1]]) * factor
    z = np.array([v1[2]+v2[2],v2[2]]) * factor
    ax.plot(x,y,z,color=color)

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

if __name__ == "__main__":

    # __________________________________________________________________________
    # Parameters

    xmin = -2
    xmax = 2

    ymin = -2
    ymax = 2

    zmin = -2
    zmax = 2


    thetat = np.pi * 0.1
    phit = np.pi * 0.1

    pt = 1.5*np.array([np.cos(phit)*np.cos(thetat), np.cos(phit)*np.sin(thetat), np.sin(phit)])

    thetad = np.pi * 0.15
    phid = np.pi * 0.15

    d = np.array([np.cos(phid)*np.cos(thetad), np.cos(phid)*np.sin(thetad), np.sin(phid)])
    wt = 1.

    pt_n = np.linalg.norm(pt)
    et   = np.sqrt(pt_n**2+1.)

    print("pt: {}".format(pt))
    print("norm(pt): {}".format(pt_n))
    print("et: {}".format(et))

    # __________________________________________________________________________
    # Orthonorme Basis

    e1 = pt / pt_n

    #e3 = -np.cross(e1,d) / np.linalg.norm(d)

    e3 = np.zeros(3)
    e3[0] = -(e1[1]*d[2] - e1[2]*d[1]) / np.linalg.norm(d)
    e3[1] = -(e1[2]*d[0] - e1[0]*d[2]) / np.linalg.norm(d)
    e3[2] = -(e1[0]*d[1] - e1[1]*d[0]) / np.linalg.norm(d)
    e3 = e3 / np.linalg.norm(e3)

    #e2 = np.cross(e1,e3)

    e2 = np.zeros(3)
    e2[0] = e1[1]**2*d[0] - e1[0]*(d[1]*e1[1] + d[2]*e1[2]) + e1[2]**2*d[0]
    e2[1] = e1[2]**2*d[1] - e1[1]*(d[2]*e1[2] + d[0]*e1[0]) + e1[0]**2*d[1]
    e2[2] = e1[0]**2*d[2] - e1[2]*(d[0]*e1[0] + d[1]*e1[1]) + e1[1]**2*d[2]
    e2 = e2 / np.linalg.norm(e2)

    print("e1: {}".format(e1))
    print("e2: {}".format(e2))
    print("e3: {}".format(e3))
    print("norm(e1): {}".format(np.linalg.norm(e1)))
    print("norm(e2): {}".format(np.linalg.norm(e2)))
    print("norm(e3): {}".format(np.linalg.norm(e3)))
    print("norm(d): {}".format(np.linalg.norm(d)))
    print("e1 . e3 = {}".format(np.dot(e1,e3)))
    print("e1 . e2 = {}".format(np.dot(e1,e2)))
    print("e2 . e3 = {}".format(np.dot(e2,e3)))

    # __________________________________________________________________________
    # Computation of pa and pb

    ea   = et / wt
    eb   = et / wt

    pa_n = (ea**2 - 1)
    pb_n = (eb**2 - 1)

    print("ea: {}".format(ea))
    print("pa_n: {}".format(pa_n))
    print("pt_n / (wt * pa_n) = {}".format(pt_n / (wt * pa_n)))

    theta = math.acos(pt_n / (wt * pa_n))

    print("theta= {}".format(theta))

    pa = pa_n * (np.cos(theta) * e1 + np.sin(theta) * e2)
    pb = pb_n * (np.cos(-theta) * e1 + np.sin(-theta) * e2)

    print("pa . e3 = {}".format(np.dot(pa,e3)))
    print("pb . e3 = {}".format(np.dot(pb,e3)))
    print("norm(pa) = {}, norm(pb) = {}".format(np.linalg.norm(pa), np.linalg.norm(pb)))

    # __________________________________________________________________________
    # Plane

    amplitude = 3.

    e1 = e1*amplitude
    e2 = e2*amplitude
    e3 = e3*amplitude

    # __________________________________________________________________________
    # Figure

    fig = figure(figsize=(8, 8))
    ax = fig.gca(projection='3d')

    # x = np.array([[e1[0],e2[0]],[-e2[0],-e1[0]]])
    # y = np.array([[e1[1],e2[1]],[-e2[1],-e1[1]]])
    # z = np.array([[e1[2],e2[2]],[-e2[2],-e1[2]]])
    # ax.plot_surface(x, y, z, alpha=0.2,color='green')

    # Perpendicular angle
    factor = 1./15.
    perpendicular_angle(e1,e2,color='r',factor=factor)
    perpendicular_angle(e1,e3,color='r',factor=factor)
    perpendicular_angle(e2,e3,color='r',factor=factor)

    plot_3d_rectangle([0,0,0],e1 ,e2 ,0.5,factor=factor)
    plot_3d_rectangle([0,0,0],e1 ,e3 ,0.5,factor=factor)
    plot_3d_rectangle([0,0,0],e2 ,e3 ,0.5,factor=factor)

    d_arrow = Arrow3D([0, d[0]], [0, d[1]],
                [0, d[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="purple")
    ax.add_artist(d_arrow)
    ax.text(d[0], d[1], d[2], r'$d$',None,color="purple")

    e1_arrow = Arrow3D([0, e1[0]], [0, e1[1]],
                [0, e1[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="r")
    ax.add_artist(e1_arrow)
    ax.text(e1[0], e1[1], e1[2], r'$e_1$',None,color="r")

    e2_arrow = Arrow3D([0, e2[0]], [0, e2[1]],
                [0, e2[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="r")
    ax.add_artist(e2_arrow)
    ax.text(e2[0], e2[1], e2[2], r'$e_2$',None,color="r")

    e3_arrow = Arrow3D([0, e3[0]], [0, e3[1]],
                [0, e3[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="r")
    ax.add_artist(e3_arrow)
    ax.text(e3[0], e3[1], e3[2], r'$e_3$',None,color="r")

    pa_arrow = Arrow3D([0, pa[0]], [0, pa[1]],
                [0, pa[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="b")
    ax.add_artist(pa_arrow)
    ax.text(pa[0], pa[1], pa[2], r'$p_a$',None,color="b")

    pb_arrow = Arrow3D([0, pb[0]], [0, pb[1]],
                [0, pb[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="b")
    ax.add_artist(pb_arrow)
    ax.text(pb[0], pb[1], pb[2], r'$p_b$',None,color="b")

    x_arrow = Arrow3D([0, 1], [0, 0],
                [0, 0], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="grey")
    ax.add_artist(x_arrow)
    ax.text(1, 0, 0, r'$x$',None,color="grey")

    z_arrow = Arrow3D([0, 0], [0, 0],
                [0, 1], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="grey")
    ax.add_artist(z_arrow)
    ax.text(0, 0, 1, r'$z$',None,color="grey")

    y_arrow = Arrow3D([0, 0], [0, 1],
                [0, 0], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="grey")
    ax.add_artist(y_arrow)
    ax.text(0, 1, 0, r'$y$',None,color="grey")

    pt_arrow = Arrow3D([0, pt[0]], [0, pt[1]],
                [0, pt[2]], mutation_scale=20,
                lw=1, arrowstyle="-|>", color="k")
    ax.text(pt[0], pt[1], pt[2], r'$p_t/w_t$',None,color="k")
    ax.add_artist(pt_arrow)

    ax.plot([pt[0],pa[0]],[pt[1],pa[1]],[pt[2],pa[2]],ls='--',color='k')

    e2 = e2 / np.linalg.norm(e2)
    ax.plot([pa[0],pa_n*e2[0]*math.sin(theta)],[pa[1],pa_n*e2[1]*math.sin(theta)],[pa[2],pa_n*e2[2]*math.sin(theta)],ls='--',color='k')

    # __________________________________________________________________________
    # Plot

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_zlim([zmin,zmax])

    ax.plot([2*xmin ,2*xmax],[0,0],[0,0],color='grey',lw=1)
    ax.plot([0 ,0],[2*ymin,2*ymax],[0,0],color='grey',lw=1)
    ax.plot([0 ,0],[0,0],[2*zmin,2*zmax],color='grey',lw=1)

    ax.set_axis_off()

    fig.tight_layout()

    show()
