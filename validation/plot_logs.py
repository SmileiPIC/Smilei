from numpy import *
from matplotlib.pyplot import *
from matplotlib.widgets import Button
from glob import glob
from json import load
from time import strptime, mktime

cases = sorted(glob("/sps2/gitlab-runner/logs/*.log"))

# Class to create a plot with buttons to switch between cases
class switchPlots:
    ind = 0
    
    def __init__(self, plotter, nplots):
        self.plotter = plotter
        self.nplots = nplots
        figure(1000)
        clf()
        self.ax = axes()
        axprev = axes([0.7, 0.01, 0.1, 0.03])
        axnext = axes([0.81, 0.01, 0.1, 0.03])
        self.bnext = Button(axnext, '>')
        self.bnext.on_clicked(self.next)
        self.bprev = Button(axprev, '<')
        self.bprev.on_clicked(self.prev)
        self.plotter(self.ax, self.ind)
    
    def next(self, event):
        self.ind += 1
        self.ind %= self.nplots
        self.plotter(self.ax, self.ind)
    
    def prev(self, event):
        self.ind -= 1
        self.ind %= self.nplots
        self.plotter(self.ax, self.ind)

# Plotter for validation logs
def myplotter(ax, i):
    # Read data
    with open(cases[i], 'r') as f:
        ax.data = load(f)
    
    # Get time
    t = array([mktime(strptime(d, "%Y_%m_%d_%H:%M:%S")) for d in ax.data["date"]], dtype=int)
    t -= t[0]
    
    # Get all timers
    y = []
    labels = []
    for k in ['Particles','Maxwell','Densities',
              'SyncParticles','SyncFields','SyncDensities','SyncSusceptibility',
              'Diagnostics','Envelope','Collisions','Movwindow',
              'Loadbalancing','Reconfiguration','PartMerging','PartInjection']:
        if k in ax.data:
            d = ax.data[k] + [None]*(t.size-len(ax.data[k]))
            d = array(d, dtype="float")
            d[isnan(d)] = 0.
            y += [d]
            labels += [k]
    
    # Get branch name to detect those which are not "develop"
    commits = ax.data["commit"]
    branches = array(["-".join(b.split("-")[1:]) for b in commits])
    branch_points = {}
    for b in unique(branches):
	if "HEAD" in b or b=="develop": continue
        ind = flatnonzero(branches == b)
        branch_points[b] = [t[ind], array(ax.data["Timeintimeloop"],dtype="float")[ind]] 
    
    # Plot
    sca(ax)
    ax.cla()
    title(cases[i])
    stackplot( t, *y, labels=labels )
    plot( t, array(ax.data["Timeintimeloop"]), '--k', label="time loop (total)")
    for i,b in enumerate(branch_points):
        plot( branch_points[b][0], branch_points[b][1], linestyle="none", marker=i+4, label=b )
    legend()
    show()

# Launch the plotter in a switchPlots window
plots = switchPlots( myplotter, len(cases) )

