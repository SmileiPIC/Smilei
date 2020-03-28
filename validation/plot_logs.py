from numpy import *
from matplotlib.pyplot import *
from matplotlib.widgets import Button
from matplotlib import dates
from glob import glob
from json import load
from time import strptime
from datetime import datetime
ion()

# RCParams parameters to tune the matplotlib style

rcParams['figure.subplot.top'] = 0.94
rcParams['figure.subplot.bottom'] = 0.05

rcParams['figure.subplot.right'] = 0.9
rcParams['figure.subplot.left'] = 0.06

rcParams['axes.linewidth'] = 1.5

rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5

# Custom list

custom_colors = ['C0','C1','C2','C3','C4','C5']
custom_markers = ['o','s','^','v','<','>']

# Obtain the data

cases = sorted(glob("/sps2/gitlab-runner/logs/*.log"))
#cases = sorted(glob("./logs/*.log"))

# Class to create a plot with buttons to switch between cases
class switchPlots:
    
    # Class parameters
    
    ind               = 0
    branch_points     = {}
    marker_properties = {}
    branches          = array([])
    marker_size       = 8
    menu_x_position   = 0.75
    
    def __init__(self, nplots):
        
        self.nplots = nplots
        self.fig = figure(figsize=(14, 6))
        clf()
        gs = GridSpec(10,10)
        self.ax = subplot(gs[0:9,0:8])
        self.pick_event = self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        
        # Button next
        axnext = axes([0.82, 0.01, 0.05, 0.03])
        self.bnext = Button(axnext, '>')
        self.bnext.on_clicked(lambda event : self.next(event,1))

        # Button previous
        axprev = axes([0.76, 0.01, 0.05, 0.03])
        self.bprev = Button(axprev, '<')
        self.bprev.on_clicked(lambda event : self.next(event,-1))
   
        # Button next +5
        axnext5 = axes([0.88, 0.01, 0.05, 0.03])
        self.bnext5 = Button(axnext5, '>>')
        self.bnext5.on_clicked(lambda event : self.next(event,5))

        # Button previous -5
        axprev5 = axes([0.70, 0.01, 0.05, 0.03])
        self.bprev5 = Button(axprev5, '<<')
        self.bprev5.on_clicked(lambda event : self.next(event,-5))
        
        #self.plotter(self.ax, self.fig, self.ind)
        self.plot()
    
    def next(self, event, jump):
        self.ind += jump
        self.ind %= self.nplots
        #self.plotter(self.ax, self.fig, self.ind)
        self.plot()

    def on_pick(self,event):
        """
        This method is used to make markers clickable.
        Information is then shown in the terminal.
        """
        artist = event.artist
        xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
        #x, y = artist.get_xdata(), artist.get_ydata()
        marker_indexes = event.ind
        print("Selected points: {}".format(marker_indexes))
        for k in marker_indexes:
            print("  > Commit {}".format(self.ax.data["commit"][k]))
            print("    Branch {}".format(self.branches[k]))
            print("    Date {}".format(self.ax.data["date"][k]))

    def plot(self):
        """
        Method to plot the log of index self.ind.
        """
        # Read data
        with open(cases[self.ind], 'r') as f:
            self.ax.data = load(f)
        
        # Get time
        t = array([dates.date2num(datetime(*strptime(d, "%Y_%m_%d_%H:%M:%S")[:6])) for d in self.ax.data["date"]])
        
        # Get all timers
        y = []
        labels = []
        for k in ['Particles','Maxwell','Densities',
                  'SyncParticles','SyncFields','SyncDensities','SyncSusceptibility',
                  'Diagnostics','Envelope','Collisions','Movwindow',
                  'Loadbalancing','Reconfiguration','PartMerging','PartInjection']:
            if k in self.ax.data:
                d = self.ax.data[k] + [None]*(t.size-len(self.ax.data[k]))
                d = array(d, dtype="float")
                d[isnan(d)] = 0.
                y += [d]
                labels += [k]
        
        # Get branch names
        commits = self.ax.data["commit"]
        self.branches = array(["-".join(b.split("-")[1:]) for b in commits])
        self.branch_points = {}
        
        self.marker_properties = {}
        print("List of unique branches: ")
        for ibranch,branch in enumerate(unique(self.branches)):
            print(" - {}".format(branch))
            self.marker_properties[branch] = {}
            self.marker_properties[branch]["marker"] = custom_markers[ibranch]
            if "HEAD" in branch or branch=="develop":
                self.marker_properties[branch]["alpha"] = 0
                self.marker_properties[branch]["color"] = custom_colors[ibranch]
            else:
                self.marker_properties[branch]["alpha"] = 1
                self.marker_properties[branch]["color"] = custom_colors[ibranch]
            index = flatnonzero(self.branches == branch)
            self.branch_points[branch] = [t[index], array(self.ax.data["Timeintimeloop"],dtype="float")[index]]
        
        colors = []
        markers = []
        for ibranch,branch in enumerate(self.branches):
           colors.append(self.marker_properties[branch]["color"])
           markers.append(self.marker_properties[branch]["marker"])
        
        # Plot
        sca(self.ax)
        self.ax.cla()
        self.ax.xaxis_date()
        self.ax.xaxis.set_major_formatter(dates.DateFormatter('%d/%m %H:%M'))
        title(cases[self.ind])
        self.ax.stackplot( t, *y, labels=labels )
        self.ax.plot( t, array(self.ax.data["Timeintimeloop"]), '--k', label="time loop (total)",zorder=0)
        
        for i,branch in enumerate(self.branch_points):
            plot( self.branch_points[branch][0], self.branch_points[branch][1], linestyle="none", ms=self.marker_size, marker=self.marker_properties[branch]["marker"], label=branch, color=self.marker_properties[branch]["color"])
            
        sc = self.ax.scatter(t,array(self.ax.data["Timeintimeloop"],dtype="float"),alpha=0,picker=10,zorder=1)
        # Set the mqrker color according to the branch only if we use the scatter plot
        #sc.set_color(colors)
        
        self.ax.legend(bbox_to_anchor=(1.05, 1.))
        self.ax.set_xlabel("Commit date")
        self.ax.set_ylabel("Time (s)")
        setp(self.ax.get_xticklabels(), rotation=30, ha="right")
        
        #pick_event = fig.canvas.callbacks.connect('pick_event', on_pick)
        
        show()

# Launch the plotter in a switchPlots window
plots = switchPlots( len(cases) )
