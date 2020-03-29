from numpy import *
from matplotlib.pyplot import *
from matplotlib.widgets import Button
from matplotlib import dates
from glob import glob
from json import load
from time import strptime
from datetime import datetime
from os.path import splitext, basename
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

# Obtain the files
cases = sorted(glob("/sps2/gitlab-runner/logs/*.log"))
#cases = sorted(glob("./logs/*.log"))

# Class to create a plot with buttons to switch between cases
class switchPlots:
    
    # _________________________________________________
    # Class parameters
    
    # Internal parameters
    ind               = 0
    
    timers = ['Particles','Maxwell','Densities',
              'SyncParticles','SyncFields','SyncDensities','SyncSusceptibility',
              'Diagnostics','Envelope','Collisions','Movwindow',
              'Loadbalancing','Reconfiguration','PartMerging','PartInjection']
    
    # Plot style
    marker_size       = 8
    menu_x_position   = 0.75
    current_marker    = 0
    current_color     = 0
    custom_colors = ['C0','C1','C2','C3','C4','C5','C6']
    custom_markers = ['o','s','^','v','<','>']
    
    # _________________________________________________
    # Class methods
    
    def __init__(self, cases):
        
        self.nplots = len(cases)
        self.fig = figure(figsize=(14, 6))
        clf()
        gs = GridSpec(10,10)
        self.ax = subplot(gs[0:9,0:8])
        self.pick_event = self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        
        # Button next
        axnext = axes([self.menu_x_position+2*0.06, 0.01, 0.05, 0.03])
        self.bnext = Button(axnext, '>')
        self.bnext.on_clicked(lambda event : self.next(event,1))

        # Button previous
        axprev = axes([self.menu_x_position+0.06, 0.01, 0.05, 0.03])
        self.bprev = Button(axprev, '<')
        self.bprev.on_clicked(lambda event : self.next(event,-1))
   
        # Button next +5
        axnext5 = axes([self.menu_x_position+3*0.06, 0.01, 0.05, 0.03])
        self.bnext5 = Button(axnext5, '>>')
        self.bnext5.on_clicked(lambda event : self.next(event,5))

        # Button previous -5
        axprev5 = axes([self.menu_x_position, 0.01, 0.05, 0.03])
        self.bprev5 = Button(axprev5, '<<')
        self.bprev5.on_clicked(lambda event : self.next(event,-5))
        
        # Get all data from all cases
        self.data = []
        for case in cases:
            D = {"name":splitext(basename(case))[0]}
            
            # Read data
            with open(case, 'r') as f:
                data = load(f)
            
            # Get time
            D["date"] = [datetime(*strptime(d, "%Y_%m_%d_%H:%M:%S")[:6]) for d in data["date"]]
            D["t"] = array([dates.date2num(d) for d in D["date"]])
            
            # Get all timers
            D["processed_data"] = []
            D["labels"] = []
            D["min_times"] = []
            D["mean_times"] = []
            D["max_times"] = []
            for k in self.timers:
                if k in data:
                    d = data[k] + [None]*(D["t"].size-len(data[k]))
                    d = array(d, dtype="float")
                    d[isnan(d)] = 0.
                    D["processed_data"] += [d]
                    D["min_times"] += [np.min(d)]
                    D["mean_times"] += [np.mean(d)]
                    D["max_times"] += [np.max(d)]
                    D["labels"] += [k]
            
            # Array time_in_timeloop that contains the full time
            D["time_in_timeloop"] = array(data["Timeintimeloop"], dtype="float")
            D["min_times" ] += [np.min (D["time_in_timeloop"])]
            D["mean_times"] += [np.mean(D["time_in_timeloop"])]
            D["max_times" ] += [np.max (D["time_in_timeloop"])]
            
            # Get branch names
            D["commits"] = data["commit"]
            D["branches"] = array(["-".join(b.split("-")[1:]) for b in D["commits"]])
            # Set parameters for branch plotting 
            D["branch_opt"] = {}
            D["branch_points"] = {}
            for branch in unique(D["branches"]):
                if "HEAD" in branch or branch=="develop":
                    D["branch_opt"][branch] = dict(
                        marker = self.get_marker(),
                        color = self.get_color(),
                        alpha = 0,
                    )
                else:
                    D["branch_opt"][branch] = dict(
                        marker = self.get_marker(),
                        color = self.get_color(),
                        alpha = 1,
                        label = branch
                    )
                index = flatnonzero(D["branches"] == branch)
                D["branch_points"][branch] = [D["t"][index], D["time_in_timeloop"][index]]
            
            # Store all info for that case
            self.data += [D]
            print("Loaded case %s"%case)
        
        self.plot()
    
    def next(self, event, jump):
        """
        This method enables to switch to another log jumping of `jump` indexes from the current log.
        """
        self.ind += jump
        self.ind %= self.nplots
        self.plot()

    def get_marker(self):
        """
        Select a marker in the custom list
        """
        marker = self.current_marker
        self.current_marker = (self.current_marker+1)%len(self.custom_markers)
        return self.custom_markers[marker]

    def get_color(self):
        """
        Select a color in the custom list
        """
        color = self.current_color
        self.current_color = (self.current_color+1)%len(self.custom_colors)
        return self.custom_colors[color]

    def on_pick(self,event):
        """
        This method is used to make markers clickable.
        Information is then shown in the terminal.
        """
        artist = event.artist
        xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
        #x, y = artist.get_xdata(), artist.get_ydata()
        marker_indexes = event.ind
        print("\n Selected points: {}".format(marker_indexes))
        D = self.data[self.ind]
        for k in marker_indexes:
            print("  > Commit: {}".format(D["commits"][k]))
            print("    Branch: {}".format(D["branches"][k]))
            print("    Date: {}".format(D["date"][k]))
            print("    ----------------------------------------------------------------------")
            print("     Timers          | Times (s)  | Min (s)    | Mean (s)   | Max (s)    |")
            print("    ----------------------------------------------------------------------")
            for label, d, min, mean, max in zip(
                D["labels"]+["Total"],
                D["processed_data"]+[D["time_in_timeloop"]],
                D["min_times"],
                D["mean_times"],
                D["max_times"]
            ):
                print("     {0:15} | {1:.4e} | {2:.4e} | {3:.4e} | {4:.4e} |".format(label, d[k], min, mean, max))

    def plot(self):
        """
        Method to plot the log of index self.ind.
        """
        
        D = self.data[self.ind]
        
        print("\n _____________________________________________________ ")
        print(" Benchmark {}\n".format(D["name"]))
        print(" List of unique branches: ")
        for branch in unique(D["branches"]):
            print(" - {}".format(branch))
        
        # Plot
        sca(self.ax)
        self.ax.cla()
        self.ax.xaxis_date()
        self.ax.xaxis.set_major_formatter(dates.DateFormatter('%d/%m %H:%M'))
        self.ax.set_title(D["name"])
        self.ax.stackplot( D["t"], *D["processed_data"], labels=D["labels"] )
        self.ax.plot( D["t"], D["time_in_timeloop"], '--k', label="time loop (total)",zorder=0)
        
        for branch in D["branch_points"]:
            plot( 
                D["branch_points"][branch][0], D["branch_points"][branch][1],
                linestyle = "none",
                ms = self.marker_size,
                **D["branch_opt"][branch]
            )
            
        sc = self.ax.scatter(D["t"], D["time_in_timeloop"], alpha=0, picker=10, zorder=1)
        
        self.ax.legend(bbox_to_anchor=(1.05, 1.))
        self.ax.set_xlabel("Commit date")
        self.ax.set_ylabel("Time (s)")
        setp(self.ax.get_xticklabels(), rotation=30, ha="right")
        
        show()

plots = switchPlots( cases )
