from numpy import *
from matplotlib.pyplot import *
from matplotlib.widgets import Button
from matplotlib import dates
from glob import glob
from json import load
from time import strptime
from datetime import datetime
from itertools import cycle
#from os.path import splitext, basename
import os
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
if len(sys.argv) > 1:
    try:
        path = sys.argv[1]
        if os.path.isdir(path):
            path = path +"/*.log"
    except:
        raise Exception("\n Please, provide a valid path to the logs.\n")
else:
    path = "/sps2/gitlab-runner/logs/*.log"
    
cases = sorted(glob(path))
if len(cases) ==0:
    raise Exception('\n No cases found for `{}`.\n The path may be invalid or not accessible.'.format(path))

# Class to create a plot with buttons to switch between cases
class switchPlots:
    
    # _________________________________________________
    # Class parameters
    
    # Internal parameters
    ind = 0
    
    timers = ['Particles','Maxwell','Densities',
              'SyncParticles','SyncFields','SyncDensities','SyncSusceptibility',
              'Diagnostics','Envelope','Collisions','Movwindow',
              'Loadbalancing','Reconfiguration','PartMerging','PartInjection']
    
    # Plot style
    marker_size       = 8
    menu_x_position   = 0.75
    custom_colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
    custom_markers = ['o','s','^','v','<','>']
    
    # _________________________________________________
    # Class methods
    
    def __init__(self, cases):
        
        self.nplots = len(cases)
        self.fig = figure(figsize=(14, 6))
        clf()
        gs = GridSpec(10,10)
        self.ax = subplot(gs[0:9,0:8])
        
        # Events
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        
        # Button previous -5
        axprev5 = axes([self.menu_x_position, 0.01, 0.05, 0.03])
        self.bprev5 = Button(axprev5, '<<')
        self.bprev5.on_clicked(lambda event : self.next(event,-5))
        
        # Button previous
        axprev = axes([self.menu_x_position+0.06, 0.01, 0.05, 0.03])
        self.bprev = Button(axprev, '<')
        self.bprev.on_clicked(lambda event : self.next(event,-1))
        
        # Button next
        axnext = axes([self.menu_x_position+2*0.06, 0.01, 0.05, 0.03])
        self.bnext = Button(axnext, '>')
        self.bnext.on_clicked(lambda event : self.next(event,1))
   
        # Button next +5
        axnext5 = axes([self.menu_x_position+3*0.06, 0.01, 0.05, 0.03])
        self.bnext5 = Button(axnext5, '>>')
        self.bnext5.on_clicked(lambda event : self.next(event,5))
        
        
        # Get all data from all cases
        self.data = []
        for case in cases:
            
            D = {"name":os.path.splitext(os.path.basename(case))[0]}
            
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
            D["commit_ids"] = array([(b.split("-")[0]) for b in D["commits"]])
            D["branches"] = array(["-".join(b.split("-")[1:]) for b in D["commits"]])
            
            # Find 8 newest branches
            recent_branches = []
            for branch in D["branches"][::-1]:
                if len(recent_branches) > 7:
                    break
                if branch not in recent_branches:
                    recent_branches.append(branch)
            
            # Set parameters for branch plotting
            D["branch_opt"] = {}
            D["branch_points"] = {}
            markers = cycle(self.custom_markers)
            colors = cycle(self.custom_colors)
            for branch in unique(D["branches"]):
                self.current_color  = 0
                self.current_marker = 0
                if branch not in recent_branches:
                    D["branch_opt"][branch] = dict(
                        marker = "s",
                        color = "k",
                        alpha = 0,
                    )
                elif "HEAD" in branch or branch=="develop":
                    D["branch_opt"][branch] = dict(
                        marker = next(markers),
                        color = "k",
                        alpha = 1,
                        label = "develop or HEAD"
                    )
                else:
                    D["branch_opt"][branch] = dict(
                        marker = next(markers),
                        color = next(colors),
                        alpha = 1,
                        label = branch
                    )
                index = flatnonzero(D["branches"] == branch)
                D["branch_points"][branch] = [D["t"][index], D["time_in_timeloop"][index]]
            
            # Store all info for that case
            self.data += [D]
            print("Loaded case %s"%case)
        
        self.plot_summary()
    
    def next(self, event, jump):
        """
        This method enables to switch to another log jumping of `jump` indexes from the current log.
        """
        self.ind += jump
        self.ind %= self.nplots
        self.plot_benchmark()
    
    def on_pick(self,event):
        """
        This method is used to make markers pickable.
        """
        if self.plottype == "summary":
            pass
        
        elif self.plottype == "benchmark":
            marker_indexes = event.ind
            print("\n Selected points: {}".format(marker_indexes))
            D = self.data[self.ind]
            for k in marker_indexes:
                print("  > Commit: {}".format(D["commits"][k]))
                print("    Branch: {}".format(D["branches"][k]))
                print("    Date: {}".format(D["date"][k]))
                print("    Link: https://llrgit.in2p3.fr/smilei/smilei/-/commit/{}".format(D["commit_ids"][k]))
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
            
        else:
            raise Exception("Impossible")
    
    def on_hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            
            if self.plottype == "summary":
                for plot in self.plots:
                    cont, _ = plot[0].contains(event)
                    if cont:
                        self.annot.set_text(plot[0].get_label())
                        self.annot.set_visible(True)
                        self.fig.canvas.draw_idle()
                        break
                else:
                    if vis:
                        self.annot.set_visible(False)
                        self.fig.canvas.draw_idle()
            
            elif self.plottype == "benchmark":
                cont, ind = self.sc.contains(event)
                if cont:
                    i = ind["ind"][0]
                    self.annot.set_text(self.data[self.ind]["commits"][i])
                    self.annot.set_visible(True)
                    self.fig.canvas.draw_idle()
                elif vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()
                
                self.update_min_mean_max_annotations()
            
            else:
                raise Exception("Impossible")

    def on_press(self,event):
        """
        """
        if event.dblclick and event.inaxes == self.ax:
            
            if self.plottype == "summary":
                for plot in self.plots:
                    cont, ind = plot[0].contains(event)
                    if cont:
                        break
                self.ind = plot[1]
                self.plot_benchmark()
            
            elif self.plottype == "benchmark":
                cont, ind = self.sc.contains(event)
                if cont:
                    D = self.data[self.ind]
                    for i in ind["ind"]:
                        os.system("xdg-open https://llrgit.in2p3.fr/smilei/smilei/-/commit/{}".format(D["commit_ids"][i]))
                
            else:
                raise Exception("Impossible")

    def display_min_mean_max_annotations(self):
        """
        Display 3 arrows with annotations for the min, max and mean value of total time
        """
        D = self.data[self.ind]
        ylim = self.ax.get_ylim()
        y_length = ylim[1] - ylim[0]
        
        if D["min_times"][-1] < ylim[1] and D["min_times"][-1] > ylim[0]:
            ypos = (D["min_times"][-1]-ylim[0])/y_length
            self.min_annot = self.ax.annotate("min", xy=(1.0, ypos),
                                            xycoords='axes fraction',
                                            xytext=(1.03, ypos),
                                            arrowprops=dict(arrowstyle="-|>",color='C0'),
                                            color = 'C0',
                                            va = 'center')
            self.min_annot_drawn = True
                    
        if D["mean_times"][-1] < ylim[1] and D["mean_times"][-1] > ylim[0]:
            ypos = (D["mean_times"][-1]-ylim[0])/y_length
            self.mean_annot = self.ax.annotate("mean", xy=(1.0, ypos),
                                            xycoords='axes fraction',
                                            xytext=(1.03, ypos),
                                            arrowprops=dict(arrowstyle="-|>",color='C2'),
                                            color = 'C2',
                                            va = 'center')
            self.mean_annot_drawn = True
        
        if D["max_times"][-1] < ylim[1] and D["max_times"][-1] > ylim[0]:
            ypos = (D["max_times"][-1]-ylim[0])/y_length
            self.max_annot = self.ax.annotate("max", xy=(1.0, ypos),
                                            xycoords='axes fraction',
                                            xytext=(1.03, ypos),
                                            arrowprops=dict(arrowstyle="-|>",color='C3'),
                                            color = 'C3',
                                            va = 'center')
            self.max_annot_drawn = True
    
    def update_min_mean_max_annotations(self):
        """
        Update the min, mean, max annotations
        """
        if self.min_annot_drawn:
            self.min_annot.remove()
            self.min_annot_drawn = False
        if self.mean_annot_drawn:
            self.mean_annot.remove()
            self.mean_annot_drawn = False
        if self.max_annot_drawn:
            self.max_annot.remove()
            self.max_annot_drawn = False
        
        self.display_min_mean_max_annotations()
    
    def plot_summary(self):
        self.plottype = "summary"
        
        # Put buttons away
        self.bprev5.ax.set_position([0.,0.,0.,0.])
        self.bprev .ax.set_position([0.,0.,0.,0.])
        self.bnext .ax.set_position([0.,0.,0.,0.])
        self.bnext5.ax.set_position([0.,0.,0.,0.])
        
        sca(self.ax)
        self.ax.cla()
        self.ax.xaxis_date()
        self.ax.xaxis.set_major_formatter(dates.DateFormatter('%d/%m %H:%M'))
        self.ax.set_title("All benchmarks")
        self.ax.set_ylabel("relative time")
        self.plots = []
        for i, D in enumerate(self.data):
            ncommits = min(40, len(D["t"]))
            med = median(D["time_in_timeloop"][-ncommits:])
            relative_time = D["time_in_timeloop"] / med
            self.plots.append( self.ax.plot( D["t"], relative_time, '--', label=D["name"]) + [i] )
        self.ax.set_ylim(0,3)
        
        self.ax.annotate(
            """
Help:
> Fly over line to see benchmark
> Double-click on line to get benchmark detail
            """,
            xy=(1.05,-0.05), xycoords='axes fraction',
            bbox=dict(boxstyle="round",ec="#aaaaaa",fc="None",pad=1.0)
        )
        
        # Hovering box
        self.annot = self.ax.annotate("", xy=(0.01,0.95), xycoords='axes fraction')
        self.annot.set_visible(False)
    
    def plot_benchmark(self):
        """
        Method to plot the log of index self.ind.
        """
        
        self.plottype = "benchmark"
        
        # Put buttons in their place
        self.bprev5.ax.set_position([self.menu_x_position, 0.01, 0.05, 0.03])
        self.bprev .ax.set_position([self.menu_x_position+0.06, 0.01, 0.05, 0.03])
        self.bnext .ax.set_position([self.menu_x_position+2*0.06, 0.01, 0.05, 0.03])
        self.bnext5.ax.set_position([self.menu_x_position+3*0.06, 0.01, 0.05, 0.03])
        
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
            self.ax.plot(
                D["branch_points"][branch][0], D["branch_points"][branch][1],
                linestyle = "none",
                ms = self.marker_size,
                **D["branch_opt"][branch]
            )
        
        self.sc = self.ax.scatter(D["t"], D["time_in_timeloop"], alpha=0, picker=10, zorder=1)
        
        self.ax.legend(bbox_to_anchor=(1.08, 1.))
        self.ax.set_xlabel("Commit date")
        self.ax.set_ylabel("Time (s)")
        setp(self.ax.get_xticklabels(), rotation=30, ha="right")
        
        # Min, mean, max annotations
        self.display_min_mean_max_annotations()
        
        # Hovering box
        self.annot = self.ax.annotate("", xy=(0.01,0.95), xycoords='axes fraction')
        self.annot.set_visible(False)
        
        self.ax.annotate(
            """
Help:
> Fly over markers to get commit id
> Click on markers to get full information
    in the terminal
> Use arrows to change log plot
            """,
            xy=(1.05,-0.05), xycoords='axes fraction',
            bbox=dict(boxstyle="round",ec="#aaaaaa",fc="None",pad=1.0)
        )
        
        show()

plots = switchPlots( cases )
