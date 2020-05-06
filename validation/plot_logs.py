from numpy import *
from matplotlib.pyplot import *
from matplotlib.widgets import Button
from matplotlib import dates
from glob import glob
from json import load
from time import strptime
from datetime import datetime
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
        
        # Events
        self.fig.canvas.callbacks.connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_hover)
        self.fig.canvas.mpl_connect("motion_notify_event", self.update_min_mean_max_annotations)
        self.fig.canvas.mpl_connect("button_press_event", self.open_commit_link)
        
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
            
            self.init_style()
            
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
            # Set parameters for branch plotting
            D["branch_opt"] = {}
            D["branch_points"] = {}
            for branch in unique(D["branches"]):
                if "HEAD" in branch or branch=="develop":
                    D["branch_opt"][branch] = dict(
                        marker = self.get_marker(),
                        color = "k",
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

    def init_style(self):
        """
        Reinitialize the plot style
        """
        self.current_color  = 0
        self.current_marker = 0

    def on_pick(self,event):
        """
        This method is used to make markers clickable.
        Information is then shown in the terminal.
        """
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

    
    def on_hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            cont, ind = self.sc.contains(event)
            if cont:
                i = ind["ind"][0]
                self.annot.set_text(self.data[self.ind]["commits"][i])
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()

    def open_commit_link(self,event):
        """
        """
        if event.dblclick:
            if event.inaxes == self.ax:
                cont, ind = self.sc.contains(event)
                if cont:
                    D = self.data[self.ind]
                    for i in  ind["ind"]:
                        os.system("xdg-open https://llrgit.in2p3.fr/smilei/smilei/-/commit/{}".format(D["commit_ids"][i]))

    def display_min_mean_max_annotations(self):
        """
        Display 3 arrows with annotations for the min, max and mean value of total time
        """
        D = self.data[self.ind]
        ylim = self.ax.get_ylim()
        y_length = ylim[1] - ylim[0]
        
        if ((D["min_times"][-1] < ylim[1]) and (D["min_times"][-1] > ylim[0])):
            ypos = (D["min_times"][-1]-ylim[0])/y_length
            self.min_annot = self.ax.annotate("min", xy=(1.0, ypos),
                                            xycoords='axes fraction',
                                            xytext=(1.03, ypos),
                                            arrowprops=dict(arrowstyle="-|>",color='C0'),
                                            color = 'C0',
                                            va = 'center')
            self.min_annot_drawn = True
                    
        if ((D["mean_times"][-1] < ylim[1]) and (D["mean_times"][-1] > ylim[0])):
            ypos = (D["mean_times"][-1]-ylim[0])/y_length
            self.mean_annot = self.ax.annotate("mean", xy=(1.0, ypos),
                                            xycoords='axes fraction',
                                            xytext=(1.03, ypos),
                                            arrowprops=dict(arrowstyle="-|>",color='C2'),
                                            color = 'C2',
                                            va = 'center')
            self.mean_annot_drawn = True
        
        if ((D["max_times"][-1] < ylim[1]) and (D["max_times"][-1] > ylim[0])):
            ypos = (D["max_times"][-1]-ylim[0])/y_length
            self.max_annot = self.ax.annotate("max", xy=(1.0, ypos),
                                            xycoords='axes fraction',
                                            xytext=(1.03, ypos),
                                            arrowprops=dict(arrowstyle="-|>",color='C3'),
                                            color = 'C3',
                                            va = 'center')
            self.max_annot_drawn = True
    
    def update_min_mean_max_annotations(self,event):
        """
        Update the min, mean, max annotations
        """
        if (self.min_annot_drawn):
            self.min_annot.remove()
            self.min_annot_drawn = False
        if (self.mean_annot_drawn):
            self.mean_annot.remove()
            self.mean_annot_drawn = False
        if (self.max_annot_drawn):
            self.max_annot.remove()
            self.max_annot_drawn = False
        
        self.display_min_mean_max_annotations()
    
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
        
        self.ax.annotate("Help:\n> Fly over markers to get commit id\n"+
                         "> Click on markers to get full information\n    in the terminal\n"+
                         "> Use arrows to change log plot", xy=(1.05,-0.05), xycoords='axes fraction',
                         bbox=dict(boxstyle="round",ec="#aaaaaa",fc="None",pad=1.0)
                         )
        
        show()

plots = switchPlots( cases )
