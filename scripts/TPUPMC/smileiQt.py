#!/usr/bin/env python 
"""
Plot fields of smilei simulaition
"""
import sys, os, random
from PyQt4 import QtCore, QtGui, uic

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

import os
import tables
import numpy as np

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class smileiQt(QtGui.QMainWindow):

    def __init__(self):
        super(smileiQt, self).__init__()

        uiFile=os.path.dirname(os.path.realpath(__file__))+'/smileiQt.ui'
        self.ui=uic.loadUi(uiFile,self)
        print uiFile
        self.ui.actionQuit.triggered.connect(QtGui.qApp.quit)

        self.figure, ax = plt.subplots()
        
        self.figure, axarr = plt.subplots(4, sharex=True)

        self.canvas = FigureCanvas(self.figure)

        self.toolbar = NavigationToolbar(self.canvas, self)

        self.plotGrid.addWidget(self.toolbar,0,0)
        self.plotGrid.addWidget(self.canvas,1,0)


        x = np.linspace(0, 2 * np.pi, 400)
        y = np.sin(x ** 2)

        plt.close('all')

        # Just a figure and one subplot
        ax.plot(x, y)
        ax.set_title('Simple plot')

        # Two subplots, the axes array is 1-d
        axarr[0].plot(x, y)
        axarr[0].set_title('Sharing X axis')
        axarr[1].scatter(x, y)
        axarr[1].set_title('pippo')
        axarr[2].scatter(x, y)
        axarr[2].set_title('pappa')
        axarr[3].scatter(x, y)
        axarr[3].set_title('pippa')

        self.canvas.draw()
    
#         self.canvas.mpl_connect('motion_notify_event', self.on_movement)

        self.show()

#     def on_movement(self, event):
#         if not (event.inaxes is None) :
#             msg = "(%f,%f)" % (event.xdata, event.ydata)
#             QtGui.QApplication.clipboard().setText(msg)
#             self.statusBar().showMessage(msg)

    def on_playStop_released(self):
        if self.playStop.isChecked() :
            self.playStop.setText("Stop")
        else:
            self.playStop.setText("Play")
        print "here"
        

def main():
    app = QtGui.QApplication(sys.argv)
    ex = smileiQt()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()  
