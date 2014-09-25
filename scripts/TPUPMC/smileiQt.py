#!/usr/bin/env python 
"""
Plot fields of smilei simulaition
"""
import sys, os, random
from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtGui import QFileDialog

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import matplotlib.pyplot as plt

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

        self.load_settings()

        uiFile=os.path.dirname(os.path.realpath(__file__))+'/smileiQt.ui'
        self.ui=uic.loadUi(uiFile,self)
        print uiFile
        self.ui.actionQuit.triggered.connect(QtGui.qApp.quit)
        self.ui.actionDir.triggered.connect(self.on_changeDir)

        self.show()
        
        self.readData()

    def on_playStop_released(self):
        if self.playStop.isChecked() :
            self.playStop.setText("Stop")
        else:
            self.playStop.setText("Play")
        print self.playStop.text()
        
    def load_settings(self):
        settings=QtCore.QSettings("smileiQt","");
        settings.beginGroup("Preferences");
        self.dirName=str(settings.value("dirName",".").toString());
        settings.endGroup();

    def doPlots(self):
        x = np.linspace(0, 2 * np.pi, 100)
        y = np.sin(x ** 2)

        plt.close('all')

        self.figure, self.axarr = plt.subplots(4, sharex=True)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.plotGrid.addWidget(self.toolbar,0,0)
        self.plotGrid.addWidget(self.canvas,1,0)

        # Two subplots, the axes array is 1-d
        self.axarr[0].plot(x, y)
        self.axarr[0].set_title(self.dirName)
        self.axarr[1].scatter(x, y)
        self.axarr[1].set_title('pippo')
        self.axarr[2].scatter(x, y)
        self.axarr[2].set_title('pappa')
        self.axarr[3].scatter(x, y)
        self.axarr[3].set_title('pippa')
        self.canvas.draw()
        print "Here"
    
    def readData(self):
        pippo= QtGui.QAction("pppp",self)
        pippo.setCheckable(True)
        self.ui.menuScalars.addAction(pippo)

        print "here"
        self.doPlots()
        
    def on_changeDir(self):
        dirName=QtGui.QFileDialog.getExistingDirectory(self,self.dirName, options=QFileDialog.ShowDirsOnly)
        if not dirName.isEmpty():
            self.dirName=dirName
            self.readData()

def main():

    app = QtGui.QApplication(sys.argv)
    ex = smileiQt()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()  
