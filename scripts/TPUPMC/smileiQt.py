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
import tables as tb
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

        self.doPlots()
        
    def load_settings(self):
        settings=QtCore.QSettings("smileiQt","");
        settings.beginGroup("Preferences");
        self.dirName=str(settings.value("dirName",".").toString());
        settings.endGroup();

    def numPlots(self):
        num=0;
        menus=[self.ui.menuScalars,self.ui.menuFields,self.ui.menuPhase_spaces,self.ui.menuProbes]
        for j in menus :
            for i in j.actions():
                if i.isChecked() :
                    num+=1
        return num
    
    def doPlots(self):

        plt.close('all')
        
        nplots=self.numPlots()
        if nplots > 0: 
            figure, axarr = plt.subplots(nplots,1, sharex=True, squeeze=False)
            canvas = FigureCanvas(figure)
            toolbar = NavigationToolbar(canvas, self)

            self.plotGrid.addWidget(toolbar,0,0)
            self.plotGrid.addWidget(canvas,1,0)
        
            nplot=0
        
            fname="scalars.txt"
            if os.path.isfile(fname) :
                file=np.loadtxt(fname)
                print file.shape
                for j in self.ui.menuScalars.actions():
                    if j.isChecked() :
                        col=j.data().toInt()[0]
                        x=file[:,0]
                        y=file[:,col]
                    
                        print x, y
                        
                        axarr[nplot][0].plot(x,y)
                        axarr[nplot][0].set_title(j.text())
                        nplot+=1


            canvas.draw()


#         self.axarr[0].plot(x, y)
#         self.axarr[0].set_title(self.dirName)
#         self.axarr[1].scatter(x, y)
#         self.axarr[1].set_title('pippo')
#         self.axarr[2].scatter(x, y)
#         self.axarr[2].set_title('pappa')
#         self.axarr[3].scatter(x, y)
#         self.axarr[3].set_title('pippa')
#         print "Here"
    
    def readData(self):
#         pippo= QtGui.QAction("pppp",self)
#         pippo.setCheckable(True)
#         self.ui.menuScalars.addAction(pippo)


        fname="scalars.txt"
        if os.path.isfile(fname) :

            names=[]
            for line in open(os.path.join(self.dirName, "scalars.txt")):
                li=line.strip()
                if li.startswith("#"):
                    list=line.split()
                    names.append(list[-1])
       
            scalars_names=names[1:-2]
            for i in range(len(scalars_names)):
                my_act= QtGui.QAction(scalars_names[i],self)
                my_act.setData(i+1)
                my_act.setCheckable(True)
                self.ui.menuScalars.addAction(my_act)
                my_act.triggered.connect(self.action_clicked)
        
        fname="Fields.h5"
        if os.path.isfile(fname) :
            f=tb.openFile(os.path.join(self.dirName, fname))
            for array in f.list_nodes("/0000000000", classname='Array'):
                name=array._v_name
                my_act= QtGui.QAction(name,self)
                my_act.setData(name)
                my_act.setCheckable(True)
                self.ui.menuFields.addAction(my_act)
                my_act.triggered.connect(self.action_clicked)
            f.close()

#         fname="PhaseSpace.h5"
#         if os.path.isfile(fname) :
#             f=tb.openFile(os.path.join(self.dirName, fname))
#             for i in f.root:
#                 for j in i:
#                     name=j._v_name
#                     my_act= QtGui.QAction(i._v_name+" "+name,self)
#                     my_act.setData(i._v_name+"/"+name)
#                     my_act.setCheckable(True)
#                     self.ui.menuPhase_spaces.addAction(my_act)
#                     my_act.triggered.connect(self.action_clicked)
#             f.close()
        
#         fname="Probes.h5"
#         if os.path.isfile(fname) :
#             f=tb.openFile(os.path.join(self.dirName, fname))
#             for i in f.root:
#                 my_menu = self.ui.menuProbes.addMenu(i._v_name)
#                 
#                 my_act= QtGui.QAction(i._v_name,self)
#                 my_act.setData(i._v_name)
#                 my_act.setCheckable(True)
#                 my_menu.addAction(my_act)
#                 
#                 my_act.triggered.connect(self.action_clicked)
#             f.close()
        
        self.doPlots()

    def action_clicked(self):
        self.doPlots()
#         print self.sender().data().toInt()
        
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
