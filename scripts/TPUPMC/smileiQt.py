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


        l = QtGui.QVBoxLayout()
        self.plot.setLayout(l)

        self.step=0
        
        self.timer=QtCore.QTimer()
        self.timer.timeout.connect(self.do_timer)
        
        self.readData()

        self.show()

    def on_back_released(self):
        self.step-=1
        self.doPlots()

    def on_forward_released(self):
        self.step+=1
        self.doPlots()

    def on_allBack_released(self):
        self.step=0
        self.doPlots()

    def on_allForward_released(self):
        self.step=len(self.fieldSteps)
        self.doPlots()

    def on_playStop_released(self):
        if self.playStop.isChecked() :
            self.playStop.setText("Stop")
            self.timer.stop()
        else:
            self.playStop.setText("Play")
            self.timer.start(100)
        
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

        self.step=max(0,min(self.step,len(self.fieldSteps)))

        plt.close('all')
        if hasattr(self, 'figure') :
            figure.clf()
        plt.clf()
        nplots=self.numPlots()
        if nplots > 0: 
            nplot=0
        
            fname="scalars.txt"
            if os.path.isfile(fname) :
                file=np.loadtxt(fname)
                for j in self.ui.menuScalars.actions():
                    if j.isChecked() :
                        col=j.data().toInt()[0]
                        x=file[:,0]
                        y=file[:,col]
                        self.axarr[nplot][0].plot(x,y)
                        self.axarr[nplot][0].set_title(j.text())
                        nplot+=1

            fname="Fields.h5"

            if os.path.isfile(fname) :
                f=tb.openFile(os.path.join(self.dirName, fname))
                nameGroup="/%010d" % self.step
                for i in self.ui.menuFields.actions() :
                    if i.isChecked() :
                        nameData=nameGroup+"/"+i.text()
                    
                        data=f.getNode(str(nameData)).read()
                        x=np.array(range(len(data)))/f.root._v_attrs.res_space[0]
                        y=data
                        if not self.autoScale.isChecked() :
                            oldRange=self.axarr[nplot][0].get_ylim()
                        
                        self.axarr[nplot][0].plot(x,y)
                        self.axarr[nplot][0].set_title(i.text() + " %.3f" % (self.step/f.root._v_attrs.res_time) )
                        if not self.autoScale.isChecked() :
                            self.axarr[nplot][0].set_ylim(oldRange)

                        nplot+=1
                    
                f.close()
                
            self.canvas.draw()

    def do_timer(self):        
        self.step+=1

        self.doPlots()
        

    def readData(self):
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
            self.fieldSteps=[]
            f=tb.openFile(os.path.join(self.dirName, fname))
            for group in f.list_nodes("/", classname='Group'):
                self.fieldSteps.append(group._v_name)

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
        
    def rescalePlots(self):
        nplots=self.numPlots()
        
        if nplots > 0: 
        
            figure, self.axarr = plt.subplots(nplots,1, squeeze=False)
            
            for i in self.axarr:
                i[0].hold(False)


            for i in reversed(range(self.plot.layout().count())): 
                self.plot.layout().itemAt(i).widget().setParent(None)

            self.canvas = FigureCanvas(figure)
            toolbar = NavigationToolbar(self.canvas, self)
            toolbar.setFixedHeight(18)
            self.plot.layout().addWidget(toolbar)
            self.plot.layout().addWidget(self.canvas)



    def action_clicked(self):
        self.rescalePlots()
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
