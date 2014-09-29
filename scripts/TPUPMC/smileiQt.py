#!/usr/bin/env python 
"""
Plot fields of smilei simulaition
"""
import sys, os, random
from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtGui import QFileDialog

import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

import os
import tables as tb
import numpy as np
import re

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
        self.timer.setInterval(50)
        self.timer.timeout.connect(self.do_timer)
        
        self.readData()

        self.show()
    
    def on_slider_valueChanged(self,step):
        self.step=step
        self.doPlots()
    
    def on_back_released(self):
        self.ui.slider.setValue(self.step-1)

    def on_forward_pressed(self):
        self.ui.slider.setValue(self.step+1)

    def on_allBack_released(self):
        self.ui.slider.setValue(0)

    def on_allForward_released(self):
        self.ui.slider.setValue(len(self.fieldSteps))

    def on_playStop_released(self):
        if self.playStop.isChecked() :
            self.playStop.setText("Stop")
            self.timer.start()
        else:
            self.playStop.setText("Play")
            self.timer.stop()
        
    def load_settings(self):
        settings=QtCore.QSettings("smileiQt","");
        settings.beginGroup("Preferences");
        self.dirName=str(settings.value("dirName",".").toString());
        settings.endGroup();
    
    def doPlots(self):
        print "doPlots"
        self.step=max(0,min(self.step,len(self.fieldSteps)-1))
        self.slider.setValue(self.step)
        
        if not hasattr(self,'figure'): return
        
        nplot=0
    
        fname=os.path.join(self.dirName, "scalars.txt")
        if os.path.isfile(fname) :
            file=np.loadtxt(fname)
            for j in self.ui.menuScalars.actions():
                if j.isChecked() :
                    col=j.data().toInt()[0]
                    x=file[:,0]
                    y=file[:,col]
                    self.axarr[nplot][0].plot(x,y)
                    self.axarr[nplot][0].set_title(j.text())
                    
                    self.axarr[nplot][0].axvline(x=self.step/self.res_time,c="red",linewidth=2,zorder=0, clip_on=False)
                    nplot+=1

        fname=os.path.join(self.dirName, "Fields.h5")
        if os.path.isfile(fname) :
            f=tb.openFile(fname)
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
                    self.axarr[nplot][0].set_title(i.text() + " %.3f" % (self.step/self.res_time) )
                    if not self.autoScale.isChecked() :
                        self.axarr[nplot][0].set_ylim(oldRange)
                        
                    self.axarr[nplot][0].set_xlim(0,self.sim_length)

                    nplot+=1
                
            f.close()

        
        fname=os.path.join(self.dirName, "PhaseSpace.h5")
        if os.path.isfile(fname) :
            f=tb.openFile(fname)

            for i in self.ui.menuPhase_spaces.actions() :
                if i.isChecked() :
                    nameData=str(i.data().toString())
                    node=f.getNode(str(nameData))
                    print "----->"
                    data=node.read()[self.step,:,:].T
                    print "<-----"
                    print data.shape
                    print data.min(), data.max()
                    extents=node._v_parent._v_attrs.extents.reshape(4).tolist()
                    self.axarr[nplot][0].imshow(data,extent=extents, aspect='auto')
#                     cb=plt.colorbar(im)
                    nplot+=1
                
                
            f.close()

        if nplot>0 : 
            self.canvas.draw()

    def do_timer(self):
        self.ui.slider.setValue(self.step+1)
        

    def readData(self):
    
        allOk=True
        
        fname=os.path.join(self.dirName, "scalars.txt")
        if os.path.isfile(fname) :
            names=[]
            for line in open(fname):
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
        else:
            allOk=False
        
        fname=os.path.join(self.dirName, "Fields.h5")
        if os.path.isfile(fname) :
            self.fieldSteps=[]
            f=tb.openFile(os.path.join(self.dirName, fname))
            
            self.res_time=f.root._v_attrs.res_time
            self.sim_length=f.root._v_attrs.sim_length
            
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
            self.slider.setRange(0,len(self.fieldSteps)-1)
        else:
            allOk=False

        fname=os.path.join(self.dirName, "PhaseSpace.h5")
        if os.path.isfile(fname) :
            f=tb.openFile(os.path.join(self.dirName, fname))
            for phaseGroup in f.walkNodes("/", classname='Array'):
                my_act= QtGui.QAction(re.sub('/',' ',phaseGroup._v_pathname),self)
                my_act.setData(phaseGroup._v_pathname)
                my_act.setCheckable(True)
                self.ui.menuPhase_spaces.addAction(my_act)
                my_act.triggered.connect(self.action_clicked)
            f.close()
        else:
            allOk=False

        if allOk:
            self.doPlots()
        
    def rescalePlots(self):
        self.nplots=0;
        menus=[self.ui.menuScalars,self.ui.menuFields,self.ui.menuPhase_spaces]
        for j in menus :
            for i in j.actions():
                if i.isChecked() :
                    self.nplots+=1
                    
        menus=[self.ui.menuPhase_spaces]

#         self.cbs=[]
#         for i in self.ui.menuPhase_spaces.actions():
#             self.cbs.append(i.isChecked())
                    
        
        
        if self.nplots > 0: 
        
            self.figure, self.axarr = plt.subplots(self.nplots,1, squeeze=False)
            
            for i in self.axarr:
                i[0].hold(False)


            for i in reversed(range(self.plot.layout().count())): 
                self.plot.layout().itemAt(i).widget().setParent(None)

            self.canvas = FigureCanvas(self.figure)
            toolbar = NavigationToolbar(self.canvas, self)
            toolbar.setFixedHeight(18)
            self.plot.layout().addWidget(toolbar)
            self.plot.layout().addWidget(self.canvas)



    def action_clicked(self):
        print "action_clicked"
        self.rescalePlots()
        self.doPlots()
#         print self.sender().data().toInt()
        
    def on_changeDir(self):
        dirName=QtGui.QFileDialog.getExistingDirectory(self,self.dirName, options=QFileDialog.ShowDirsOnly)
        if not dirName.isEmpty():
            self.dirName=str(dirName)
            self.readData()

def main():

    app = QtGui.QApplication(sys.argv)
    ex = smileiQt()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()  
