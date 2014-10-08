#!/usr/bin/env python 
"""
Plot fields of smilei simulaition
"""
import sys, os, random
from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtGui import QFileDialog

import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib import gridspec

import matplotlib.pyplot as plt

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
        self.ui.actionQuit.triggered.connect(QtGui.qApp.quit)
        
        
        l = QtGui.QVBoxLayout()
        self.plot.setLayout(l)

        self.step=0
        
        self.timer=QtCore.QTimer()
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.do_timer)
        
        self.ui.autoScale.toggled.connect(self.doPlots)
        
        self.fig = plt.figure()
        self.canvas = FigureCanvas(self.fig)
        
        self.canvas.mpl_connect('motion_notify_event', self.on_movement)
        
        toolbar = NavigationToolbar(self.canvas, self)
        toolbar.setFixedHeight(18)
        self.plot.layout().addWidget(toolbar)
        self.plot.layout().addWidget(self.canvas)

        self.scalarData = None
        self.fieldFile = None
        self.phaseFile = None

        self.createActions()

        self.show()        
    
    def on_movement(self, event):
        if not (event.inaxes is None) :
            msg = "%.3f %.3f" % (event.xdata, event.ydata)
            self.ui.mouse.setText(msg)
        
    def closeEvent(self,event):
        if not self.fieldFile is None :
            self.fieldFile.close()
        if not self.phaseFile is None :
            self.phaseFile.close()
    
    def on_slider_valueChanged(self,step):
        self.step=step
        self.doPlots()
    
    @QtCore.pyqtSignature("int")
    def on_spinStep_valueChanged(self,my_step):
        self.ui.slider.setValue(my_step)
        
    def on_back_released(self):
        self.ui.slider.setValue(self.step-1)

    def on_forward_pressed(self):
        self.ui.slider.setValue(self.step+1)

    def on_allBack_released(self):
        self.ui.slider.setValue(0)

    def on_allForward_released(self):
        self.ui.slider.setValue(len(self.fieldSteps)-1)

    def on_playStop_toggled(self):
        if self.playStop.isChecked() :
            self.playStop.setText("Stop")
            self.timer.start()
        else:
            self.playStop.setText("Play")
            self.timer.stop()
        
    def load_settings(self):
        settings=QtCore.QSettings("smileiQt","")
        settings.beginGroup("Preferences")
        self.dirName=str(settings.value("dirName",".").toString())
        settings.endGroup()
    
    def doPlots(self):
        self.step %= len(self.fieldSteps)
        self.slider.setValue(self.step)
        
        if not hasattr(self,'fig'): return
        
        self.ui.spinStep.setValue(self.step)
        time=self.step/self.res_time*self.fieldEvery
        self.fig.suptitle("Time: %.3f" % time)
        nplot=0
    
        if not self.scalarData is None :
            for j in self.ui.menuScalars.actions():
                if j.isChecked() :
                    col=j.data().toInt()[0]
                    x=self.scalarData[:,0]
                    y=self.scalarData[:,col]
                    ax=plt.subplot2grid((self.nplots,10),(nplot, 0),colspan=10)
                    ax.plot(x,y)
                    ax.set_ylabel(j.text())
                    ax.axvline(x=time,c="red",linewidth=2,zorder=0, clip_on=False)
                    nplot+=1

        if not self.fieldFile is None :
            nameGroup="/%010d" % (self.step*self.fieldEvery)
            for i in self.ui.menuFields.actions() :
                if i.isChecked() :
                    nameData=nameGroup+"/"+i.text()
                
                    data=self.fieldFile.getNode(str(nameData))
                    
                    ax=plt.subplot2grid((self.nplots,10),(nplot, 0),colspan=9)
                    ax.set_ylabel(i.text())

                    if len(data.shape) == 1 :
                        x=np.array(range(data.shape[0]))/self.res_space
                        y=data
                        ax.plot(x,y)

                        if self.autoScale.isChecked() or self.lims[nplot]==None :
                            self.lims[nplot]=ax.get_ylim()

                        ax.set_ylim(self.lims[nplot])
                        ax.set_xlim(0,self.sim_length)
                    elif len(data.shape) == 2 :
                        data=np.array(data)
                        im=ax.imshow(data.T,extent=(0,self.sim_length[0],0,self.sim_length[1]), aspect='auto',origin='lower')
                        if self.autoScale.isChecked() or self.lims[nplot]==None :
                            self.lims[nplot]=(data.min(),data.max())
                        
                        im.set_clim(self.lims[nplot])
                        axcb=plt.subplot2grid((self.nplots,10),(nplot, 9)) 
                        cb=plt.colorbar(im, cax=axcb)

                    nplot+=1
                
        if not self.phaseFile is None :
            for i in self.ui.menuPhase_spaces.actions() :
                if i.isChecked() :
                    nameData=str(i.data().toString())
                    node=self.phaseFile.getNode(str(nameData))
                    data=node[self.step].T
                    
                    ax=plt.subplot2grid((self.nplots,10),(nplot, 0),colspan=9)

                    
                    im=ax.imshow(data,extent=node._v_parent._v_attrs.extents.reshape(4).tolist(), aspect='auto',origin='lower')

                    if self.autoScale.isChecked() or self.lims[nplot]==None :
                        self.lims[nplot]=(data.min(),data.max())
                        
                    im.set_clim(self.lims[nplot])

                    ax.set_ylabel(i.text())
                    axcb=plt.subplot2grid((self.nplots,10),(nplot, 9))
                    
                    cb=plt.colorbar(im, cax=axcb)


                    nplot+=1
                
                
        self.canvas.draw()
        if self.ui.actionSave_images.isChecked():
            plt.savefig('smilei-%06d.png' % self.step)

    def do_timer(self):
        self.ui.slider.setValue(self.step+1)
    

    def createActions(self):
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

        fname=os.path.join(self.dirName, "Fields.h5")
        if os.path.isfile(fname) :
            self.fieldSteps=[]
            f=tb.openFile(os.path.join(self.dirName, fname))
            
            self.res_time=f.root._v_attrs.res_time
            self.sim_length=f.root._v_attrs.sim_length
            self.fieldEvery=f.root._v_attrs.every
            
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
            self.ui.slider.setRange(0,len(self.fieldSteps)-1)

        self.ui.spinStep.setSuffix("/"+str(len(self.fieldSteps)-1))
        self.ui.spinStep.setMaximum(len(self.fieldSteps)-1)
        
        fname=os.path.join(self.dirName, "PhaseSpace.h5")
        if os.path.isfile(fname) :
            f=tb.openFile(os.path.join(self.dirName, fname))
            for phaseData in f.walkNodes("/", classname='Array'):
                namephase= phaseData._v_pathname + " " + phaseData._v_parent._v_attrs.species
                my_act= QtGui.QAction(namephase,self)
                my_act.setData(phaseData._v_pathname)
                my_act.setCheckable(True)                
                self.ui.menuPhase_spaces.addAction(my_act)
                my_act.triggered.connect(self.action_clicked)
            f.close()

        self.doPlots()
        
    def action_clicked(self):
        self.nplots=0
        menus=[self.ui.menuScalars,self.ui.menuFields,self.ui.menuPhase_spaces]
        for j in menus :
            for i in j.actions():
                if i.isChecked() :
                    self.nplots+=1
                    
                
        if self.nplots > 0:
            self.lims = [None] * self.nplots
            
            self.scalarData = None
            fname=os.path.join(self.dirName, "scalars.txt")
            if os.path.isfile(fname) :
                self.scalarData = np.loadtxt(fname)

            fname=os.path.join(self.dirName, "Fields.h5")
            if self.fieldFile != None:
                self.fieldFile.close()
                
            self.fieldFile = None
            if os.path.isfile(fname) :
                self.fieldFile=tb.openFile(os.path.join(self.dirName, fname))
                self.res_space=self.fieldFile.root._v_attrs.res_space[0]

            fname=os.path.join(self.dirName, "PhaseSpace.h5")
            if self.phaseFile != None:
                self.phaseFile.close()
                
            self.phaseFile = None
            if os.path.isfile(fname) :
                self.phaseFile=tb.openFile(os.path.join(self.dirName, fname))
            
            
            self.doPlots()



        
    def on_actionDir_triggered(self):
        print "on_actionDir_triggered", self.sender()
        dirName=QtGui.QFileDialog.getExistingDirectory(self,self.dirName, options=QFileDialog.ShowDirsOnly)
        if not dirName.isEmpty():
            self.dirName=str(dirName)
            self.createActions()

def main():

    app = QtGui.QApplication(sys.argv)
    ex = smileiQt()
    
    if sys.platform == "darwin":
        ex.raise_()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()  
