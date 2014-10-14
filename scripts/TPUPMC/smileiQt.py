#!/usr/bin/env python 
"""
Plot fields of smilei simulaition
"""
import sys, os, random

from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends import qt_compat

if qt_compat.QT_API == qt_compat.QT_API_PYSIDE:
    from PySide.QtCore import *
    from PySide.QtGui import *
    from PyQt4 import uic
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    from PyQt4 import uic

import tables as tb
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class smileiQtPlot(QWidget):

    def __init__(self,parent,dirName):
        super(smileiQtPlot, self).__init__()    
        self.setParent(parent)
        uiFile=os.path.dirname(os.path.realpath(__file__))+'/smileiQtPlot.ui'
        self.ui=uic.loadUi(uiFile,self)
        
        self.parent=parent

        self.parent.timer.timeout.connect(self.next)
        self.parent.ui.next.released.connect(self.next)
        self.parent.ui.previous.released.connect(self.previous)
        self.parent.ui.first.released.connect(self.first)
        self.parent.ui.last.released.connect(self.last)

        self.setWindowFlags(Qt.Window)
        self.setWindowTitle(dirName)

        self.step=0
        
        self.ui.autoScale.toggled.connect(self.doPlots)
        self.ui.tabWidget.currentChanged.connect(self.changeTab)
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setFixedHeight(18)
        self.ui.plotLayout.addWidget(self.canvas)
        self.ui.plotLayout.addWidget(self.toolbar)

        fname=os.path.join(dirName, "scalars.txt")
        if os.path.isfile(fname) :
            self.scalarData = np.loadtxt(fname)
            names=[]
            for line in open(fname):
                li=line.strip()
                if li.startswith("#"):
                    list=line.split()
                    names.append(list[-1])
       
            scalars_names=names[1:-2]
            for i in range(len(scalars_names)):
                my_button= QCheckBox(scalars_names[i])
                self.ui.layoutScalars.addWidget(my_button)
                
            self.ui.layoutScalars.addStretch()
        else :
            self.deleteLater()

        self.fieldSteps=[]
        fname=os.path.join(dirName, "Fields.h5")
        if os.path.isfile(fname) :

            self.fieldFile=tb.openFile(fname)
            self.res_space=self.fieldFile.root._v_attrs.res_space[0]
            self.res_time=self.fieldFile.root._v_attrs.res_time
            self.sim_length=self.fieldFile.root._v_attrs.sim_length
            self.fieldEvery=self.fieldFile.root._v_attrs.every
            
            first=True
            for group in self.fieldFile.list_nodes("/", classname='Group'):
                self.fieldSteps.append(group._v_name)
                if first:
                    first=False
                    for array in group:
                        my_button= QCheckBox(array._v_name)
                        self.ui.layoutFields.addWidget(my_button)

            self.ui.layoutFields.addStretch()
            self.ui.slider.setRange(0,len(self.fieldSteps)-1)
        else :
            self.deleteLater()

        self.ui.spinStep.setSuffix("/"+str(len(self.fieldSteps)-1))
        self.ui.spinStep.setMaximum(len(self.fieldSteps)-1)
        
        fname=os.path.join(dirName, "PhaseSpace.h5")
        if os.path.isfile(fname) :
            self.phaseFile=tb.openFile(fname)
            for phaseData in self.phaseFile.walkNodes("/", classname='Array'):
                my_button= QCheckBox(phaseData._v_pathname)
                self.ui.layoutPhase.addWidget(my_button)
            self.ui.layoutPhase.addStretch()
        else :
            self.deleteLater()

        if sys.platform == "darwin":
            self.raise_()
        self.show()

    def next(self):
        self.ui.slider.setValue(self.step+1)
    
    def previous(self):
        self.ui.slider.setValue(self.step-1)
    
    def first(self):
        self.ui.slider.setValue(0)
    
    def last(self):
        self.ui.slider.setValue(len(self.fieldSteps)-1)
               
    def on_slider_valueChanged(self,step):
        self.step=step
        self.doPlots()
    
    @pyqtSignature("int")
    def on_spinStep_valueChanged(self,my_step):
        self.ui.slider.setValue(my_step)

    @pyqtSignature("int")
    def changeTab(self, tabNum):
        if  self.ui.tabWidget.currentIndex()==0: 

            self.nplots=0
            menus=[self.ui.layoutScalars,self.ui.layoutFields,self.ui.layoutPhase]
            for j in self.ui.Diags.findChildren(QCheckBox) :
                if j.isChecked() : self.nplots+=1
                    
            if self.nplots > 0:
                self.lims = [None] * self.nplots            
                self.doPlots()
    
    def doPlots(self):

        if len(self.fieldSteps) == 0 : return
        
        self.step %= len(self.fieldSteps)
        self.slider.setValue(self.step)
        
        self.ui.spinStep.setValue(self.step)
        time=self.step/self.res_time*self.fieldEvery
        self.fig.suptitle("Time: %.3f" % time)
        nplot=0
        self.fig.clear()
        col=0
        for i in self.ui.scalars.findChildren(QCheckBox):
            col+=1
            if i.isChecked() :          
                x=self.scalarData[:,0]
                y=self.scalarData[:,col]
                ax=self.fig.add_subplot(self.nplots,1,nplot+1)
                ax.plot(x,y)
                ax.set_ylabel(i.text())
                ax.axvline(x=time,c="red",linewidth=2,zorder=0, clip_on=False)
                nplot+=1

        nameGroup="/%010d" % (self.step*self.fieldEvery)
        for i in self.ui.fields.findChildren(QCheckBox):
            if i.isChecked() :
                nameData=nameGroup+"/"+i.text()
            
                data=self.fieldFile.getNode(str(nameData))
                
                ax=self.fig.add_subplot(self.nplots,1,nplot+1)
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
                    divider = make_axes_locatable(ax)
                    cax=divider.append_axes("right", "5%", pad="3%")
                    cb=plt.colorbar(im, cax=cax)

                nplot+=1
                    
        for i in self.ui.phase.findChildren(QCheckBox):
            if i.isChecked() :
                nameData=i.text()
                node=self.phaseFile.getNode(str(nameData))
                data=node[self.step].T
                
                ax=self.fig.add_subplot(self.nplots,1,nplot+1)
                im=ax.imshow(data,extent=node._v_parent._v_attrs.extents.reshape(4).tolist(), aspect='auto',origin='lower')
                if self.autoScale.isChecked() or self.lims[nplot]==None :
                    self.lims[nplot]=(data.min(),data.max())                    
                im.set_clim(self.lims[nplot])
                ax.set_ylabel(i.text())
                

                divider = make_axes_locatable(ax)
                cax = divider.new_horizontal(size="2%", pad=0.05)
                self.fig.add_axes(cax)
                
                cb=plt.colorbar(im, cax=cax)


                nplot+=1
                
                
        self.canvas.draw()
        if self.ui.saveImages.isChecked():
            plt.savefig('smilei-%06d.png' % self.step)
               
    def closeEvent(self,event):
        self.fieldFile.close()
        self.phaseFile.close()
        self.parent.plots.remove(self)
        self.deleteLater()


class smileiQt(QMainWindow):        

    next = pyqtSignal()
    previous = pyqtSignal()
    first = pyqtSignal()
    last = pyqtSignal()
    timer=QTimer()
    
    plots=[]
    def __init__(self, args):
        super(smileiQt, self).__init__()
                
        self.setWindowFlags(self.windowFlags() | Qt.Tool | Qt.Dialog)

        uiFile=os.path.dirname(os.path.realpath(__file__))+'/smileiQt.ui'
        self.ui=uic.loadUi(uiFile,self)

        self.ui.dir.setIcon(self.ui.style().standardIcon(QStyle.SP_DirOpenIcon))
        self.ui.playStop.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaPlay))
        self.ui.previous.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaSeekBackward))
        self.ui.next.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaSeekForward))
        self.ui.first.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaSkipBackward))
        self.ui.last.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaSkipForward))        
        
        for i in args : 
            self.addDir(i)

        self.timer.setInterval(500)
        self.timer.timeout.connect(self.timeout)

        if sys.platform == "darwin":
            self.raise_()
        self.show()
        
    def timeout(self):
        if len(self.plots) == 0:
            self.ui.playStop.setChecked(False)
            
    def addDir(self,name):
        self.plots.append(smileiQtPlot(self,name))
             
    def on_playStop_toggled(self):
        if self.playStop.isChecked() :
            self.ui.playStop.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaStop))
            self.timer.start()
        else:
            self.ui.playStop.setIcon(self.ui.style().standardIcon(QStyle.SP_MediaPlay))
            self.timer.stop()

    def on_dir_released(self):
        dirName=QFileDialog.getExistingDirectory(self,options=QFileDialog.ShowDirsOnly)
        if not dirName.isEmpty():
            self.addDir(str(dirName))
        

def main():

    app = QApplication(sys.argv)
    args = ["."] if len(sys.argv) == 1 else sys.argv[1:]

    smilei=smileiQt(args)

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()  
