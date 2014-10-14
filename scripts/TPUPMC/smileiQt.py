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
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    from PyQt4 import uic

import tables as tb
import numpy as np

import matplotlib.pyplot as plt

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class smileiQtPlot(QWidget):

    def __init__(self,dirName):
        super(smileiQtPlot, self).__init__()
        

        uiFile=os.path.dirname(os.path.realpath(__file__))+'/smileiQtPlot.ui'
        self.ui=uic.loadUi(uiFile,self)

        self.setWindowFlags(Qt.Window)
        self.setWindowTitle(dirName)

        self.step=0
        
        self.ui.autoScale.toggled.connect(self.doPlots)
        self.ui.tabWidget.currentChanged.connect(self.doPlots)
        
        self.fig = plt.figure(dirName)
        self.canvas = FigureCanvas(self.fig)
        
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setFixedHeight(18)
        self.ui.plotLayout.addWidget(self.toolbar)
        self.ui.plotLayout.addWidget(self.canvas)

        self.scalarData = None
        self.fieldFile = None
        self.phaseFile = None

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
                my_button.stateChanged.connect(self.checkBox_clicked)
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
                        my_button.stateChanged.connect(self.checkBox_clicked)
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
                my_button.stateChanged.connect(self.checkBox_clicked)
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
                    
    def doPlots(self):
        if  self.ui.tabWidget.currentIndex()==1: return

        print "doPlots", self
        if len(self.fieldSteps) == 0 : return
        
        self.step %= len(self.fieldSteps)
        self.slider.setValue(self.step)
        
        self.ui.spinStep.setValue(self.step)
        time=self.step/self.res_time*self.fieldEvery
        self.fig.suptitle("Time: %.3f" % time)
        nplot=0
    
        if not self.scalarData is None :
            col=0
            for i in self.ui.scalars.findChildren(QCheckBox):
                col+=1
                if i.isChecked() :          
                    x=self.scalarData[:,0]
                    y=self.scalarData[:,col]
                    ax=self.fig.add_subplot(nplot,1)
                    ax.plot(x,y)
                    ax.set_ylabel(i.text())
                    ax.axvline(x=time,c="red",linewidth=2,zorder=0, clip_on=False)
                    nplot+=1

        if not self.fieldFile is None :
            nameGroup="/%010d" % (self.step*self.fieldEvery)
            for i in self.ui.fields.findChildren(QCheckBox):
                if i.isChecked() :
                    nameData=nameGroup+"/"+i.text()
                
                    data=self.fieldFile.getNode(str(nameData))
                    
                    ax=self.fig.add_subplot(nplot,1)
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
                        axcb= self.fig.add_subplot(nplot,2)
                        cb=plt.colorbar(im, cax=axcb)

                    nplot+=1
                    
        if not self.phaseFile is None :
            for i in self.ui.phase.findChildren(QCheckBox):
                if i.isChecked() :
                    nameData=i.text()
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
        if self.ui.saveImages.isChecked():
            plt.savefig('smilei-%06d.png' % self.step)
        
    def checkBox_clicked(self):
        self.nplots=0
        
        menus=[self.ui.layoutScalars,self.ui.layoutFields,self.ui.layoutPhase]
        for j in self.ui.Diags.findChildren(QCheckBox) :
            if j.isChecked() : self.nplots+=1
                    
        if self.nplots > 0:
            self.lims = [None] * self.nplots            
            self.doPlots()
       
    def closeEvent(self,event):
        self.fieldFile.close()
        self.phaseFile.close()
    

class smileiQt(QMainWindow):        

    next = pyqtSignal()
    previous = pyqtSignal()
    first = pyqtSignal()
    last = pyqtSignal()
    timer=QTimer()
    
    plots=[]
    def __init__(self, args):
        super(smileiQt, self).__init__()
                
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

        self.fieldFile=None
        self.phaseFile=None

        self.timer.setInterval(500)

        if sys.platform == "darwin":
            self.raise_()
        self.show()
        
    def addDir(self,name):
        plot=smileiQtPlot(name)
        self.plots.append(plot)
        self.timer.timeout.connect(plot.next)
        self.ui.next.released.connect(plot.next)
        self.ui.previous.released.connect(plot.previous)
        self.ui.first.released.connect(plot.first)
        self.ui.last.released.connect(plot.last)
             
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
