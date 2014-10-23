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
    scalarDict=dict()
    fieldDict=dict()
    phaseDict=dict()
    fieldFile=None
    phaseFile=None
    nplots=0
    ax={}
    dirName=None
    someCheckBoxChanged=False
    
    def __init__(self,parent,dirName):
        super(smileiQtPlot, self).__init__()   
        self.setParent(parent)
        uiFile=os.path.dirname(os.path.realpath(__file__))+'/smileiQtPlot.ui'
        self.ui=uic.loadUi(uiFile,self)
        
        self.dirName=dirName
        self.parent=parent

        self.parent.timer.timeout.connect(self.next)
        self.parent.ui.next.released.connect(self.next)
        self.parent.ui.previous.released.connect(self.previous)
        self.parent.ui.first.released.connect(self.first)
        self.parent.ui.last.released.connect(self.last)

        self.setWindowFlags(Qt.Window)
        self.setWindowTitle(dirName)

        self.step=0
        
        self.ui.tabWidget.currentChanged.connect(self.changeTab)
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.canvas.setFocus()
        self.canvas.mpl_connect('motion_notify_event', self.on_movement)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)

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
                my_button.stateChanged.connect(self.checkBoxChanged)

                self.ui.layoutScalars.addWidget(my_button)
                
            self.ui.layoutScalars.addStretch()
        else :
            print "Problem reading ",fname
#             self.deleteLater()

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
                        my_button.stateChanged.connect(self.checkBoxChanged)
                        self.ui.layoutFields.addWidget(my_button)

            self.ui.layoutFields.addStretch()
            self.ui.slider.setRange(0,len(self.fieldSteps)-1)
        else :
            print "Problem reading ",fname
#             self.deleteLater()

        self.ui.spinStep.setSuffix("/"+str(len(self.fieldSteps)-1))
        self.ui.spinStep.setMaximum(len(self.fieldSteps)-1)
        
        fname=os.path.join(dirName, "PhaseSpace.h5")
        if os.path.isfile(fname) :
            self.phaseFile=tb.openFile(fname)
            for phaseData in self.phaseFile.walkNodes("/", classname='Array'):
                my_button= QCheckBox(phaseData._v_pathname)
                my_button.stateChanged.connect(self.checkBoxChanged)
                self.ui.layoutPhase.addWidget(my_button)
            self.ui.layoutPhase.addStretch()
        else :
            print "Problem reading ",fname
#             self.deleteLater()

        self.load_settings()

        if sys.platform == "darwin":
            self.raise_()
        self.show()

    def checkBoxChanged(self):
        self.someCheckBoxChanged=True

    def load_settings(self):
        settings=QSettings(QFileInfo(__file__).fileName(),"")
        settings.beginGroup(QDir(self.dirName).dirName())
        frames=[self.ui.scalars, self.ui.fields, self.ui.phase]
        for frame in [self.ui.scalars, self.ui.fields, self.ui.phase] :
            settings.beginGroup(frame.objectName())            
            for chkbox in frame.findChildren(QCheckBox):
                chkbox.setChecked(settings.value(chkbox.text()).toBool())
            settings.endGroup()
        settings.endGroup()
        self.ui.tabWidget.setCurrentIndex(0)

    def save_settings(self):
        settings=QSettings(QFileInfo(__file__).fileName(),"")
        settings.beginGroup(QDir(self.dirName).dirName())
        frames=[self.ui.scalars, self.ui.fields, self.ui.phase]
        for frame in [self.ui.scalars, self.ui.fields, self.ui.phase] :
            settings.beginGroup(frame.objectName())            
            for chkbox in frame.findChildren(QCheckBox):
                settings.setValue(chkbox.text(),chkbox.isChecked())
            settings.endGroup()
        settings.endGroup()

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
        if self.ui.tabWidget.currentIndex()==0:
            self.doPlots()


    def preparePlots(self):
        self.someCheckBoxChanged=False
        
        self.scalarDict=dict()
        self.fieldDict=dict()
        self.phaseDict=dict()

        self.nplots=0
        frames=[self.ui.scalars, self.ui.fields, self.ui.phase]
        for frame in [self.ui.scalars, self.ui.fields, self.ui.phase] :
            for chkbox in frame.findChildren(QCheckBox):
                if chkbox.isChecked() :
                    self.nplots+=1

        if self.nplots > 0:
            self.fig.clear()
            self.ax={}
              
            plot=0
            col=0
            for i in self.ui.scalars.findChildren(QCheckBox):
                col+=1
                if i.isChecked() :
                    name=str(i.text())
                    x=self.scalarData[:,0]
                    y=self.scalarData[:,col]
                    self.scalarDict[name]=(x,y)
                    ax=self.fig.add_subplot(self.nplots,1,plot+1)
                    ax.xaxis.grid(True)
                    ax.yaxis.grid(True)
                    ax.plot(x,y)
                    ax.set_xlim(x.min(),x.max())

                    ax.set_ylabel(name)
                    ax.axvline(x=0,c="red",linewidth=2,zorder=0, clip_on=False)
                    self.ax[name]=ax
                    plot+=1

            for i in self.ui.fields.findChildren(QCheckBox):
                if i.isChecked() :
                    data=[]
                    name=str(i.text())
                    for d in self.fieldFile.root:
                        data.append(d._f_getChild(name))
                    data=np.array(data)
                    self.fieldDict[name]=data
                    ax=self.fig.add_subplot(self.nplots,1,plot+1)
                    ax.xaxis.grid(True)
                    ax.yaxis.grid(True)
                    
                    if len(data.shape) == 2 :
                        ax.set_xlim(0,self.sim_length)
                        ax.set_ylim(data.min(), data.max())
                        ax.set_ylabel(name)
                        x=np.array(range(data.shape[1]))/self.res_space
                        y=data[0]
                        ax.plot(x,y)
                        self.ax[name]=ax
                    elif len(data.shape) == 3 :
                        divider = make_axes_locatable(ax)
                        cax = divider.new_horizontal(size="2%", pad=0.05)
                        self.fig.add_axes(cax)
                        ax.set_ylabel(name)

                        im=ax.imshow([[0]],extent=(0,self.sim_length[0],0,self.sim_length[1]), aspect='auto',origin='lower')
                        im.set_clim(data.min(),data.max())
                        cb=plt.colorbar(im, cax=cax)
                        self.ax[name]=ax

                    plot+=1

            for i in self.ui.phase.findChildren(QCheckBox):
                if i.isChecked() :
                    data=dict()
                    name=str(i.text())
                    node=self.phaseFile.getNode(name)
                    data['extent']=node._v_parent._v_attrs.extents.reshape(4).tolist()
                    data['data']=node

                    self.phaseDict[name]=data
                    ax=self.fig.add_subplot(self.nplots,1,plot+1)
                    ax.xaxis.grid(True)
                    ax.yaxis.grid(True)
                    ax.set_ylabel(name)
                    divider = make_axes_locatable(ax)
                    cax = divider.new_horizontal(size="2%", pad=0.05)
                    self.fig.add_axes(cax)

                    im=ax.imshow([[0]],extent=self.phaseDict[name]['extent'],aspect='auto',origin='lower')
                    cb=plt.colorbar(im, cax=cax)


                    self.ax[name]=ax
                    plot+=1
                
            self.doPlots()
    
    def on_movement(self, event):
        if not (event.inaxes is None) :
            msg = "%G %G" % (event.xdata, event.ydata)
            self.ui.pos.setText(msg)

    def on_click(self,event):
        if event.button==3 :
            self.ui.logger.moveCursor (QTextCursor.End);
            self.ui.logger.insertPlainText('\n%g %g'%(event.xdata, event.ydata));
            self.ui.logger.moveCursor (QTextCursor.End);

    def on_key_press(self,event):
        if event.key == 'a':
            self.doPlots(True)
        
    def doPlots(self, autoscale=False):
            
        if len(self.fieldSteps) == 0 : return

        if self.someCheckBoxChanged==True:
            self.preparePlots()
        
        self.step %= len(self.fieldSteps)
        self.slider.setValue(self.step)
        
        self.ui.spinStep.setValue(self.step)
        time=self.step/self.res_time*self.fieldEvery
        
        for name in self.scalarDict:
            self.ax[name].lines[-1].set_xdata(time)
           
        for name in self.fieldDict:
            data=self.fieldDict[name][self.step]
            if len(data.shape) == 1 :
                self.ax[name].lines[-1].set_ydata(data)
                if autoscale:
                    self.ax[name].set_ylim(data.min(),data.max())
            elif len(data.shape) == 2 :
                im=self.ax[name].images[-1]
                im.set_data(data.T)
                im.set_clim(data.min(),data.max())

                    
        for name in self.phaseDict:
            data=self.phaseDict[name]['data'][self.step].T
            im=self.ax[name].images[-1]
            im.set_data(data)
            im.set_clim(data.min(),data.max())
                
                
        self.fig.suptitle("Time: %.3f" % time)       
        self.canvas.draw()
        if self.ui.saveImages.isChecked():
            self.fig.savefig(self.dirName+'-%06d.png' % self.step)
               
    def closeEvent(self,event):
        print "Closing window ",self.windowTitle()
        self.save_settings()
        if self.fieldFile is not None : self.fieldFile.close()
        if self.phaseFile is not None : self.phaseFile.close()
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

        self.timer.setInterval(100)
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
    
    
    def closeEvent(self,event):
        for plot in self.plots:
            plot.deleteLater()
        self.deleteLater()


def main():

    app = QApplication(sys.argv)
    args = ["."] if len(sys.argv) == 1 else sys.argv[1:]

    smilei=smileiQt(args)

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()  
