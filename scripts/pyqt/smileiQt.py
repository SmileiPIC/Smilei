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

import os
import tables
import numpy as np

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class smileiQt(QtGui.QMainWindow):
    
    def __init__(self):
        super(smileiQt, self).__init__()
            
        self.ui=uic.loadUi(os.path.dirname(os.path.realpath(__file__))+'/smileiQt.ui',self)
        self.ui.timeStep.currentIndexChanged.connect(self.on_draw)
        self.ui.field1.currentIndexChanged.connect(self.on_draw)
        self.ui.field2.currentIndexChanged.connect(self.on_draw)
        
        validator=QtGui.QDoubleValidator()
        self.ui.mini.setValidator(validator)
        self.ui.maxi.setValidator(validator)
        
        self.ui.mini.editingFinished.connect(self.on_draw)
        self.ui.maxi.editingFinished.connect(self.on_draw)
        
        self.ui.actionOpen_File.triggered.connect(self.on_file_change)
        self.logBox.stateChanged.connect(self.on_draw)

        self.h5data1 = []
        self.h5data2 = []
        self.filename = "" 
        
        self.fig1 = Figure()
        self.canvas1 = FigureCanvas(self.fig1)
        self.canvas1.mpl_connect('motion_notify_event', self.on_movement)
        self.ui.grid1.addWidget(self.canvas1,0,0)

        self.fig2 = Figure()
        self.canvas2 = FigureCanvas(self.fig2)
        self.canvas2.mpl_connect('motion_notify_event', self.on_movement)
        self.ui.grid2.addWidget(self.canvas2,0,0)

        self.load_settings()
        self.update_files()
        self.show()

    def load_settings(self):
        settings=QtCore.QSettings("smilePy","");
        settings.beginGroup("Preferences");
        self.filename=str(settings.value("filename","").toString());
        settings.endGroup();
        self.update_files()

    def save_settings(self):
        settings=QtCore.QSettings("smilePy","");
        settings.beginGroup("Preferences");
        settings.setValue("filename",self.filename);
        settings.endGroup();

    def on_movement(self, event):
        if not (event.inaxes is None) :
            zval1=self.h5data1[int(event.ydata),int(event.xdata)]
            zval2=self.h5data2[int(event.ydata),int(event.xdata)]
            msg = "(%d,%d) %.3f %.3f" % (int(event.xdata), int(event.ydata), zval1, zval2)
            self.statusBar().showMessage(msg)

    def on_file_change(self):
        filename = QtGui.QFileDialog.getOpenFileName(self,"Open File", self.filename, "HDF5 Files (*.h5)")
        if os.path.isfile(filename) :
            self.filename=str(filename)
            self.update_files()
    
    def update_files (self):
        print self.filename

        self.ui.timeStep.currentIndexChanged.disconnect(self.on_draw)
        self.ui.field1.currentIndexChanged.disconnect(self.on_draw)
        self.ui.field2.currentIndexChanged.disconnect(self.on_draw)

        if os.path.isfile(self.filename) : 
            self.setWindowTitle("Smilei "+self.filename);
            self.save_settings()
            self.field1.clear()
            self.field2.clear()
            self.timeStep.clear()
            
            fieldlist = []
            h5file=tables.open_file(self.filename, mode = "r")
            
            first=True
            for time in h5file.root:
                self.timeStep.addItem(time._v_name)
                if first:
                    first=False
                    for field in time :
                        self.field1.addItem(str(field._v_name))
                        self.field2.addItem(str(field._v_name))
                
            h5file.close()
        
            self.ui.timeStep.currentIndexChanged.connect(self.on_draw)
            self.ui.field1.currentIndexChanged.connect(self.on_draw)
            self.ui.field2.currentIndexChanged.connect(self.on_draw)
    
            self.on_draw()

    def on_draw(self):
        """display dir
        """        
        h5file=tables.open_file(self.filename, mode = "r")
        
        name1=str("/"+self.timeStep.currentText()+"/"+self.field1.currentText())
        name2=str("/"+self.timeStep.currentText()+"/"+self.field2.currentText())
        self.h5data1 = h5file.getNode(name1).read()
        self.h5data2 = h5file.getNode(name2).read()

        h5file.close()
        
        print name1,self.h5data1.shape
        print name2,self.h5data2.shape
        
        self.fig1.clear()
        self.fig2.clear()
        self.axes1 = self.fig1.add_subplot(111)
        self.axes2 = self.fig2.add_subplot(111)
        
        mini=self.mini.text().toDouble()
        maxi=self.maxi.text().toDouble()
        if mini[1] and maxi[1] :
            if self.ui.logBox.isChecked() and mini[0]>0 and maxi[0]>0:
                self.img1 = self.axes1.imshow(self.h5data1,norm=LogNorm(vmin=mini[0],vmax=maxi[0]))
                self.img2 = self.axes2.imshow(self.h5data2,norm=LogNorm(vmin=mini[0],vmax=maxi[0]))
            else:
                self.img1 = self.axes1.imshow(self.h5data1,vmin=mini[0],vmax=maxi[0])
                self.img2 = self.axes2.imshow(self.h5data2,vmin=mini[0],vmax=maxi[0])
        else :
            self.img1 = self.axes1.imshow(self.h5data1)    
            self.img2 = self.axes2.imshow(self.h5data2)    

        self.fig1.colorbar(self.img1)
        self.fig2.colorbar(self.img2)

        self.canvas1.draw()
        self.canvas2.draw()

        
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = smileiQt()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()  