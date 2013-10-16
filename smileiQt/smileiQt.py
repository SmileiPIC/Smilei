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
        self.ui=uic.loadUi('smileiQt.ui',self)
        self.ui.comboField.currentIndexChanged.connect(self.on_draw)
        self.ui.comboPart.currentIndexChanged.connect(self.on_part_change)
        self.ui.comboPart2.currentIndexChanged.connect(self.on_draw)
        
        validator=QtGui.QDoubleValidator()
        self.ui.mini.setValidator(validator)
        self.ui.maxi.setValidator(validator)
        
        self.ui.mini.editingFinished.connect(self.on_draw)
        self.ui.maxi.editingFinished.connect(self.on_draw)
        
        self.ui.actionChange_Dir.triggered.connect(self.on_dir_cange)
        self.logBox.stateChanged.connect(self.on_draw)

        self.h5data = []

        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.ui.plot)
        self.canvas.mpl_connect('motion_notify_event', self.on_movement)

        self.ui.grid.addWidget(self.canvas,0,0)
        self.load_settings()
        self.update_files()
        self.show()

    def load_settings(self):
        settings=QtCore.QSettings("smilePy","");
        settings.beginGroup("Preferences");
        self.dirname=QtCore.QDir(settings.value("dirname",".").toString());
        settings.endGroup();

    def save_settings(self):
        settings=QtCore.QSettings("smilePy","");
        settings.beginGroup("Preferences");
        settings.setValue("dirname",self.dirname.path());
        settings.endGroup();

    def on_movement(self, event):
        if not (event.inaxes is None) :
            zval=self.h5data[int(event.ydata),int(event.xdata)]
            msg = "(%d,%d) %.3f" % (int(event.xdata), int(event.ydata), zval)
            self.statusBar().showMessage(msg)

    def on_dir_cange(self):
        dirname = QtGui.QFileDialog.getExistingDirectory(self,"Open Directory",self.dirname.absolutePath(),QtGui.QFileDialog.ShowDirsOnly | QtGui.QFileDialog.DontResolveSymlinks);
        self.update_files(QtCore.QDir(dirname))
    
    def update_files (self, dirname=QtCore.QDir('..')):
        if dirname.exists(): 
            self.dirname=dirname
            self.save_settings()
            self.comboField.clear()
            self.comboPart.clear()
        
            directory = QtCore.QDir(self.dirname);
            for file in directory.entryList(QtCore.QStringList("*.h5")):
                h5file=tables.open_file(str(directory.absoluteFilePath(file)), mode = "r")
                if h5file.root.__contains__("Field"):
                    self.comboField.addItem(file,directory.absoluteFilePath(file))
                else:                    
                    self.comboPart.addItem(file,directory.absoluteFilePath(file))
                h5file.close()
        
            self.on_draw()

    def on_draw(self):
        """display dir
        """        
        file = str(self.comboField.itemData(self.comboField.currentIndex()).toString())
        if  os.path.isfile(file):        
            self.fig.clear()
            self.axes = self.fig.add_subplot(111)
            self.axes.set_title(os.path.basename(file))
            h5file=tables.open_file(file, mode = "r")
            self.h5data = h5file.root.Field.read()
            h5file.close()
            mini=self.mini.text().toDouble()
            maxi=self.maxi.text().toDouble()
            if mini[1] and maxi[1] :
                if self.ui.logBox.isChecked() and mini[0]>0 and maxi[0]>0:
                    self.img = self.axes.imshow(self.h5data,norm=LogNorm(vmin=mini[0],vmax=maxi[0]))
                else:
                    self.img = self.axes.imshow(self.h5data,vmin=mini[0],vmax=maxi[0])
            else :
                self.img = self.axes.imshow(self.h5data)    

            cbar = self.fig.colorbar(self.img)
        name = str(self.comboPart2.itemData(self.comboPart2.currentIndex()).toString())
        if name:
            file=str(self.comboPart.itemData(self.comboPart.currentIndex()).toString())
            if  os.path.isfile(file):
                h5file=tables.open_file(file, mode = "r")
                print file+" / "+name
                method=getattr(h5file.root,name)
#                 partData=method.read()
#                 np.savetxt('test.txt',partData)
                h5file.close()
        self.canvas.draw()

    def on_part_change(self, number = 0):
        file = str(self.comboPart.itemData(number).toString())
        h5file=tables.open_file(file, mode = "r")
        self.comboPart2.clear()
        self.comboPart2.addItem('None','')
        for node in h5file.root:
            self.comboPart2.addItem(node.name,node.name)
        h5file.close()


        
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = smileiQt()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()  