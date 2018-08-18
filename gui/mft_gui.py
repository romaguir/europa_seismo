#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
gui for dispersion measurement with the 
multiple frequency technique (Dziewonski, 1969)

Some of this code was copied and pasted from 
the instaseis gui written by Lion Krischer
"""

import os
import sys
import imp
import obspy
import inspect
import instaseis
import pyqtgraph as pg
import matplotlib.pyplot as plt
from glob import iglob
from PyQt4 import QtGui, QtCore
from PyQt4 import uic

class Window(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = layout.Ui_MainWindow()
        self.ui.setupUi(self)

#borrowed from instaseis gui... still learning how this works
def use_ui_layout():
    DIR = os.path.dirname(os.path.abspath(__file__)) 
    print DIR
    ui_file = DIR+'/layout.ui'
    py_file = DIR+'/layout.py'
    with open(py_file, 'w') as open_file:
        uic.compileUi(ui_file, open_file)

    import_name = os.path.splitext(os.path.basename(py_file))[0]
    print "import name", import_name
    globals()[import_name] = imp.load_source(import_name, py_file)

def launch():
    
    use_ui_layout()
    app = QtGui.QApplication(sys.argv, QtGui.QApplication.GuiClient)
    window = Window()
    window.show()
    app.installEventFilter(window)
    window.raise_()
    os._exit(app.exec_()) 
