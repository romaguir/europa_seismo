# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/romaguir/Tools/europa_seismo/gui/layout.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1065, 869)
        MainWindow.setInputMethodHints(QtCore.Qt.ImhFormattedNumbersOnly)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.open_sac_button = QtGui.QPushButton(self.centralwidget)
        self.open_sac_button.setObjectName(_fromUtf8("open_sac_button"))
        self.horizontalLayout_2.addWidget(self.open_sac_button)
        self.pushButton = QtGui.QPushButton(self.centralwidget)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.horizontalLayout_2.addWidget(self.pushButton)
        self.open_instaseis_button = QtGui.QPushButton(self.centralwidget)
        self.open_instaseis_button.setObjectName(_fromUtf8("open_instaseis_button"))
        self.horizontalLayout_2.addWidget(self.open_instaseis_button)
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_2.addWidget(self.label_2)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.horizontalLayout_2.addWidget(self.line)
        self.horizontalSlider = QtGui.QSlider(self.centralwidget)
        self.horizontalSlider.setMaximum(100)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName(_fromUtf8("horizontalSlider"))
        self.horizontalLayout_2.addWidget(self.horizontalSlider)
        self.line_2 = QtGui.QFrame(self.centralwidget)
        self.line_2.setFrameShape(QtGui.QFrame.VLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.horizontalLayout_2.addWidget(self.line_2)
        self.integrate_button = QtGui.QPushButton(self.centralwidget)
        self.integrate_button.setObjectName(_fromUtf8("integrate_button"))
        self.horizontalLayout_2.addWidget(self.integrate_button)
        self.differentiate_button = QtGui.QPushButton(self.centralwidget)
        self.differentiate_button.setObjectName(_fromUtf8("differentiate_button"))
        self.horizontalLayout_2.addWidget(self.differentiate_button)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.line_3 = QtGui.QFrame(self.centralwidget)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.verticalLayout.addWidget(self.line_3)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_3 = QtGui.QGridLayout()
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.label_17 = QtGui.QLabel(self.centralwidget)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.horizontalLayout_3.addWidget(self.label_17)
        self.corners = QtGui.QComboBox(self.centralwidget)
        self.corners.setObjectName(_fromUtf8("corners"))
        self.corners.addItem(_fromUtf8(""))
        self.corners.addItem(_fromUtf8(""))
        self.corners.addItem(_fromUtf8(""))
        self.corners.addItem(_fromUtf8(""))
        self.horizontalLayout_3.addWidget(self.corners)
        self.gridLayout_3.addLayout(self.horizontalLayout_3, 3, 2, 1, 1)
        self.prefilter_min_period = QtGui.QDoubleSpinBox(self.centralwidget)
        self.prefilter_min_period.setMaximum(500.0)
        self.prefilter_min_period.setProperty("value", 1.0)
        self.prefilter_min_period.setObjectName(_fromUtf8("prefilter_min_period"))
        self.gridLayout_3.addWidget(self.prefilter_min_period, 2, 1, 1, 1)
        self.label_15 = QtGui.QLabel(self.centralwidget)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.gridLayout_3.addWidget(self.label_15, 3, 0, 1, 1)
        self.label_13 = QtGui.QLabel(self.centralwidget)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_3.addWidget(self.label_13, 1, 0, 1, 1)
        self.label_14 = QtGui.QLabel(self.centralwidget)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_3.addWidget(self.label_14, 2, 0, 1, 1)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.planet = QtGui.QComboBox(self.centralwidget)
        self.planet.setObjectName(_fromUtf8("planet"))
        self.planet.addItem(_fromUtf8(""))
        self.planet.addItem(_fromUtf8(""))
        self.planet.addItem(_fromUtf8(""))
        self.verticalLayout_2.addWidget(self.planet)
        self.map = Qt4MplCanvas(self.centralwidget)
        self.map.setObjectName(_fromUtf8("map"))
        self.verticalLayout_2.addWidget(self.map)
        self.gridLayout_3.addLayout(self.verticalLayout_2, 5, 0, 1, 3)
        self.prefilter_max_period = QtGui.QDoubleSpinBox(self.centralwidget)
        self.prefilter_max_period.setMinimum(1.0)
        self.prefilter_max_period.setMaximum(1000.0)
        self.prefilter_max_period.setProperty("value", 200.0)
        self.prefilter_max_period.setObjectName(_fromUtf8("prefilter_max_period"))
        self.gridLayout_3.addWidget(self.prefilter_max_period, 3, 1, 1, 1)
        self.line_4 = QtGui.QFrame(self.centralwidget)
        self.line_4.setFrameShape(QtGui.QFrame.HLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_3.addWidget(self.line_4, 0, 0, 1, 1)
        self.zerophase = QtGui.QCheckBox(self.centralwidget)
        self.zerophase.setObjectName(_fromUtf8("zerophase"))
        self.gridLayout_3.addWidget(self.zerophase, 2, 2, 1, 1)
        self.line_5 = QtGui.QFrame(self.centralwidget)
        self.line_5.setFrameShape(QtGui.QFrame.HLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout_3.addWidget(self.line_5, 4, 0, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_3, 1, 0, 1, 1)
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.evlo = QtGui.QDoubleSpinBox(self.centralwidget)
        self.evlo.setMinimum(-180.0)
        self.evlo.setMaximum(180.0)
        self.evlo.setObjectName(_fromUtf8("evlo"))
        self.gridLayout_2.addWidget(self.evlo, 2, 3, 1, 1)
        self.m_tp = ScientificDoubleSpinBox(self.centralwidget)
        self.m_tp.setObjectName(_fromUtf8("m_tp"))
        self.gridLayout_2.addWidget(self.m_tp, 5, 1, 1, 1)
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_2.addWidget(self.label_4, 4, 0, 1, 1)
        self.label_10 = QtGui.QLabel(self.centralwidget)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_2.addWidget(self.label_10, 2, 2, 1, 1)
        self.stlo = QtGui.QDoubleSpinBox(self.centralwidget)
        self.stlo.setMinimum(-180.0)
        self.stlo.setMaximum(180.0)
        self.stlo.setProperty("value", 30.0)
        self.stlo.setObjectName(_fromUtf8("stlo"))
        self.gridLayout_2.addWidget(self.stlo, 4, 3, 1, 1)
        self.label_5 = QtGui.QLabel(self.centralwidget)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_2.addWidget(self.label_5, 3, 0, 1, 1)
        self.label_6 = QtGui.QLabel(self.centralwidget)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_2.addWidget(self.label_6, 2, 0, 1, 1)
        self.label_8 = QtGui.QLabel(self.centralwidget)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_2.addWidget(self.label_8, 0, 0, 1, 1)
        self.label_7 = QtGui.QLabel(self.centralwidget)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_2.addWidget(self.label_7, 1, 0, 1, 1)
        self.evdp = QtGui.QDoubleSpinBox(self.centralwidget)
        self.evdp.setMaximum(700.0)
        self.evdp.setObjectName(_fromUtf8("evdp"))
        self.gridLayout_2.addWidget(self.evdp, 1, 3, 1, 1)
        self.m_rr = ScientificDoubleSpinBox(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.m_rr.sizePolicy().hasHeightForWidth())
        self.m_rr.setSizePolicy(sizePolicy)
        self.m_rr.setObjectName(_fromUtf8("m_rr"))
        self.gridLayout_2.addWidget(self.m_rr, 0, 1, 1, 1)
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_2.addWidget(self.label, 3, 2, 1, 1)
        self.m_tt = ScientificDoubleSpinBox(self.centralwidget)
        self.m_tt.setObjectName(_fromUtf8("m_tt"))
        self.gridLayout_2.addWidget(self.m_tt, 1, 1, 1, 1)
        self.label_11 = QtGui.QLabel(self.centralwidget)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_2.addWidget(self.label_11, 4, 2, 1, 1)
        self.m_rt = ScientificDoubleSpinBox(self.centralwidget)
        self.m_rt.setObjectName(_fromUtf8("m_rt"))
        self.gridLayout_2.addWidget(self.m_rt, 3, 1, 1, 1)
        self.m_rp = ScientificDoubleSpinBox(self.centralwidget)
        self.m_rp.setObjectName(_fromUtf8("m_rp"))
        self.gridLayout_2.addWidget(self.m_rp, 4, 1, 1, 1)
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_2.addWidget(self.label_3, 5, 0, 1, 1)
        self.label_12 = QtGui.QLabel(self.centralwidget)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_2.addWidget(self.label_12, 5, 2, 1, 1)
        self.m_pp = ScientificDoubleSpinBox(self.centralwidget)
        self.m_pp.setObjectName(_fromUtf8("m_pp"))
        self.gridLayout_2.addWidget(self.m_pp, 2, 1, 1, 1)
        self.label_9 = QtGui.QLabel(self.centralwidget)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_2.addWidget(self.label_9, 1, 2, 1, 1)
        self.evla = QtGui.QDoubleSpinBox(self.centralwidget)
        self.evla.setMinimum(-90.0)
        self.evla.setMaximum(90.0)
        self.evla.setObjectName(_fromUtf8("evla"))
        self.gridLayout_2.addWidget(self.evla, 3, 3, 1, 1)
        self.stla = QtGui.QDoubleSpinBox(self.centralwidget)
        self.stla.setMinimum(-90.0)
        self.stla.setMaximum(90.0)
        self.stla.setObjectName(_fromUtf8("stla"))
        self.gridLayout_2.addWidget(self.stla, 5, 3, 1, 1)
        self.component = QtGui.QComboBox(self.centralwidget)
        self.component.setObjectName(_fromUtf8("component"))
        self.component.addItem(_fromUtf8(""))
        self.component.addItem(_fromUtf8(""))
        self.component.addItem(_fromUtf8(""))
        self.gridLayout_2.addWidget(self.component, 0, 3, 1, 1)
        self.label_16 = QtGui.QLabel(self.centralwidget)
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.gridLayout_2.addWidget(self.label_16, 0, 2, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 0, 0, 1, 1)
        self.seismogram = Qt4MplCanvas(self.centralwidget)
        self.seismogram.setMouseTracking(True)
        self.seismogram.setObjectName(_fromUtf8("seismogram"))
        self.gridLayout.addWidget(self.seismogram, 0, 1, 1, 1)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.window_start = QtGui.QSlider(self.centralwidget)
        self.window_start.setMaximum(10000)
        self.window_start.setProperty("value", 0)
        self.window_start.setOrientation(QtCore.Qt.Horizontal)
        self.window_start.setObjectName(_fromUtf8("window_start"))
        self.verticalLayout_3.addWidget(self.window_start)
        self.window_end = QtGui.QSlider(self.centralwidget)
        self.window_end.setMaximum(10000)
        self.window_end.setProperty("value", 10000)
        self.window_end.setSliderPosition(10000)
        self.window_end.setOrientation(QtCore.Qt.Horizontal)
        self.window_end.setObjectName(_fromUtf8("window_end"))
        self.verticalLayout_3.addWidget(self.window_end)
        self.calc_dispersion_button = QtGui.QPushButton(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.calc_dispersion_button.sizePolicy().hasHeightForWidth())
        self.calc_dispersion_button.setSizePolicy(sizePolicy)
        self.calc_dispersion_button.setObjectName(_fromUtf8("calc_dispersion_button"))
        self.verticalLayout_3.addWidget(self.calc_dispersion_button)
        self.gridLayout_4 = QtGui.QGridLayout()
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.label_19 = QtGui.QLabel(self.centralwidget)
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.gridLayout_4.addWidget(self.label_19, 1, 0, 1, 1)
        self.max_period = ScientificDoubleSpinBox(self.centralwidget)
        self.max_period.setMinimum(5.0)
        self.max_period.setMaximum(500.0)
        self.max_period.setProperty("value", 150.0)
        self.max_period.setObjectName(_fromUtf8("max_period"))
        self.gridLayout_4.addWidget(self.max_period, 1, 3, 1, 1)
        self.label_20 = QtGui.QLabel(self.centralwidget)
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.gridLayout_4.addWidget(self.label_20, 1, 2, 1, 1)
        self.label_18 = QtGui.QLabel(self.centralwidget)
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.gridLayout_4.addWidget(self.label_18, 1, 6, 1, 1)
        self.min_period = ScientificDoubleSpinBox(self.centralwidget)
        self.min_period.setMinimum(1.0)
        self.min_period.setMaximum(100.0)
        self.min_period.setProperty("value", 20.0)
        self.min_period.setObjectName(_fromUtf8("min_period"))
        self.gridLayout_4.addWidget(self.min_period, 1, 1, 1, 1)
        self.line_6 = QtGui.QFrame(self.centralwidget)
        self.line_6.setFrameShape(QtGui.QFrame.VLine)
        self.line_6.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.gridLayout_4.addWidget(self.line_6, 1, 4, 1, 1)
        self.alpha = ScientificDoubleSpinBox(self.centralwidget)
        self.alpha.setMinimum(1.0)
        self.alpha.setMaximum(500.0)
        self.alpha.setProperty("value", 100.0)
        self.alpha.setObjectName(_fromUtf8("alpha"))
        self.gridLayout_4.addWidget(self.alpha, 1, 7, 1, 1)
        self.mft_type = QtGui.QComboBox(self.centralwidget)
        self.mft_type.setObjectName(_fromUtf8("mft_type"))
        self.mft_type.addItem(_fromUtf8(""))
        self.mft_type.addItem(_fromUtf8(""))
        self.gridLayout_4.addWidget(self.mft_type, 1, 5, 1, 1)
        self.press_me = QtGui.QPushButton(self.centralwidget)
        self.press_me.setObjectName(_fromUtf8("press_me"))
        self.gridLayout_4.addWidget(self.press_me, 1, 8, 1, 1)
        self.gridLayout_4.setColumnStretch(0, 1)
        self.gridLayout_4.setColumnStretch(1, 2)
        self.gridLayout_4.setColumnStretch(2, 1)
        self.gridLayout_4.setColumnStretch(3, 2)
        self.gridLayout_4.setColumnStretch(4, 2)
        self.gridLayout_4.setColumnStretch(5, 2)
        self.gridLayout_4.setColumnStretch(6, 2)
        self.gridLayout_4.setColumnStretch(7, 2)
        self.verticalLayout_3.addLayout(self.gridLayout_4)
        self.gabormatrix = Qt4MplCanvas(self.centralwidget)
        self.gabormatrix.setObjectName(_fromUtf8("gabormatrix"))
        self.verticalLayout_3.addWidget(self.gabormatrix)
        self.verticalLayout_3.setStretch(0, 1)
        self.verticalLayout_3.setStretch(1, 1)
        self.verticalLayout_3.setStretch(4, 2)
        self.gridLayout.addLayout(self.verticalLayout_3, 1, 1, 1, 1)
        self.gridLayout.setColumnStretch(0, 1)
        self.gridLayout.setColumnStretch(1, 2)
        self.gridLayout.setRowStretch(0, 1)
        self.gridLayout.setRowStretch(1, 2)
        self.horizontalLayout.addLayout(self.gridLayout)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout.setStretch(2, 10)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1065, 27))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen_sac_file = QtGui.QAction(MainWindow)
        self.actionOpen_sac_file.setObjectName(_fromUtf8("actionOpen_sac_file"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.Save_dispersion_curve = QtGui.QAction(MainWindow)
        self.Save_dispersion_curve.setObjectName(_fromUtf8("Save_dispersion_curve"))
        self.menuFile.addAction(self.actionOpen_sac_file)
        self.menuFile.addAction(self.actionExit)
        self.menuFile.addAction(self.Save_dispersion_curve)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.open_sac_button.setText(_translate("MainWindow", "Open sac file", None))
        self.pushButton.setText(_translate("MainWindow", "Open CMT", None))
        self.open_instaseis_button.setText(_translate("MainWindow", "Open instaseis database", None))
        self.label_2.setText(_translate("MainWindow", "Noise", None))
        self.integrate_button.setText(_translate("MainWindow", "Integrate", None))
        self.differentiate_button.setText(_translate("MainWindow", "Differentiate", None))
        self.label_17.setText(_translate("MainWindow", "corners", None))
        self.corners.setItemText(0, _translate("MainWindow", "1", None))
        self.corners.setItemText(1, _translate("MainWindow", "2", None))
        self.corners.setItemText(2, _translate("MainWindow", "3", None))
        self.corners.setItemText(3, _translate("MainWindow", "4", None))
        self.label_15.setText(_translate("MainWindow", "max period", None))
        self.label_13.setText(_translate("MainWindow", "pre filter properties", None))
        self.label_14.setText(_translate("MainWindow", "min period", None))
        self.planet.setItemText(0, _translate("MainWindow", "Earth", None))
        self.planet.setItemText(1, _translate("MainWindow", "Mars", None))
        self.planet.setItemText(2, _translate("MainWindow", "Europa", None))
        self.zerophase.setText(_translate("MainWindow", "zerophase", None))
        self.label_4.setText(_translate("MainWindow", "m_rp", None))
        self.label_10.setText(_translate("MainWindow", "source longitude", None))
        self.label_5.setText(_translate("MainWindow", "m_rt", None))
        self.label_6.setText(_translate("MainWindow", "m_pp", None))
        self.label_8.setText(_translate("MainWindow", "m_rr", None))
        self.label_7.setText(_translate("MainWindow", "m_tt", None))
        self.label.setText(_translate("MainWindow", "source latitude", None))
        self.label_11.setText(_translate("MainWindow", "receiver longtiude", None))
        self.label_3.setText(_translate("MainWindow", "m_tp", None))
        self.label_12.setText(_translate("MainWindow", "receiver latitude", None))
        self.label_9.setText(_translate("MainWindow", "source depth", None))
        self.component.setItemText(0, _translate("MainWindow", "Z", None))
        self.component.setItemText(1, _translate("MainWindow", "R", None))
        self.component.setItemText(2, _translate("MainWindow", "T", None))
        self.label_16.setText(_translate("MainWindow", "component", None))
        self.calc_dispersion_button.setText(_translate("MainWindow", "Calculate Dispersion", None))
        self.label_19.setText(_translate("MainWindow", "Tmin", None))
        self.label_20.setText(_translate("MainWindow", "Tmax", None))
        self.label_18.setText(_translate("MainWindow", "alpha", None))
        self.mft_type.setItemText(0, _translate("MainWindow", "gaussian", None))
        self.mft_type.setItemText(1, _translate("MainWindow", "butterworth", None))
        self.press_me.setText(_translate("MainWindow", "press_me", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionOpen_sac_file.setText(_translate("MainWindow", "Open", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.Save_dispersion_curve.setText(_translate("MainWindow", "Save Dispersion Curve As", None))

from instaseis.gui.qt4mplcanvas import Qt4MplCanvas
from instaseis.gui.scientific_double_spin_box import ScientificDoubleSpinBox
