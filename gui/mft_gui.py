#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
gui for dispersion measurement with the 
multiple frequency technique (Dziewonski, 1969)
"""

import os
import sys
import imp
import obspy
import inspect
import instaseis
import numpy as np
import pyqtgraph as pg
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from PyQt4 import uic
from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import QVBoxLayout
from obspy import geodetics
from mpl_toolkits.basemap import Basemap
from glob import iglob
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from europa_seismo.europa_seismo.utils import gauss_filter

class Window(QtGui.QMainWindow):

    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = layout.Ui_MainWindow()
        self.ui.setupUi(self)
        self.instaseis = False
        self.instaseis_db = None
        self.sac_file = None
        self.stream = None
        self.stream_copy = None
        self.stream_slice = None
        self.source = None
        self.window_start = None
        self.window_end = None
        self.nbands = 100
        self.gabor_matrix = np.zeros((self.nbands,self.nbands))
        self.periods = np.linspace(self.ui.min_period.value(),
                                   self.ui.max_period.value(),
                                   self.nbands)
        self.period_pick = None
        self.vel_pick = None
        self.dist_km = None

        if self.ui.planet.currentText() == 'Earth':
            self.planet_radius = 6371.0
        elif self.ui.planet.currentText() == 'Mars':
            self.planet_radius = 3385.5
        elif self.ui.planet.currentText() == 'Europa':
            self.planet_radius = 1565.0
        
        m_rr = 3.98E13
        m_tt = 0.0
        m_pp = -3.98E13
        m_rt = 0.0
        m_rp = 0.0
        m_tp = 0.0
        self.ui.m_rr.setValue(m_rr)
        self.ui.m_tt.setValue(m_tt)
        self.ui.m_pp.setValue(m_pp)
        self.ui.m_rt.setValue(m_rt)
        self.ui.m_rp.setValue(m_rp)
        self.ui.m_tp.setValue(m_tp)
        self.ui.stlo.setValue(30.0)
        self.ui.stla.setValue(0.0)
        self.ui.evlo.setValue(0.0)
        self.ui.evla.setValue(0.0)
        self.ui.evdp.setValue(0.0)
        self.ui.window_start.setMinimum(0.0)
        self.ui.window_start.setMaximum(1.0)
        self.ui.window_start.setValue(0.0)
        self.ui.window_end.setMinimum(0.0)
        self.ui.window_end.setMaximum(1.0)
        self.ui.window_end.setValue(1.0)

        #TODO move this eventually so it only gets plotted after MFT analysis
        #self.plot_gabor()
        self.plot_map()
   
    def on_open_sac_button_released(self):
        cwd = os.getcwd()
        self.sac_file = str(QtGui.QFileDialog.getOpenFileName(
            self, "choose sac file", cwd))
        if not self.sac_file:
            return

        self.stream = obspy.read(self.sac_file)
        self.stream_copy = self.stream.copy()
        time = self.get_time_axis 
        self.ui.window_start.setMaximum(time[-1])
        self.ui.window_end.setMaximum(time[-1])
        self.ui.window_start.setValue(time[0])
        self.ui.window_end.setValue(time[-1])
        self.instaseis = False
        self.ui.stlo.setValue(self.stream[0].stats.sac['stlo'])
        self.ui.stla.setValue(self.stream[0].stats.sac['stla'])
        self.ui.evlo.setValue(self.stream[0].stats.sac['evlo'])
        self.ui.evla.setValue(self.stream[0].stats.sac['evla'])
        self.ui.evdp.setValue(self.stream[0].stats.sac['evdp'])
        self.time_window = None
        self.update()

    def on_open_instaseis_button_released(self):
        cwd = os.getcwd()
        self.folder = str(QtGui.QFileDialog.getExistingDirectory(
            self, "choose instaseis database folder", cwd))
        if not self.folder:
            return

        self.instaseis_db = instaseis.open_db(self.folder)
        self.source = instaseis.Source(latitude=self.ui.evla.value(),
                                       longitude=self.ui.evlo.value(),
                                       depth_in_m=self.ui.evdp.value()*1000.,
                                       m_rr = self.ui.m_rr.value(),
                                       m_tt = self.ui.m_tt.value(),
                                       m_pp = self.ui.m_pp.value(),
                                       m_rt = self.ui.m_rt.value(),
                                       m_rp = self.ui.m_rp.value(),
                                       m_tp = self.ui.m_pp.value())
        self.receiver = instaseis.Receiver(latitude=self.ui.stla.value(),
                                           longitude=self.ui.stlo.value())
        self.stream = self.instaseis_db.get_seismograms(source=self.source,
                                       receiver=self.receiver,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
        self.stream[0].stats.sac = {}
        self.stream[0].stats.sac['o'] = 0.0
        self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
           float(self.ui.evla.value()),float(self.ui.evlo.value()),
           float(self.ui.stla.value()),float(self.ui.stlo.value())) 
        self.stream_copy = self.stream.copy()

        time = self.get_time_axis
        self.ui.window_start.setMaximum(time[-1])
        self.ui.window_end.setMaximum(time[-1])
        self.ui.window_start.setValue(time[0])
        self.ui.window_end.setValue(time[-1])
        self.instaseis = True
        self.plot_map()
        self.update()

    def on_integrate_button_released(self):
        if not self.stream:
            return
        self.stream = self.stream.integrate()
    
        self.update()

    def on_differentiate_button_released(self):
        if not self.stream:
            return
        self.stream = self.stream.differentiate()
        self.update()

    def on_calc_dispersion_button_released(self):
        #self.gabor_matrix = np.zeros((self.stream_slice.stats.npts,
                                      #self.stream_slice.stats.npts))
        #self.gabor_matrix[:,0] = self.stream_slice[0].data
        #for i in range(0,self.gabor_matrix.shape[0]):
        #    self.gabor_matrix[i,:] = self.stream_slice.data 
        self.multiple_filter()
        self.plot_gabor()
        self.update()

    def multiple_filter(self,Tmin=20.0,Tmax=200.0):
        #time = self.get_time_axis
        self.km_per_deg = (2.*np.pi*self.planet_radius / 360.0)
        self.dist_km = self.stream[0].stats.sac['gcarc'] * self.km_per_deg

        try:
            self.time_window = np.linspace(self.ui.window_start.value(),
                                           self.ui.window_end.value(),
                                           len(self.stream_slice.data))
        except TypeError:
            time = np.linspace(100,200,100)

        #self.periods = np.linspace(Tmin,Tmax,nbands)
        self.gabor_matrix = np.zeros((len(self.time_window),len(self.periods))) 
        print 'GABOR MATRIX SHAPE', self.gabor_matrix.shape
        self.period_pick = []
        self.vel_pick = []

        if self.ui.mft_type.currentText() == 'gaussian':
            for i,period in enumerate(self.periods):
                self.stream_slice_copy = self.stream_slice.copy()
                tr = self.stream_slice_copy
                dcol = gauss_filter(tr.data,tr.stats.sampling_rate,
                    w_0=(1/period),alpha=self.ui.alpha.value())
                env = np.abs(hilbert(dcol.real))
                self.gabor_matrix[:,i] = env

                if np.max(env) > 0.0:
                    self.vel_pick.append(self.dist_km/self.time_window[np.argmax(env)])
                    self.period_pick.append(period)

        elif self.ui.mft_type.currentText() == 'butterworth':
            for i,period in enumerate(self.periods):
                Tstart = period / 1.2
                Tend = period + (period / 2.0)
                freqmin = 1. / Tend
                freqmax = 1. / Tstart
                self.stream_slice_copy = self.stream_slice.copy()

                if freqmin > 0:
                    tr = self.stream_slice_copy.filter('bandpass',
                        freqmin=freqmin,freqmax=freqmax,corners=4,
                        zerophase=True)
                else:
                    tr = self.stream_slice_copy.filter('lowpass',
                        freq=freqmax)

                dcol = tr.data
                env = np.abs(hilbert(dcol.real))
                self.gabor_matrix[:,i] = env

                if np.max(env) > 0.0:
                    self.vel_pick.append(self.dist_km/self.time_window[np.argmax(env)])
                    self.period_pick.append(period)
        else:
            raise ValueError('filter type',kind,' not implemented')

        self.update()

    @property
    def mt_source(self):
        m_rr = float(self.ui.m_rr.value())
        m_tt = float(self.ui.m_tt.value())
        m_pp = float(self.ui.m_pp.value())
        m_rt = float(self.ui.m_rt.value())
        m_rp = float(self.ui.m_rp.value())
        m_tp = float(self.ui.m_tp.value())
        source = instaseis.Source(latitude=self.ui.evla.value(),
                                  longitude=self.ui.evlo.value(),
                                  depth_in_m=self.ui.evdp.value(),
                                  m_rr=m_rr,
                                  m_tt=m_tt,
                                  m_pp=m_pp,
                                  m_rt=m_rt,
                                  m_rp=m_rp,
                                  m_tp=m_tp)
        return source

    @property
    def instaseis_receiver(self):
        longitude=float(self.ui.stlo.value())
        latitude=float(self.ui.stla.value())
        rec = instaseis.Receiver(latitude=latitude,
                                 longitude=longitude)
        return rec

    @property
    def get_time_axis(self):
        time = np.linspace(self.stream[0].stats.sac['o'],
                           self.stream[0].stats.npts*self.stream[0].stats.delta,
                           self.stream[0].stats.npts)
        return time

    def on_m_rr_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream_copy = self.stream.copy()
        self.update()

    def on_m_tt_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_m_pp_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_m_rt_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_m_rp_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_m_tp_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_evdp_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_evlo_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_evla_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_stlo_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_stla_valueChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_min_period_valueChanged(self, *args):
        self.periods=np.linspace(self.ui.min_period.value(),
                                 self.ui.max_period.value(),
                                 self.nbands)
        self.update()

    def on_max_period_valueChanged(self, *args):
        self.periods=np.linspace(self.ui.min_period.value(),
                                 self.ui.max_period.value(),
                                 self.nbands)
        self.update()

    def on_component_currentIndexChanged(self, *args):
        if self.instaseis:
            src = self.mt_source
            rec = self.instaseis_receiver
            self.stream = self.instaseis_db.get_seismograms(source=src,
                                       receiver=rec,
                                       components=str(self.ui.component.currentText()),
                                       kind='displacement',
                                       remove_source_shift=True)
            self.stream[0].stats.sac = {}
            self.stream[0].stats.sac['o'] = 0.0
            self.stream[0].stats.sac['gcarc'] = geodetics.locations2degrees(
            float(self.ui.evla.value()),float(self.ui.evlo.value()),
            float(self.ui.stla.value()),float(self.ui.stlo.value())) 
            self.stream_copy = self.stream.copy()
        self.update()

    def on_mft_type_currentIndexChanged(self, *args):
        self.update()

    def on_planet_currentIndexChanged(self, *args):

        if self.ui.planet.currentText() == 'Earth':
            self.planet_radius = 6371.0
        elif self.ui.planet.currentText() == 'Mars':
            self.planet_radius = 3385.5
        elif self.ui.planet.currentText() == 'Europa':
            self.planet_radius = 1565.0
        print self.planet_radius
        self.update()

    def on_window_start_valueChanged(self, *args):
        t_start = float(self.ui.window_start.value())
        t_end = float(self.ui.window_end.value())
        self.stream_slice = self.stream[0].slice(
            self.stream[0].stats.starttime + t_start,
            self.stream[0].stats.starttime + t_end)
        self.stream = self.stream_copy
        self.update()

    def on_window_end_valueChanged(self, *args):
        t_start = float(self.ui.window_start.value())
        t_end = float(self.ui.window_end.value())
        self.stream_slice = self.stream[0].slice(
            self.stream[0].stats.starttime + t_start,
            self.stream[0].stats.starttime + t_end)
        self.stream = self.stream_copy
        self.update()

    def on_alpha_valueChanged(self, *args):
        self.update()

    def update(self, force=False):
        #plot_widget = self.ui.seismogram
        #plot_widget.clear()
        #time = np.linspace(self.stream[0].stats.sac['o'],
        #                   self.stream[0].stats.npts*self.stream[0].stats.delta,
        #                   self.stream[0].stats.npts)
        #plot_widget.plot(time,self.stream[0].data,pen='b')
        self.plot_seismogram()
        #self.plot_gabor()

    def plot_seismogram(self):

        time=self.get_time_axis

        self.mpl_seis_figure = self.ui.seismogram.fig
        if hasattr(self,'mpl_seis_ax'):
            self.mpl_seis_ax.clear()
        self.mpl_seis_ax = self.mpl_seis_figure.add_axes([0.1,0.15,0.90,0.70]) 
        self.mpl_seis_ax.plot(time,self.stream[0].data,c='b')
        self.mpl_seis_ax.set_xlabel('time (s)')

        self.mpl_seis_ax.axvline(float(self.ui.window_start.value()),c='k')
        self.mpl_seis_ax.axvline(float(self.ui.window_end.value()),c='k')
        self.mpl_seis_ax.set_xlim([self.ui.window_start.value(),
                                    self.ui.window_end.value()])
        self.mpl_seis_figure.canvas.draw()

    def plot_map(self):
        self.mpl_map_figure = self.ui.map.fig
        self.mpl_map_ax = self.mpl_map_figure.add_axes([0.01,0.01,0.99,0.99])
        self.map = Basemap(projection='moll',lon_0=0,resolution='c',
                           ax = self.mpl_map_ax)
        self.map.drawcoastlines()
        self.mpl_map_figure.canvas.draw()

    def plot_gabor(self):

        if hasattr(self,'mpl_gabor_ax'):
            self.mpl_gabor_ax.clear()

        self.mpl_gabor_figure = self.ui.gabormatrix.fig
        self.mpl_gabor_ax = self.mpl_gabor_figure.add_axes([0.1,0.1,0.8,0.6])

        if self.stream is not None:
            self.km_per_deg = (2.*np.pi*self.planet_radius / 360.0)
            self.dist_km = self.stream[0].stats.sac['gcarc'] * self.km_per_deg
            vel_min = self.dist_km / float(self.ui.window_start.value())
            vel_max = self.dist_km / float(self.ui.window_end.value())

        #extent = [20.,200.,vel_max,vel_min]
        #self.mpl_gabor_ax.imshow(self.gabor_matrix.T,
        #    aspect='auto',extent=extent,cmap='jet')
        #self.mpl_gabor_ax.set_xlabel('period (s)')
        #self.mpl_gabor_ax.set_ylabel('velocity (km/s)')

        self.gabor_matrix /= np.max(self.gabor_matrix)
        self.mpl_gabor_ax.contourf(self.periods,self.dist_km/self.time_window,self.gabor_matrix,
                                   cmap='jet',levels=np.linspace(0,1,50))

        if self.vel_pick is not None:
            self.mpl_gabor_ax.scatter(self.period_pick,self.vel_pick,
                                      c='k',marker='+')

        veloc = self.dist_km / self.time_window
        #ftest = np.loadtxt('/home/romaguir/Tools/europa_seismo/data/prem_grpvel_minos.txt')
        #self.mpl_gabor_ax.plot(ftest[:,0],ftest[:,1],c='k',alpha=0.75)
        self.mpl_gabor_ax.set_xlim([self.ui.min_period.value(),self.ui.max_period.value()])
        self.mpl_gabor_ax.set_ylim([np.min(veloc),min(5,np.max(veloc))])
        self.mpl_gabor_ax.set_xlabel('period (s)')
        self.mpl_gabor_ax.set_ylabel('velocity (km/s)')
        self.mpl_gabor_figure.canvas.draw()

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
