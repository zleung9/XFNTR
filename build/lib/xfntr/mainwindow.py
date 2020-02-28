import sys
import os
# Use absolute path instead of relative path ('./') to avoid trouble when installed by pip
dir_path = os.path.dirname(os.path.realpath(__file__))
UI_path = dir_path + '/GUI/'
import time
import multiprocessing


######################################################
# This block of code is needed for properly working with PyInstaller

# import the following three modules in order to work with PyInstaller

# mplwidget is imported explicitly here because PyInstaller needs to find it.

# Define function to import external files when using PyInstaller.
# https://stackoverflow.com/questions/37888581/pyinstaller-ui-files-filenotfounderror-errno-2-no-such-file-or-directory
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        # when not bundled by PyInstaller, normal method is used.
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


#######################################################
from collections import OrderedDict

from PyQt5 import uic
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from scipy.interpolate import interp1d
from scipy.special import *
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
default_cycler = (cycler(color=['b','g','c','m','y','k']))
plt.rc('axes', prop_cycle=default_cycler)
plt.rc('lines',markersize=3, linewidth=1)

import lmfit as lm # fitting module
import periodictable as pdtb # a module for periodic table
r_e = pdtb.constants.electron_radius * 1e10  # classical electron radius, in A
N_A = pdtb.constants.avogadro_number  # Avogadro number, unitless
k_B = 1.38065e-23  # Boltzman constant, in J/K

# user defined module
import flu_routines_new as fl

# Here the absolute path is used because PyInstaller needs to find it.
(Ui_MainWindow, QMainWindow) = uic.loadUiType(resource_path(UI_path + 'mainwindow.ui'))


class myThread(QThread):
    def __init__(self,func):
        QThread.__init__(self)
        self.func = func

    def run(self):
        self.func
        
def fluErrorFitSingleCore2(i, value_list, sh, qz, pname, flu_par, flucal_par, data_to_fit):
        flu_par[pname].value = value_list[i]  # change value of the chosen parameter
        fluerr_result = lm.minimize(fl.flu2min, flu_par,
                                    args=((sh, qz), flucal_par),
                                    kws={'data': data_to_fit[:, 1], 'eps': data_to_fit[:, 2]})
        return ([i, value_list[i], fluerr_result.nfree, fluerr_result.redchi])

def multiCore(func, iterable):
    pool = multiprocessing.Pool()
    result = pool.map(func, range(len(iterable)))
    pool.close()
    pool.join()
    return result


class MainWindow (QMainWindow):
    """MainWindow inherits QMainWindow"""

    def __init__(self, parent = None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.tabWidget.setCurrentIndex(1)
        self.directory=os.getcwd()

        self.halftab = '    '
        self.flusavefitindex = 0
        self.flufiles = []
        self.flufitfiles = []
        self.fludata = []
        self.flufitdata = []
        self.beam = 'uniform'
        self.flucal_par = OrderedDict()
        self.flu = 0
        self.selectedflufiles_rows = []
        self.selectedflufitfiles_rows = []
        self.eleradius = pdtb.constants.electron_radius*1e10
        self.avoganum = pdtb.constants.avogadro_number
        self.boltzmann = 1.38065e-23
        self.errorlist = np.array([[1, 1.074], [2, 1.204], [3, 1.222],
                                [4, 1.220], [5, 1.213], [6, 1.205], 
                                [7, 1.198], [8, 1.191], [9, 1.184], 
                                [10, 1.178], [11, 1.173], [12, 1.168], 
                                [13, 1.163], [14, 1.159], [15, 1.155], 
                                [16, 1.151], [17, 1.148], [18, 1.145], 
                                [19, 1.142], [20, 1.139], [22, 1.134], 
                                [24, 1.129], [26, 1.125], [28, 1.121], 
                                [30, 1.118], [32, 1.115], [34, 1.112], 
                                [36, 1.109], [38, 1.106], [40, 1.104], 
                                [42, 1.102], [44, 1.100], [46, 1.098], 
                                [48, 1.096], [50, 1.094], [60, 1.087], 
                                [70, 1.081], [80, 1.076], [90, 1.072], 
                                [100, 1.069], [120, 1.063], [140, 1.059], 
                                [160, 1.055], [180, 1.052]]) #, [3000, 1.050]])

        self.setupUI()
        self.updatePar()
        self.debugErr()

    def setupUI(self):
        self.ui.addflufilePB.clicked.connect(self.addFluFile)
        self.ui.flufileLW.itemSelectionChanged.connect(self.updateSelectedFluFile)
        self.ui.rmflufilePB.clicked.connect(self.removeFluFile)
        self.ui.addflufitfilePB.clicked.connect(self.addFluFitFile)
        self.ui.flufitfileLW.itemSelectionChanged.connect(self.updateSelectedFluFitFile)
        self.ui.rmflufitfilePB.clicked.connect(self.removeFluFitFile)
        self.ui.fluxaxisCB.currentIndexChanged.connect(self.updateUI)
        self.ui.fluqcCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flulineCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flulegendCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flulegendlocCoB.currentIndexChanged.connect(self.updateFluPlot)
        self.ui.flulogyCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flugridCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flushowCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flucompCB.stateChanged.connect(self.updateFluPlot)
        self.ui.flusimuPB.clicked.connect(self.updateFluCal)
        self.ui.flufitPB.clicked.connect(self.fitFlu)
        self.ui.flusaveCB.activated.connect(self.saveFlu)
        self.ui.fluloadCB.activated.connect(self.loadFlu)
        self.ui.insflusubPB.clicked.connect(self.insFluIon)
        self.ui.rmflusubPB.clicked.connect(self.rmFluIon)
        self.ui.fluErrPB.clicked.connect(self.fluErrorInit)
        self.ui.flulimitPB.clicked.connect(self.setupLimitsUI)
        # connect system parameter signals
        self.ui_syspar = OrderedDict(
            [('E_inc',self.ui.fluIncEnLE),
             ('E_emt',self.ui.fluEmtEnLE),
             ('mu_top_inc',self.ui.flumutopincLE),
             ('mu_top_emt', self.ui.flumutopemtLE),
             ('mu_bot_inc', self.ui.flumubotincLE),
             ('mu_bot_emt',self.ui.flumubotemtLE),
             ('rho_top', self.ui.flurhotopLE),
             ('rho_bot', self.ui.flurhobotLE),
             ('width', self.ui.fluwidthLE),
             ('det_len',self.ui.fludetLE)]
        )
        for p,u in self.ui_syspar.items():
            u.returnPressed.connect(self.updatePar)
            u.returnPressed.connect(self.updateFluCal)
        self.ui.flubmpfCombo.currentIndexChanged.connect(self.updatePar)

        # connect fitting parameter signals
        self.ui_params = OrderedDict(
            [('losc', [self.ui.flubotscaleLE, self.ui.flubotscaleCB]),
             ('hisc', [self.ui.flutopscaleLE, self.ui.flutopscaleCB]),
             ('lobk', [self.ui.flubotbulkLE, self.ui.flubotbulkCB]),
             ('upbk', [self.ui.flutopbulkLE, self.ui.flutopbulkCB]),
             ('surd', [self.ui.flusurdLE, self.ui.flusurdCB]),
             ('bg',   [self.ui.flubgLE, self.ui.flubgCB]),
             ('qoff', [self.ui.fluqoffLE, self.ui.fluqoffCB]),
             ('curv', [self.ui.flucurvLE, self.ui.flucurvCB]),
             ('loff', [self.ui.fluloffLE, self.ui.fluloffCB]),
             ('soff', [self.ui.flusoffLE, self.ui.flusoffCB])]
        )
        for p,u in self.ui_params.items():
            u[0].returnPressed.connect(self.updatePar)
            u[0].returnPressed.connect(self.updateFluCal)
            u[1].stateChanged.connect(self.updatePar)

        # set up a container for parameter limits
        self.par_limits = OrderedDict()
        for p,u in self.ui_params.items():
            if p in ['lobk','upbk','surd','bg','curv']:
                self.par_limits[p] = [0, 1000]
            elif p in ['hisc','losc']:
                self.par_limits[p] = [1e-10, 1000]
            else:
                self.par_limits[p] = [None, None]

    def updateUI(self,index=None, fresh=True):
        '''
        This function also works as a slot for qcombobox signal, so the first arg/kws by default
        is reserved for receiving signal argument, which can be either index or text.
        '''
        print("{}, {}".format(index,fresh))
        # set text for fitting parameters
        for p, u in self.ui_params.items():
            if p in ['curv', 'loff']:
                u[0].setText(format(self.flu_par[p].value, '.4f'))
            else:
                u[0].setText(format(self.flu_par[p].value, '.2e'))

        # set text for system parameters
        for p, u in self.ui_syspar.items():
            u.setText(format(self.sys_par[p], '.4f'))

        # set beam profile
        index = self.ui.flubmpfCombo.findText(self.beam)
        self.ui.flubmpfCombo.setCurrentIndex(index)
        self.xaxis = self.ui.fluxaxisCB.currentText()

        if fresh == True:
            if self.xaxis == 'Qz':
                _xrange = '0.005:0.016'
                self.ui.fluloffCB.setCheckable(True)
                self.ui.fluloffCB.setText('L2 offset')
                self.ui.fluloffqzlabel.setText('mm')
            elif self.xaxis == 'Sh':
                _xrange = '-0.1:0.1'
                self.ui.fluloffCB.setText('Qz')
                self.ui.fluloffLE.setText(format(0.006, '.4f'))
                self.ui.fluloffqzlabel.setText(self.ui.fluqofflabel.text())
            self.ui.fluSimuRangeLE.setText(_xrange)
            self.ui.fluFitRangeLE.setText(_xrange)
        else:
            try:
                _simu_range = ':'.join(str(self.flu_simu_range)[1:-1].split(','))
                _fit_range = ':'.join(str(self.flu_fit_range)[1:-1].split(','))
                self.ui.fluSimuRangeLE.setText(_simu_range)
                self.ui.fluFitRangeLE.setText(_fit_range)
                self.ui.fluloffLE.setText(str(self.qz[0]))
            except Exception as e:
                self.updateUI(fresh=True)

        self.updatePar()

    def updatePar(self):  #initialize the flu parameters

        # system parameters for fluorescence
        self.beam = str(self.ui.flubmpfCombo.currentText()) # beam profile
        self.sys_par = OrderedDict()
        for p, u in self.ui_syspar.items():
            self.sys_par[p] = float(u.text())
        self.sys_par['beam'] = self.ui.flubmpfCombo.currentText()
        self.sys_par['span'] = 75.6 # the length of sample cell, in mm.

        # fitting parameters for fluorescence
        self.flu_par = lm.Parameters()
        try:
            for name,par in self.ui_params.items():
                self.flu_par.add(name, value=float(par[0].text()), vary=par[1].isChecked(),
                                       min=self.par_limits[name][0], max=self.par_limits[name][1],
                                       expr=None, brute_step=None)
        except ValueError as VE:
            print("ValueError: ", VE)
        else:
            # fitting (if any) parameters for reflectivity
            self.ref_par = lm.Parameters()
                # add tuples:       (NAME       VALUE   VARY MIN  MAX  EXPR BRUTE_STEP)
            self.ref_par.add_many( ('rho_t', self.sys_par['rho_top'], 0, None, None, None, None),
                                   ('rho_b', self.sys_par['rho_bot'], 0, None, None, None, None),
                                   ('mu_t', self.sys_par['mu_top_inc'], 0, None, None, None, None),
                                   ('mu_b', self.sys_par['mu_bot_inc'], 0, None, None, None, None),
                                   ('sigma0', 3.0, 0, None, None, None, None),
                                   ('q_off', 0, 0, None, None, None, None ))
            # info of element in the system
            self.flu_elements = [['Eu', 1, 0.947]]  # name, composition, Ionic Radius(A)

        # update the list of sh and qz
        try:
            self.flu_simu_range = [float(i) for i in str(self.ui.fluSimuRangeLE.text()).split(':')]
            self.flu_fit_range = [float(i) for i in str(self.ui.fluFitRangeLE.text()).split(':')]
        except:
            print("Error: Check if the range is right.")
        self.xaxis = self.ui.fluxaxisCB.currentText()
        if self.xaxis == 'Qz':
            self.sh = np.array([0])
            self.qz = np.linspace(self.flu_simu_range[0], self.flu_simu_range[1], 200)
        elif self.xaxis == 'Sh':
            self.qz = np.array([float(self.ui.fluloffLE.text())])
            self.sh = np.linspace(self.flu_simu_range[0], self.flu_simu_range[1], 200)

        # update parameters for flurescence calculation
        self.flucal_par = fl.update_flu_parameters(self.flucal_par,
                                                   self.flu_par,
                                                   self.sys_par,
                                                   self.flu_elements)

    def updateLimits(self):
        ui = self.ui_limits
        try:
            for p, u in self.ui_par_limits.items():
                if u[0].isChecked():
                    try:
                        self.par_limits[p][0] = float(u[1].text())
                        self.par_limits[p][1] = float(u[2].text())
                        assert self.par_limits[p][1] >= self.par_limits[p][0]
                    except ValueError as e:
                        print("{} Please provide a valid limit for '{}'. ".format(e,u[0].text()))
                        raise
                    except AssertionError:
                        print("Max should be larger than Min for {}".format(u[0].text()))
                        raise
                else:
                    self.par_limits[p][0] = None
                    self.par_limits[p][1] = None
        except:
            print("An Error occurs... See above")
        else:
            self.updatePar()
            ui.close()

    def setupLimitsUI(self):

        ui = uic.loadUi(UI_path + 'err4.ui', QDialog(self))
        ui.cancelPB.clicked.connect(ui.close)
        ui.confirmPB.clicked.connect(self.updateLimits)
        self.ui_par_limits = OrderedDict(
            [('losc', [ui.flubotscaleCB, ui.minbotscaleLE, ui.maxbotscaleLE]),
             ('hisc', [ui.flutopscaleCB, ui.mintopscaleLE, ui.maxtopscaleLE]),
             ('lobk', [ui.flubotbulkCB, ui.minbotbulkLE, ui.maxbotbulkLE]),
             ('upbk', [ui.flutopbulkCB, ui.mintopbulkLE, ui.maxtopbulkLE]),
             ('surd', [ui.flusurdCB, ui.minsurdLE, ui.maxsurdLE]),
             ('bg', [ui.flubgCB, ui.minbgLE, ui.maxbgLE]),
             ('qoff', [ui.fluqoffCB, ui.minqoffLE, ui.maxqoffLE]),
             ('curv', [ui.flucurvCB, ui.mincurvLE, ui.maxcurvLE]),
             ('loff', [ui.fluloffCB, ui.minloffLE, ui.maxloffLE]),
             ('soff', [ui.flusoffCB, ui.minsoffLE, ui.maxsoffLE])]
        )
        for p, u in self.ui_par_limits.items():
            if self.par_limits[p][0] is None and self.par_limits[p][1] is None:
                u[0].setChecked(False)
                continue
            else:
                u[0].setChecked(True)
            try:
                u[1].setText(format(self.par_limits[p][0], '.2e'))
            except TypeError:
                pass  # leave empty if limit is None
            try:
                u[2].setText(format(self.par_limits[p][1], '.2e'))
            except TypeError:
                pass
        self.ui_limits = ui
        ui.show()


    def addFluFile(self): #add flu files into the listwidget and deselect all flu files in the listwidget

        f, _ = QFileDialog.getOpenFileNames(
            caption='Select Multiple Fluorescence Files to import',
            directory=self.directory,
            filter='Flu Files (*.flu*;*_flu.txt)'
        )
        self.flufiles = self.flufiles + f
        self.directory = str(QFileInfo(self.flufiles[0]).absolutePath())
        self.updateFluFile()

    def updateFluFile(self): # update flu files in the listwidget
        self.ui.flufileLW.clear()
        for i, f in enumerate(self.flufiles):
            try:
                self.ui.flufileLW.addItem('#'+str(i+1)+self.halftab+str(f.split('\\')[-2])+'\\'+str(f.split('\\')[-1]))
            except:
                self.ui.flufileLW.addItem('#'+str(i+1)+self.halftab+str(f.split('/')[-2])+'/'+str(f.split('/')[-1]))

    def updateSelectedFluFile(self): #update the selected flu files in the listwidget
        self.fludata = []
        selectedflufiles = self.ui.flufileLW.selectedItems()
        self.selectedflufiles_rows = [self.ui.flufileLW.row(item) for item in selectedflufiles]
        self.selectedflufiles_rows.sort()
        if len(selectedflufiles) != 0:
            try:
                for i, r in enumerate(self.selectedflufiles_rows):
                    data = np.loadtxt(str(self.flufiles[r]),comments='#')
                    for d in data: # replace zero error bar with 10% error
                        if d[2]==0:
                            print('Error bar replaced with 10% of value for entry {}'.format(d))
                            d[2] = float(d[1]) / 10
                    print('\n')

                    self.fludata.append(data)
            except OSError as e:
                print(e)
        self.updateFluPlot()

    def removeFluFile(self): #remove flu files in the listwidget and deselect all flu files in the listwidget

        to_del = [self.ui.flufileLW.row(item) for item in self.ui.flufileLW.selectedItems()]
        self.flufiles = [f for i,f in enumerate(self.flufiles) if i not in to_del]

        self.ui.flufileLW.clear()
        self.updateFluFile()

    def addFluFitFile(self): #add flu fit files into the listwidget and deselect flu fit files in the listwidget
        try:
            f, _ = QFileDialog.getOpenFileNames(
                caption = 'Select Multiple Fluorescence Fit Files to import',
                directory = self.directory,
                filter = 'FIT Files (*.fit*; *_fit.txt)'
            )
            self.flufitfiles = self.flufitfiles + f
            self.directory = str(QFileInfo(self.flufitfiles[0]).absolutePath())
            self.updateFluFitFile()
        except IndexError as IE:
            pass # ignore IndexError
        except:
            print("Something went wrong when reading fit files!")
    def updateFluFitFile(self): #update flu fit files in the listwidget
        self.ui.flufitfileLW.clear()
        for i, f in enumerate(self.flufitfiles):
            try:
                self.ui.flufitfileLW.addItem(
                    '#' + str(i + 1) + self.halftab + str(f.split('\\')[-2]) + '\\' + str(f.split('\\')[-1]))
            except:
                self.ui.flufitfileLW.addItem(
                    '#' + str(i + 1) + self.halftab + str(f.split('/')[-2]) + '/' + str(f.split('/')[-1]))

    def updateSelectedFluFitFile(self): #update the selected flu fit files in the listwidget
        self.flufitdata = []
        selectedflufitfiles = self.ui.flufitfileLW.selectedItems()
        self.selectedflufitfiles_rows = [self.ui.flufitfileLW.row(item) for item in selectedflufitfiles]
        self.selectedflufitfiles_rows.sort()
        try:
            if len(selectedflufitfiles) != 0:
                for i, r in enumerate(self.selectedflufitfiles_rows):
                    data = np.loadtxt(str(self.flufitfiles[r]),comments='#')
                    self.flufitdata.append(data)
        except OSError as OE:
            print(OE)
        self.updateFluPlot()

    def removeFluFitFile(self):  #remove flu fit files in the listwidget and deselect all flu fit files in the listwidget

        to_del = [self.ui.flufitfileLW.row(item) for item in self.ui.flufitfileLW.selectedItems()]
        self.flufitfiles = [f for i, f in enumerate(self.flufitfiles) if i not in to_del]

        self.ui.flufitfileLW.clear()
        self.updateFluFitFile()

    def updateFluPlot(self): #update the plot in the flu plotwidget

        ax1 = self.ui.fluPW.canvas.ax
        ax1.clear()

        if self.ui.flulineCB.isChecked():
            ls = '-'
        else:
            ls = ''

        if len(self.fludata) != 0: #plot flu files
            for i,d in enumerate(self.fludata):
                ax1.errorbar(d[:,0],d[:,1],yerr=d[:,2],
                             marker='o', ls=ls,
                             label='#'+str(i+1))

        if len(self.flufitdata) != 0: #plot flu fit files
            for i, d in enumerate(self.flufitdata):
                ax1.plot(d[:, 0], d[:, 1],marker='', ls='-', label=' fit #'+str(i + 1))


        self.xaxis = self.ui.fluxaxisCB.currentText()
        if self.ui.flushowCB.isChecked():
            if np.all(self.flu==0):
                print('Please print simulate button first!!')
                return
            else:
                if self.xaxis == 'Qz':
                    x = self.qz
                    y = self.flu[0, :, 2:]
                    x_range = [x[0]-0.001, x[-1]+0.001]
                    x_label = r'$Q_z$' + ' ' + r'$[\AA^{-1}]$'
                    if self.ui.fluqcCB.isChecked():
                        ax1.axvline(self.qc, color='black', alpha=0.5)
                elif self.xaxis == 'Sh':
                    x = self.sh
                    y = self.flu[:, 0, 2:]
                    x_range = [x[0]-0.001, x[-1]+0.001]
                    x_label = r'$\Delta sh$' + ' ' + r'$[mm]$'
            try:
                ax1.set_xlabel(x_label)
                ax1.set_ylabel(r'$Intensity [a.u.]$')
                ax1.set_xlim(x_range)
                ax1.plot(x, y[:,0], ls='-', label='total', color='r')
                if self.ui.flucompCB.isChecked():
                    ax1.plot(x, y[:,1], ls='-', label='water',color='b', alpha=0.5)
                    ax1.plot(x, y[:,2], ls='-', label='interface',color='purple', alpha=0.5)
                    ax1.plot(x, y[:,3], ls='-', label='oil', color='g', alpha=0.5)

            except ValueError as VE:
                print(VE)
        if self.ui.flulegendCB.isChecked():
            ax1.legend(loc = str(self.ui.flulegendlocCoB.currentText()),
                       frameon=False,
                       scatterpoints=0,
                       numpoints=1)
        if self.ui.flugridCB.isChecked():
            ax1.grid(1)
        if self.ui.flulogyCB.isChecked():
            ax1.set_yscale('log')
        else:
            ax1.set_yscale('linear')

        self.ui.fluPW.canvas.draw()

    def setFluPlotScale(self): #set the scale of each data in the flu plot
        if len(self.selectedflufiles_rows)+len(self.selectedflufitfiles_rows)==0:
            self.messageBox('Warning:: No Fluorescence or Fit files selected!')
        else:
            row_flu=len(self.selectedflufiles_rows)
            row_fit=len(self.selectedflufitfiles_rows)
            row=row_flu+row_fit
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi(UI_path + 'plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Fluorescence Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header
            for i in range(row_flu):
                vlabel.append("Flu #"+str(self.selectedflufiles_rows[i]+1))
            for i in range(row_fit):
                vlabel.append("Fit #"+str(self.selectedflufitfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_flu):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.fluscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            for i in range(row_fit):
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i+row_flu,j,QTableWidgetItem(str(self.flufitscale[i][j])))
                    self.uiplotscale.scaleTW.item(i+row_flu,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW, SIGNAL('cellChanged(int,int)'), self.updateFluPlotScale) #update the flu scale and plot
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale) #close the scale setup window

    def updateFluPlotScale(self): #update the scale of each data in the flu plot
        row_flu=len(self.selectedflufiles_rows)
        row_fit=len(self.selectedflufitfiles_rows)
        self.fluscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_flu)]
        self.flufitscale=[[float(str(self.uiplotscale.scaleTW.item(i+row_flu,j).text())) for j in range(4)] for i in range(row_fit)]
        self.updateFluPlot()

    def updateFluCal(self): # calculate the flu  based on current parameters.

        self.updatePar()
        # p is the parameter set taken by the core function, must be initialized before updatding
        p = OrderedDict()
        p = fl.update_flu_parameters(p, self.flu_par, self.sys_par,self.flu_elements)
        if not self.ui.flushowCB.isChecked():
            return # if show is not checked, do nothing.
        if self.xaxis == 'Qz':
            self.qc = np.sqrt(2*(p['ibDt']-p['itDt'])) * 2 * p['k0']

        self.flu = fl.flu2min(self.flu_par, (self.sh, self.qz), p)

        self.updateFluPlot()

    def fitFlu(self, uncertainty_calculation=False):

        selectedflufiles = self.ui.flufileLW.selectedItems()
        try:
            assert len(selectedflufiles) == 1
        except:
            print("Error: please select one data to fit.")
            return
        data = self.fludata[0]

        self.updatePar()

        self.data_to_fit = data[(data[:,0] >= self.flu_fit_range[0]) * (data[:, 0] <= self.flu_fit_range[1])]

        if self.xaxis == 'Qz':
            self.qz = self.data_to_fit[:,0]
        elif self.xaxis == 'Sh':
            self.sh = self.data_to_fit[:,0]

        # for uncertainty calculation, this function is only used to get preparaed for multicore processing
        if uncertainty_calculation == True: return

        self.flu_result = lm.minimize(fl.flu2min, self.flu_par,
                                      args=((self.sh, self.qz), self.flucal_par),
                                      kws={'data':self.data_to_fit[:,1], 'eps':self.data_to_fit[:,2]}
                                      )
        self.flu_par = self.flu_result.params

        # for uncertainty calculation, the following is not needed

        tb = self.ui.fluparaTB
        tb.clear()
        tb.append(lm.fit_report(self.flu_result))

        self.updateUI(fresh=False)  # it has to be before updateFluCal()
        self.updateFluCal() # self.updateUI() has to be excucated before it reads parameters from GUI

    def saveFlu(self):
        if str(self.ui.flusaveCB.currentText())=='Save Fit':
            self.saveFluFitDig()
        elif str(self.ui.flusaveCB.currentText())=='Save Para':
            self.saveFluPara()

    def saveFluPara(self):

        self.updatePar()

        self.saveFileName = QFileDialog.getSaveFileName(caption='Save Fluorescence Fitting Parameters',
                                                        directory=self.directory,
                                                        filter='Par Files (*.par*;*_par.txt)')
        with open(self.saveFileName[0] + '_par.txt','w') as fid:
            try:
                try:
                    fid.write('Chi_Square\t' + format(self.flu_result.redchi, '.3f') + '\n')  # chisquare
                except:
                    fid.write('Chi_Square\tNA\n')

                fid.write('Fitting_Parameters\n')
                for p, u in self.ui_params.items():
                    fid.write(p + '\t\t' + format(float(u[0].text()),'.3e') + '\n')

                fid.write('\nSystem_Parameters\n')
                fid.write('Beam_Profile\t\t' + self.beam +'\n')
                for p, u in self.ui_syspar.items():
                    fid.write(p + '\t\t' + format(float(u.text()), '.4f') + '\n')

                print("Parameters saved!")

            except:
                print('Oops! Something went wrong, please check your parameters!')

    def loadFlu(self):
        if str(self.ui.fluloadCB.currentText())=='Load Para':
            self.loadFluPara()

    def loadFluPara(self):

        try:
            filename, _ = QFileDialog.getOpenFileName(caption='Select Parameter File to read',
                                                   directory=self.directory,
                                                   filter='Par Files (*.par*;*_par.txt)')
            self.directory = str(QFileInfo(filename).absolutePath())
            with open(str(filename)) as fid:
                fdata=fid.readlines()
        except IOError: # if the dialog is canceled.
            return
        # set ui values with loaded value
        line_num = 0
        line_type = 0 # 0: not a parameter line; 1: fitting parameter line; 2: system parameter line
        while True:
            try:
                line = fdata[line_num].split()
            except IndexError: # end of file
                break
            if line == []:
                line_type = 0
            else:
                if line[0] == 'Fitting_Parameters':
                    line_type = 1
                elif line[0] == 'Beam_Profile':
                    self.beam = line[1]
                    line_type = 2
                elif line_type == 1:
                    self.flu_par[line[0]].value = float(line[1])
                elif line_type == 2:
                    self.sys_par[line[0]] = float(line[1])
            line_num += 1

        # update parameter with new ui values
        self.updateUI(fresh=False)
        self.updateFluCal()

    def saveFluFitDig(self):

        Dialog=QDialog(self)
        self.uiflusavefit = uic.loadUi(resource_path(UI_path+'refsave.ui'), Dialog)
        self.uiflusavefit.label.setText('Save Fluorescence Fit/Calcualtion!')
        try:
            self.uiflusavefit.xminLE.setText(str(self.flu_simu_range[0]))
            self.uiflusavefit.xmaxLE.setText(str(self.flu_simu_range[1]))
        except:
            pass
        self.uiflusavefit.numpointLE.setText(str(200))

        self.uiflusavefit.cancelPB.clicked.connect(self.cancelSaveFluFit)
        self.uiflusavefit.okPB.clicked.connect(self.saveFluFit)

        self.uiflusavefit.show()

    def cancelSaveFluFit(self):
        self.uiflusavefit.close()
        self.flusavefitindex=0

    def saveFluFit(self):

        try:
            self.flusavefitindex = 1
            self.flunp = float(self.uiflusavefit.numpointLE.text())
            self.fluxmin = float(self.uiflusavefit.xminLE.text())
            self.fluxmax = float(self.uiflusavefit.xmaxLE.text())
            assert (self.fluxmin < self.fluxmax), "Maximum smaller than Minimum"

            self.saveFileName = QFileDialog.getSaveFileName(caption='Save Fluorescence Fit Data',
                                                            directory=self.directory)
            fname = self.saveFileName[0] + '_fit.txt'
            if self.xaxis == 'Qz':
                fit_to_save = self.flu[0,:,(1,2)].transpose()
            elif self.xaxis == 'Sh':
                fit_to_save = self.flu[:,0,(0,2)]
            np.savetxt(fname, fit_to_save, fmt='%.4e\t%.4e')

            self.flusavefitindex=0
            self.uiflusavefit.close()
        except AssertionError as AE:
            print("Error: {0}".format(AE))
        except IndexError as IE:
            pass # Ignore IndexError
        except:
            print("An error happens while saving fit file!")

    def debugErr(self):
        flufile = ['/Users/zhuzi/work/data/201912Dec/sh_sample03_318_50mMEu(NO3)3_s1h0.2_qz0.0015_flu.txt',
                   '/Users/zhuzi/work/data/201912Dec/sh_sample03_320_50mMEu(NO3)3_s1h0.2_qz0.0015_flu.txt',
                   '/Users/zhuzi/work/data/201912Dec/sh_sample03_494_50mMEu(NO3)3_s1h0.2_abs8_qz0.015_flu.txt']
        parfile = '/Users/zhuzi/work/data/201912Dec/sh_sample03_318_50mMEu(NO3)3_s1h0.2_qz0.0015_par.txt'
        self.flufiles = self.flufiles + flufile
        self.directory = str(QFileInfo(self.flufiles[0]).absolutePath())
        self.updateFluFile()

        with open(str(parfile)) as fid:
            fdata = fid.readlines()
        # set ui values with loaded value
        line_num = 0
        line_type = 0  # 0: not a parameter line; 1: fitting parameter line; 2: system parameter line
        while True:
            try:
                line = fdata[line_num].split()
            except IndexError:  # end of file
                break
            if line == []:
                line_type = 0
            else:
                if line[0] == 'Fitting_Parameters':
                    line_type = 1
                elif line[0] == 'Beam_Profile':
                    self.beam = line[1]
                    line_type = 2
                elif line_type == 1:
                    self.flu_par[line[0]].value = float(line[1])
                elif line_type == 2:
                    self.sys_par[line[0]] = float(line[1])
            line_num += 1
        # update parameter with new ui values
        self.updateUI()
        self.updateFluCal()

    def fluErrorInit(self):

        # choose the parameter for which the chisq is calculated
        self.fluerr_pname = [] # initialize # of the chosen parameters
        try:
            self.fluerr_pname = [p for p,u in self.ui_params.items() if u[1].isChecked()]
            if len(self.fluerr_pname) != 1:
                raise ValueError
            print("Calculating Chi-square for:", *self.fluerr_pname)
        except ValueError:
            print(" Did u pick the right number of parameters to fit?\n\n")
            # if multiple para's r checked, uncheck all and raise error
            for name in self.fluerr_pname:
                self.ui_params[name][1].setChecked(False)
            return

        self.uifluerr1=uic.loadUi(UI_path + 'err1.ui',QDialog(self))
        self.uifluerr1.label.setText('Uncertainty Calculation for Parameter:' + self.fluerr_pname[0])

        best_value = float(self.ui_params[self.fluerr_pname[0]][0].text())
        half_range_to_fit = abs(best_value*0.1)
        # the length of left and right half of range for the chosen values.
        self.uifluerr1.bestvalLE.setText(format(best_value, '.2e'))
        self.uifluerr1.leftLimitLE.setText(  # set left limit
            format(best_value - half_range_to_fit, '.2e'))
        self.uifluerr1.rightLimitLE.setText( # set right limit
            format(best_value + half_range_to_fit, '.2e'))

        self.uifluerr1.numIntervalLE.setText(format(10  ,'d'))

        # connect the pushbutton to next step
        # self.uifluerr1.cancelPB.clicked.connect(lambda x: self.uifluerr1.close())
        self.uifluerr1.cancelPB.clicked.connect(self.uifluerr1.close)
        self.uifluerr1.nextPB.clicked.connect(self.fluErrorPara)
        self.uifluerr1.show()

    def fluErrorPara(self):

        self.uifluerr1.close()
        # calculate a list of values the parameter should take where the chisq is calculated.
        self.fluerr_best_value = float(self.uifluerr1.bestvalLE.text())
        self.fluerr_left_limit = float(self.uifluerr1.leftLimitLE.text())
        self.fluerr_right_limit = float(self.uifluerr1.rightLimitLE.text())
        self.fluerr_num_points = int(self.uifluerr1.numIntervalLE.text())+1
        # append the fittted value for that parameter for displaying that
        # value in the chisq plot as the red dot.
        self.fluerr_fit_range = np.append(self.fluerr_best_value,
                                          np.linspace(self.fluerr_left_limit,
                                                      self.fluerr_right_limit,
                                                      self.fluerr_num_points))
        self.fluerr_chisq_list = np.zeros(self.fluerr_fit_range.shape)

        # automatically toggle the state of fiting and fixed parameters
        for p, u in self.ui_params.items():  u[1].toggle()

        # close the first dialog and open a new dialog
        self.uifluerr2 = uic.loadUi(UI_path + 'err2.ui', QDialog(self))
        self.uifluerr2.label.setText('Please check parameters to fit')
        self.uifluerr2.fluErrorProgress.setValue(0)
        # self.fluErroFit = myThread()
        # self.fluErrorFit.started.connect()
        # self.fluErrorFit.finished.connect(self.fluErrorResult)
        self.uifluerr2.fluErrorProgress.setMaximum(len(self.fluerr_fit_range))
        self.uifluerr2.cancelPB.clicked.connect(self.uifluerr2.close)
        self.uifluerr2.nextPB.clicked.connect(self.fluErrorFit)

        self.uifluerr2.show()

    def fluErrorFitSingleCore(self, q, i, flu_par,flucal_par,data_to_fit):
        fluerr_result = lm.minimize(fl.flu2min, flu_par,
                                    args=((self.sh, self.qz), flucal_par),
                                    kws={'data': data_to_fit[:, 1], 'eps': data_to_fit[:, 2]})
        q.put([i, fluerr_result.nfree, fluerr_result.redchi])



    def fluErrorFit(self):
        self.uifluerr2.label.setText('Calculating the uncertainty for ' + self.fluerr_pname[0])
        self.uifluerr2.nextPB.setEnabled(False) # unable the next push button
        # create a Parameter() object for fitting
        self.updatePar()  # update parameters to fit with GUI

        # make a copy of parameter just for uncertainty calculation
        fluerr_pname = self.fluerr_pname[0]
        fluerr_par = self.flu_par
        fluerr_cal_par = self.flucal_par
        fluerr_par[fluerr_pname].vary = False # make sure the chosen parameter does not vary

        # time the calculation
        start_time = time.time()
        # self.ErrorFit = myThread(self.fitFlu,errorbar=True,args={})



        # self.fitFlu(uncertainty_calculation=True)
        # fluErrorFitSingleCore_i = partial(fluErrorFitSingleCore2,
        #                                   value_list = self.fluerr_fit_range,
        #                                   sh = self.sh,
        #                                   qz = self.qz,
        #                                   pname = fluerr_pname,
        #                                   flu_par = fluerr_par,
        #                                   flucal_par = fluerr_cal_par,
        #                                   data_to_fit = self.data_to_fit)
        #
        #
        # results = multiCore(fluErrorFitSingleCore_i, range(len(self.fluerr_fit_range)))
        # for pp in results: print(pp)
        # print("Uncertainty calculation takes:", time.time() - start_time, "seconds")
        # for result in results:
        #     self.fluerr_chisq_list[result[0]] = result[3]
        # self.fluerr_nfree = results[-1][2]
        # print(self.fluerr_chisq_list)


        self.fitFlu(uncertainty_calculation=True)
        processes = []  # create a pool for processes
        q = multiprocessing.Queue()
        for i,value in enumerate(self.fluerr_fit_range):
            fluerr_par[fluerr_pname].value = value # change value of the chosen parameter
            p = multiprocessing.Process(target=self.fluErrorFitSingleCore,
                                        args=(q, i, fluerr_par,fluerr_cal_par,self.data_to_fit))
            processes.append(p)
            p.start()
        for process in processes:
            process.join()
        results = [q.get(True) for process in processes]

        for pp in results: print(pp)
        print("Uncertainty calculation takes:", time.time()-start_time, "seconds")
        for result in results:
            self.fluerr_chisq_list[result[0]] = result[2]
        self.fluerr_nfree = results[-1][1]
        print(self.fluerr_chisq_list)

        # # fit data and calculate chisq at each grid point
        # for i,para_value in enumerate(self.fluerr_fit_range):
        #     self.fluerr_parameters[self.fluerr_pname_to_fit].value = para_value
        #     fluresult=lm.minimize(self.flu2min, self.fluerr_parameters, args=(x,y,yerr))
        #     self.fluerr_chisq_list[i] = fluresult.redchi
        #     # update progress
        #
        # self.progressDialog.hide()

        self.uifluerr2.close()
        self.fluErrorResult()

    def fluErrorResult(self):
        # calculate the left/right error for the parameter
        funChisqFactor=interp1d(self.errorlist[:,0],self.errorlist[:,1],kind='cubic')
        chisq_factor = funChisqFactor(self.fluerr_nfree) # chisq_factor corresponding to degree of freedom
        idx_min_chisq = np.argmin(self.fluerr_chisq_list[1:]) + 1
        min_chisq = np.min(self.fluerr_chisq_list[1:])
        self.target_chisq = min_chisq * chisq_factor
        try: # interpolate function of left values against various chisq's
            funChisqListLeft = interp1d(self.fluerr_chisq_list[1:idx_min_chisq+1],
                                        self.fluerr_fit_range[1:idx_min_chisq+1],
                                        kind='linear')
            left_err = self.fluerr_best_value - funChisqListLeft(self.target_chisq)
            left_err_str = format(float(left_err),'.2e')
        except:
            left_err_str = "not found"
        try: # interpolate function of right values against various chisq's
            funChisqListRight = interp1d(self.fluerr_chisq_list[idx_min_chisq:],
                                         self.fluerr_fit_range[idx_min_chisq:],
                                         kind='linear')
            right_err = funChisqListRight(self.target_chisq) - self.fluerr_best_value
            right_err_str = format(float(right_err),'.2e')
        except:
            right_err_str = "not found"

        self.uifluerr3=uic.loadUi(UI_path + 'err3.ui',QDialog(self))
        self.uifluerr3.label.setText('Plot for Chi-square vs Parameter:'+self.fluerr_pname[0])
        self.uifluerr3.minchiLE.setText(format(min_chisq,'.2f'))
        self.uifluerr3.tarchiLE.setText(format(self.target_chisq,'.2f'))
        self.uifluerr3.lefterrLE.setText(left_err_str)
        self.uifluerr3.righterrLE.setText(right_err_str)
        self.uifluerr3.logyCB.stateChanged.connect(self.fluErrorPlot)
        self.uifluerr3.closePB.clicked.connect(lambda x: self.uifluerr3.close())
        self.uifluerr3.closePB.clicked.connect(self.uifluerr3.close)
        self.uifluerr3.savePB.clicked.connect(self.fluErrorSave)
        self.uifluerr3.show()
        self.fluErrorPlot()

    def fluErrorPlot(self):
        the_ax = self.uifluerr3.plotWidget.canvas.ax
        the_ax.clear()
        the_ax.set_xlabel(self.fluerr_pname[0])
        the_ax.set_ylabel('Chi-square')
        # check if y axis is logscale
        if self.uifluerr3.logyCB.checkState()!=0:
            the_ax.set_yscale('log')
        else:
            the_ax.set_yscale('linear')

        # plot the calculated chisq
        the_ax.plot(self.fluerr_fit_range[1:], self.fluerr_chisq_list[1:],
                    marker='o',ls='-')

        # plot the fitted parameter value and corresponding chisq
        the_ax.plot(self.fluerr_fit_range[0], self.fluerr_chisq_list[0],
                    marker='o',color='red')

        # plot the target chisq
        the_ax.plot(self.fluerr_fit_range[[1,-1]],
                    self.target_chisq * np.array([1,1]),
                    ls='-',color='green')

        self.uifluerr3.plotWidget.canvas.draw()

    def fluErrorSave(self):
        print("Save function to be released...")

    def insFluIon(self):  # add one ion in the subphase
        insrows=self.ui.flusubTW.selectionModel().selectedRows()
        insrows=[self.ui.flusubTW.row(self.ui.flusubTW.itemFromIndex(insrows[i])) for i in range(len(insrows))]
        if len(insrows)!=1:
            self.messageBox('Warning:: Only one row can be seleted!')
        else:
            self.ui.flusubTW.insertRow(insrows[0])
        for i in range(3):
            self.ui.flusubTW.setItem(insrows[0],i,QTableWidgetItem('Cl/2/1.80'.split('/')[i]))

    def rmFluIon(self): #remove one ion in the subphase
        rmrows=self.ui.flusubTW.selectionModel().selectedRows()
        removerows=[]
        for rmrow in rmrows:
            removerows.append(self.ui.flusubTW.row(self.ui.flusubTW.itemFromIndex(rmrow)))
            removerows.sort(reverse=True)
        if len(removerows)==0:
            self.messageBox('Warning:: No ion is selected!!')
        else:
            for i in range(len(removerows)):
                self.ui.flusubTW.removeRow(removerows[i])
