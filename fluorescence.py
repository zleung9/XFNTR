from PyQt4 import uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *
from PyQt4 import QtCore, QtGui
from mplwidget import MplWidget
from matplotlib.widgets import MultiCursor
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy import special
from scipy.special import *
from scipy.integrate import quad
import pylab as pl
import numpy as np
import time
import os
import math
import cmath
import glob
import sys
import xr_ref as xr
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
import matplotlib as mpl
import periodictable 
from periodictable import *

def cancelErrCal():
    pass
    
def nextErrCal():
    pass

class uiErr1(QDialog):
    def __init__(self, para, num_interval=20): 
        self.iniUI()
        
    def iniUI(self):
        self.ui = uic.loadUi('err1.ui')
        print self.ui   
        self.ui.bestvalLE.setText(format(para.value,'.2e'))
        self.ui.leftranLE.setText(format(para.value*0.9, '.2e')) # left range: 0.1*value
        self.ui.rightranLE.setText(format(para.value*1.1, '.2e')) # right range: 0.1*value
        self.ui.numIntervalLE.setText(format(num_interval+1),'i') # number of intervals for the sample

        self.ui.CancelPB.clicked.connect(self.cancelErrCal)
        self.ui.NextPB.clicked.connect(self.nextErrCal)
        
        
    def paraFitRange(self,para,num_samples):

        best_value, left_limit, right_limit, num_samples = tuple(
            float(ui_err1.bestvalLE.text()),
            float(ui_err1.leftLimitLE.text()),
            float(ui_err1.rightLimitLE.text()),
            float(ui_err1.num_samples.text())
        )
        fit_range = np.linspace(lift_limit,right_limit,num_samples)