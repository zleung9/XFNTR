from PyQt4 import uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *
import numpy as np
import time
import matplotlib as mpl


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