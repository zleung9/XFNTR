<<<<<<< HEAD
from PyQt4 import QtGui
import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib.backends.backend_qt4agg \
=======
# from PyQt4 import QtGui

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib.backends.backend_qt5agg \
>>>>>>> master
 import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

class MplCanvas(FigureCanvas):
	def __init__(self):
		self.fig = Figure()
		self.ax = self.fig.add_subplot(111)
		FigureCanvas.__init__(self, self.fig)
<<<<<<< HEAD
		FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)

class MplWidget(QtGui.QWidget):
	def __init__(self, parent = None):
		QtGui.QWidget.__init__(self, parent)
		self.canvas = MplCanvas()
		self.toolbar = NavigationToolbar(self.canvas, parent)
		self.vbl = QtGui.QVBoxLayout()
=======
		FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)

class MplWidget(QWidget):
	def __init__(self, parent = None):
		QWidget.__init__(self, parent)
		self.canvas = MplCanvas()
		self.toolbar = NavigationToolbar(self.canvas, parent)
		self.vbl = QVBoxLayout()
>>>>>>> master
		self.vbl.addWidget(self.canvas)
		self.vbl.addWidget(self.toolbar)
		self.setLayout(self.vbl)

