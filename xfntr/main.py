import sys
import os

#######################################################
# define an exception hook to prevent the app from crashing on exception
# reference: https://stackoverflow.com/questions/38020020/pyqt5-app-exits-on-error-where-pyqt4-app-would-not
sys._excepthook = sys.excepthook
def exception_hook(exctype, value, traceback):
    print(exctype, value, traceback)
    sys._excepthook(exctype, value, traceback)
sys.excepthook = exception_hook

# Use absolute path instead of relative path ('./') to avoid trouble when installed by pip
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)

from PyQt5.QtCore import pyqtRemoveInputHook
from PyQt5.QtWidgets import QApplication
from mainwindow import MainWindow

def main():
    # create application
    pyqtRemoveInputHook()
    app = QApplication(sys.argv)
    app.setApplicationName('15ID-C XFNTR Analyzer')

    # create widget
    w = MainWindow()
    w.setWindowTitle('15ID-C XFNTR Analyzer')
    # w.setWindowIcon(QIcon('logo.png'))
    w.show()

    # connection
    app.lastWindowClosed.connect(app.quit)

    # execute application
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()

