import sys

from PyQt5.QtCore import pyqtRemoveInputHook
from PyQt5.QtWidgets import QApplication


from mainwindow3 import MainWindow

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

