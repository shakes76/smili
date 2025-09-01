'''
Test the Python SMILI bindings with polydata
'''
import sys
from PySMILI import milxQtModel, milxQtFile
# from PySide2 import QtWidgets
from PySide6 import QtWidgets

filename = "Bunny.vtp"

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("PolyData Viewer")

    fileIO = milxQtFile()
    model = milxQtModel()

    fileIO.openModel(filename, model)

    model.generateModel()
    model.colourMapToJet()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)
    mainWindow.resize(512, 512)
    mainWindow.show()

    app.exec()
    print("Done")
