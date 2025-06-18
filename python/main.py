'''
Test the Python SMILI bindings
'''
import sys
from PySMILI import milxQtModel, milxQtFile
from PySide2 import QtWidgets

filename = "femur.vtk"

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("Model Viewer")

    fileIO = milxQtFile()
    model = milxQtModel(mainWindow)
    
    fileIO.openModel(filename, model)
    
    model.generateModel()
    model.colourMapToJet()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)
    mainWindow.resize(256, 256)
    mainWindow.show()

    app.exec_()
    print("Done")

