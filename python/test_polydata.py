'''
Test the Python SMILI bindings with polydata
'''
import sys
from PySMILI import vtkPolyData, milxQtModel, milxQtFile
from PySide2 import QtWidgets

filename = "femur.vtk"

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("PolyData Viewer")

    fileIO = milxQtFile()
    polydata = vtkPolyData()

    fileIO.openModel(filename, polydata)
  
    model = milxQtModel(mainWindow)
    model.SetInput(polydata)
    model.generateModel()
    model.colourMapToJet()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)
    mainWindow.resize(256, 256)
    mainWindow.show()

    app.exec_()
    print("Done")

