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

    model = milxQtModel(mainWindow)
    model.InsertNextPoint(1, 0, 0)
    model.InsertNextPoint(0, 1, 0)
    model.InsertNextPoint(0, 0, 1)
    model.generateModel()
    model.generatePointModel()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)
    mainWindow.resize(256, 256)
    mainWindow.show()

    app.exec_()
    print("Done")
