'''
Test the Python SMILI bindings with polydata
'''
import sys
from PySMILI import vtkPoints, vtkPolyData, milxQtModel, milxQtFile
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

    # 4. Create vtkPolyData and set its points and cells
    polydata = vtkPolyData()
    fileIO.openModel(filename, polydata)

    model = milxQtModel(mainWindow)
    model.SetInput(polydata)
    model.generateModel()
    model.colourMapToJet()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)
    mainWindow.resize(512, 512)
    mainWindow.show()

    app.exec()
    print("Done")
