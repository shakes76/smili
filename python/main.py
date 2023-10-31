'''
Test the Python SMILI bindings
'''
import sys
from PySMILI import milxQtWindow, milxQtModel, milxQtFile
from PySide2 import QtWidgets

filename = "/home/uqscha22/data/SyngoData_focused/AverageNonrigid_R_femur.vtk"

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
    #model.setWindowTitle("Model")
    model.resize(256, 256)

    mainWindow.setCentralWidget(model)
    mainWindow.show()

    app.exec_()
    print("Done")

