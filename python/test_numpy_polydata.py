'''
Test the Python SMILI bindings with polydata
'''
import sys
import numpy as np
from PySMILI import milxQtModel
# from PySide2 import QtWidgets
from PySide6 import QtWidgets

filename = "Bunny.vtp"

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("PolyData Viewer")

    # Create a 3x3 NumPy array for 3 points.
    # Need float64 for double cast in SMILI
    points_array = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ], dtype=np.float64)

    model = milxQtModel(mainWindow)
    model.SetInput(points_array)
    model.generatePointModel()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)

    mainWindow.resize(512, 512)
    mainWindow.show()

    app.exec()
    print("Done")
