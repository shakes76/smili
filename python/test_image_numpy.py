'''
Test the Python SMILI bindings with polydata
'''
import sys
import numpy as np
from PySMILI import milxQtImage
# from PySide2 import QtWidgets
from PySide6 import QtWidgets

filename = "vase_1comp.vti"

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("Image Viewer")

    # create 2D array
    N = 128
    # image_array = np.ones((N,N), dtype=np.float32)
    # image_array = np.random.randint(0, 256, size=(N, N, N), dtype=np.uint8)
    image_array = np.random.normal(loc=0, scale=1.0, size=(N, N, N)).astype(np.float32)
    print("image array shape:", image_array.shape)

    image = milxQtImage(mainWindow)
    image.setData(image_array.ravel(), N, N, N, 1.0)
    image.generateImage()
    image.setWindowTitle("Image")

    mainWindow.setCentralWidget(image)
    mainWindow.resize(512, 512)
    mainWindow.show()

    app.exec()
    print("Done")
