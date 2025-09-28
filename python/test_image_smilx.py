'''
Test the Python SMILI bindings with polydata
'''
import sys
import numpy as np
from PySMILI import milxQtImage, milxQtFile, milxQtMain
# from PySide2 import QtWidgets
from PySide6 import QtWidgets

# filename = "vase_1comp.vti"
filename = "brainweb165a10f17.mha" #ITKData

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = milxQtMain()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("Image Viewer")

    fileIO = milxQtFile()

    # load vtkImageData
    image = milxQtImage()
    fileIO.openImage(filename, image)
    image.generateImage()
    image.setWindowTitle("Image")
    image_shape = image.shape()
    print("image shape:", image_shape)

    mainWindow.addImage(image)

    mainWindow.show()

    app.exec()
    print("Done")
