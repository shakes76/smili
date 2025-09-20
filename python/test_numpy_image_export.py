'''
Test the Python SMILI bindings with polydata
'''
import sys
import numpy as np
from PySMILI import milxQtImage, milxQtFile
# from PySide2 import QtWidgets
from PySide6 import QtWidgets
import matplotlib.pyplot as plt

filename = "vase_1comp.vti"

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)

    app.setOrganizationName("PySMILI")
    app.setApplicationName("Image Viewer")

    fileIO = milxQtFile()

    # load vtkImageData
    image = milxQtImage()
    fileIO.openImage(filename, image)
    image.imageInformation()
    image.generateImage()
    image.setWindowTitle("Image")

    image_array = image.getData()
    # image_array = image.get8BitData()
    image_shape = image.shape()
    print("image shape:", image_shape)
    image_array = np.array(image_array)
    image_array = image_array.reshape(image_shape)
    print("image array shape:", image_array.shape, "dtype:", image_array.dtype)

    plt.imshow(image_array[:,:,image_shape[2]//2])
    plt.title("Image")
    plt.show()

    print("Done")
