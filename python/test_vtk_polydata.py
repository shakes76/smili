'''
Test the Python SMILI bindings with polydata
'''
import sys
import vtk
from PySMILI import milxQtModel, milxQtFile
# from PySide2 import QtWidgets
from PySide6 import QtWidgets

if __name__ == "__main__":
    #setup Qt app
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()

    app.setOrganizationName("PySMILI")
    app.setApplicationName("PolyData Viewer")

    model = milxQtModel()

    points = vtk.vtkPoints()
    points.InsertNextPoint(0.0, 0.0, 0.0)  # Vertex 0
    points.InsertNextPoint(1.0, 0.0, 0.0)  # Vertex 1
    points.InsertNextPoint(0.5, 1.0, 0.0)  # Vertex 2

    # 2. Create a cell (in this case, a triangle) to connect the points
    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, 0) # Connect point 0 to vertex 0 of the triangle
    triangle.GetPointIds().SetId(1, 1) # Connect point 1 to vertex 1 of the triangle
    triangle.GetPointIds().SetId(2, 2) # Connect point 2 to vertex 2 of the triangle

    # 3. Create a vtkCellArray to store the cell(s)
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(triangle)

    # 4. Create vtkPolyData and set its points and cells
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    model = milxQtModel(mainWindow)
    model.SetInputData(polydata)
    model.generateModel()
    model.colourMapToJet()
    model.setWindowTitle("Model")

    mainWindow.setCentralWidget(model)
    mainWindow.resize(512, 512)
    mainWindow.show()

    app.exec()
    print("Done")
