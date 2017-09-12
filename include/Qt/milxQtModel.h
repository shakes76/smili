/*=========================================================================
  The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO)
  ABN 41 687 119 230.
  All rights reserved.

  Licensed under the CSIRO BSD 3-Clause License
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License in the file LICENSE.md or at

  https://stash.csiro.au/projects/SMILI/repos/smili/browse/license.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
=========================================================================*/
#ifndef MILXQTMODEL_H
#define MILXQTMODEL_H

#include "milxQtRenderWindow.h"
#include <QFormLayout>

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkAlgorithmOutput.h>
#include <vtkTransform.h>
#include <vtkAbstractTransform.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkOutlineFilter.h>
#include <vtkCubeAxesActor.h>
//Labelling
#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkImageActor.h>
//Rendering
#include <vtkPolyDataMapper.h>
#include <vtkLODActor.h>

#include "milxModel.h"

static const int defaultChannelValue = 192;
static const float defaultColour = defaultChannelValue/255.0; //HTML Gray

/*!
    \class milxQtModel
    \brief This class represents the MILX Qt Model/Mesh Display object using VTK.
    \author Shekhar S. Chandra, 2013

    The class displays models from polygonal values using OpenGL via the VTK library. The rendering is encapsulated within
    a QVTK widget. It is important to generateModel(), after setting the polydata or when loading the model in a different way.
    This ensures everything needed to display has been updated.

    Controls:
        - Left Mouse (hold): adjust the colour level.
        - Wheel: zoom the camera.
        - Middle Mouse (hold): translate image.
        - Right Mouse: context menu and zoom toggle.
        - Keypress f: fly to the picked point.
        - Keypress Shift+r: reset the camera view along the current view direction. Centers the actors and moves the camera so that all actors are visible.
        - Keypress r: reset the colour level to defaults.

    Usage Examples:
    Loading a model
    \code
    ///Load model
    const QString filename = "model.vtk";
    QPointer<milxQtFile> reader = new milxQtFile;
    QPointer<milxQtModel> model = new milxQtModel;
    bool success = reader->openModel(filename, model);

    //If loaded, set attributes and generate.
    if(success)
    {
        model->setName(filename);
        model->setConsole(console);
        model->setDefaultView(defaultViewBox->currentIndex());
        model->generateModel();
    }
    \endcode
    Displaying a sphere
    \code
    ///Define Test Polygon
    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetThetaResolution(32);
        sphere->SetPhiResolution(32);
        sphere->SetRadius(10.0);
        sphere->Update();

    model->setName("Sphere");
        model->SetInput(sphere->GetOutput());
        model->generateModel(0.5, 0.5, 0.5); //Grey model
        model->show();
    \endcode
    Delaunay Graph from points
    \code
    QPointer< milxQtModel > model = new milxQtModel;
        model->setName("Original Model");
        model->SetPoints(modelPoints);
        model->generateDelaunayGraph();
        model->generateLabels();
        model->show();
    \endcode
    Point Model from points
    \code
    QPointer< milxQtModel > centroidModel = new milxQtModel;
        centroidModel->SetPoints(centroidPoints);
        centroidModel->generatePointModel(0.25, 0.5, 0.5, 1.0); //!< Graph in grey
    \endcode
*/
class MILXQT_EXPORT milxQtModel : public milxQtRenderWindow
{
    Q_OBJECT

public:
    /*!
        \fn milxQtModel::milxQtModel(QWidget *parent = 0, bool contextSystem = true)
        \brief The standard constructor
    */
    milxQtModel(QWidget *theParent = 0, bool contextSystem = true);
    /*!
        \fn milxQtModel::~milxQtModel()
        \brief The standard destructor
    */
    virtual ~milxQtModel();

    /*!
        \fn milxQtModel::strippedNamePrefix()
        \brief Returns the stripped (path removed) name of the data with "Model" prefix.
    */
    inline QString strippedNamePrefix()
    {
        return prefix + QFileInfo(name).fileName();
    }

    /*!
        \fn milxQtModel::AddInput(vtkSmartPointer<vtkPolyData> mesh)
        \brief Assigns and coninually appends meshes provided to the class, preparing for display. Call generateModel() and then show() to display.
    */
    void AddInput(vtkSmartPointer<vtkPolyData> mesh);
    /*!
        \fn milxQtModel::AddArray(vtkSmartPointer<vtkDataArray> array)
        \brief Adds an array to the model, appending it to the current list of arrays.
    */
    void AddArray(vtkSmartPointer<vtkDataArray> array);
    /*!
        \fn milxQtModel::SetInput(vtkSmartPointer<vtkPolyData> mesh)
        \brief Assigns the mesh provided to the class, preparing for display. Call generateModel() and then show() to display.
    */
    void SetInput(vtkSmartPointer<vtkPolyData> mesh);
    /*!
        \fn milxQtModel::SetInputPointSet(vtkSmartPointer<vtkPointSet> mesh)
        \brief Assigns the pointset mesh provided to the class, preparing for display. No need to call generateModel(), just show() to display.
    */
    //void SetInputPointSet(vtkSmartPointer<vtkPointSet> mesh);
    /*!
        \fn milxQtModel::SetNumberOfInputs(const int inputs)
        \brief Number of inputs to expect when using SetInput() member with indices.
    */
//    void SetNumberOfInputs(const int inputs);
    /*!
        \fn milxQtModel::SetInput(const int index, vtkSmartPointer<vtkPolyData> mesh)
        \brief Assigns the mesh provided to the class to index, preparing for display. Call generateModel() and then show() to display. UNTESTED

        Ensure the number of inputs (via SetNumberOfInputs()) is set before calling this member.
    */
//    void SetInput(const int index, vtkSmartPointer<vtkPolyData> mesh);
    /*!
        \fn milxQtModel::SetPoints(vtkSmartPointer<vtkPoints> modelPoints)
        \brief Sets the points for the model to be generated. Must pass a vtkPoints objects, which is easy to use.
    */
    void SetPoints(vtkSmartPointer<vtkPoints> modelPoints);
    /*!
        \fn milxQtModel::SetPolys(vtkSmartPointer<vtkCellArray> modelPolys)
        \brief Sets the polygons for the model to be generated. Must pass a vtkCellArray objects, which is easy to use.
    */
    void SetPolys(vtkSmartPointer<vtkCellArray> modelPolys);
    /*!
        \fn milxQtModel::SetScalars(vtkSmartPointer<vtkDataArray> modelScalars)
        \brief Sets the scalars for the model to be generated. Must pass a vtkDataArray or its derivatives (such as a vtkFloatArray) objects, which is easy to use.
    */
    void SetScalars(vtkSmartPointer<vtkDataArray> modelScalars);
    /*!
        \fn milxQtModel::SetActiveScalars(vtkSmartPointer<vtkDataArray> modelScalars)
        \brief Sets the array of name as scalars for the model to be generated.
    */
    inline void SetActiveScalars(std::string nameOfArray)
    { model.SetActiveScalars(nameOfArray);  }
    /*!
        \fn milxQtModel::RemoveScalars()
        \brief Removes the scalars of the model.
    */
    void RemoveScalars();
    /*!
        \fn milxQtModel::SetVectors(vtkSmartPointer<vtkDataArray> modelVectors)
        \brief Sets the vectors for the model to be generated. Must pass a vtkFloatArray objects, which is easy to use.
    */
    void SetVectors(vtkSmartPointer<vtkDataArray> modelVectors);
    /*!
        \fn milxQtModel::SetTransform(vtkSmartPointer<vtkTransform> transform)
        \brief Sets the transform for the model that will be generated. Must pass a vtkTransform objects, which is easy to use.
    */
    void SetTransform(vtkSmartPointer<vtkTransform> transform);
    /*!
        \fn milxQtModel::SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph)
        \brief Assigns the graph provided to the class, preparing for display. Call generateModel() and then show() to display the graph in 3D.
    */
    void SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph);
    /*!
        \fn milxQtModel::SetOpacity(double opacity)
        \brief Adjusts the opacity of the model in the display. 1.0 is totally opaque and 0.0 is completely transparent.

        WARNING: The model should have been generated for this member to have any effect.
    */
    inline void SetOpacity(double opacity)
    {
        if(modelled) modelActor->GetProperty()->SetOpacity(opacity);
    }

    //VTK Interfacing
    /*!
        \fn milxQtModel::GetPolyDataInput()
        \brief Returns the original mesh data object (PolyData) used internally VTK style.

        This is different to VTK in the sense that it may return a PolyData set by SetGraph() or SetInput().
    */
    vtkSmartPointer<vtkPolyData> GetPolyDataInput();
    /*!
        \fn milxQtModel::GetOutput()
        \brief Returns the mesh data object (PolyData) used internally VTK style.
    */
    vtkSmartPointer<vtkPolyData> GetOutput();
#if VTK_MAJOR_VERSION <=5
    /*!
        \fn milxQtModel::GetOutputPort()
        \brief Returns the connection to the mesh data object (PolyData) used internally VTK style. Call after generate functions. EXPERIMENTAL
    */
    inline vtkSmartPointer<vtkAlgorithmOutput> GetOutputPort()
    {
        return model.GetProducerPort();
    }
#endif
    /*!
        \fn milxQtModel::GetActor()
        \brief Returns the VTK actor used.
    */
    inline vtkSmartPointer<vtkLODActor> GetActor()
    {
        return modelActor;
    }
    /*!
        \fn milxQtModel::GetActor2D()
        \brief Returns the VTK actor 2D used.
    */
    inline vtkSmartPointer<vtkActor2D> GetActor2D()
    {
        return modelLabelsActor;
    }
    /*!
        \fn milxQtModel::GetMapper()
        \brief Returns the VTK mapper used.
    */
    inline vtkSmartPointer<vtkPolyDataMapper> GetMapper()
    {
        return modelMapper;
    }
    /*!
        \fn milxQtModel::GetNumberOfPoints()
        \brief Returns the total number of points of the model currently held.
    */
    inline vtkIdType GetNumberOfPoints()
    {
        return model.GetNumberOfPoints();
    }
    /*!
        \fn milxQtModel::GetNumberOfArrays()
        \brief Returns the total number of arrays in the model.
    */
    inline vtkIdType GetNumberOfArrays()
    {
        return model.GetNumberOfArrays();
    }
    /*!
        \fn milxQtModel::GetPoints()
        \brief Returns the points of the model currently held.
    */
    inline vtkSmartPointer<vtkPoints> GetPoints()
    {
        return model.GetPoints();
    }
    /*!
        \fn milxQtModel::GetScalars()
        \brief Returns the scalars of the model currently held.
    */
    inline vtkSmartPointer<vtkDataArray> GetScalars()
    {
        return model.GetScalars();
    }
    /*!
        \fn milxQtModel::GetVectors()
        \brief Returns the vectors of the model currently held.
    */
    inline vtkSmartPointer<vtkDataArray> GetVectors()
    {
        return model.GetVectors();
    }
    /*!
        \fn milxQtModel::GetNormals()
        \brief Returns the normals of the model currently held.
    */
    inline vtkSmartPointer<vtkDataArray> GetNormals()
    {
        return model.GetNormals();
    }
    /*!
        \fn milxQtModel::GetTensors()
        \brief Returns the tensors (3x3 matrix per vertex) of the model currently held.
    */
    inline vtkSmartPointer<vtkDataArray> GetTensors()
    {
        return model.GetTensors();
    }
    /*!
      \brief Get the polydata data, return the vtkImageData object downcast to vtkDataSet, useful for getting scalar range etc.
    */
    inline virtual vtkDataSet* GetDataSet()
    {
        return model.Result();
    }
    /*!
        \fn milxQtModel::SetScalarRange(double range[2])
        \brief Set the scalar range of the model currently held.
    */
    inline void SetScalarRange(double range[2])
    {
        modelMapper->SetScalarRange(range);
        modelMapper->SetLookupTable(lookupTable);
        modelMapper->Update();
    }
    /*!
        \fn milxQtModel::GetScalarRange(double range[2])
        \brief Get the scalar range of the model currently held.
    */
    inline void GetScalarRange(double range[2])
    {
        model.Result()->GetScalarRange(range);
    }
    /*!
        \fn milxQtModel::resetScalarRange()
        \brief Restores the scalar range of the model currently held to current max range for display.

        Affects display only.
    */
    inline void resetScalarRange()
    {
        if(GetScalars())
            modelMapper->SetScalarRange( model.Result()->GetScalarRange() );
    }
    /*!
        \fn milxQtModel::GetPoint(vtkIdType id)
        \brief Returns the point with id of the model currently held.
    */
    inline double * GetPoint(vtkIdType id)
    {
        return model.Result()->GetPoints()->GetPoint(id);
    }
    /*!
        \fn milxQtModel::GetTransform()
        \brief Returns the transform of the model currently held.
    */
    inline vtkSmartPointer<vtkAbstractTransform> GetTransform()
    {
        return vtkAbstractTransform::SafeDownCast(model.GetTransform());
    }
    /*!
        \fn milxQtModel::GetOpacity()
        \brief Returns the opacity of the model in the display. 1.0 is totally opaque and 0.0 is completely transparent.

        The model should have been generated for this member to have any effect.
    */
    inline double GetOpacity()
    {
        if(modelled) return modelActor->GetProperty()->GetOpacity();
        return 1.0;
    }
    /*!
        \fn milxQtModel::Update()
        \brief Update the model currently held.
    */
    inline void Update()
    {
        model.Update();
    }

    //Rendering
    /**
        \fn milxQtModel::setLargeDataSetMode(bool large)
        \brief Set display mode to memory conserving and fast rendering (in most cases).
    */
    inline void setLargeDataSetMode(bool large)
    {
        largeMode = large;
    }
    /**
        \fn milxQtModel::ImmediateModeRenderingOn()
        \brief Improve rendering performance for large datasets. Assumes generateModel() has already been called.
    */
    inline void ImmediateModeRenderingOn()
    {
        if(modelled) modelMapper->ImmediateModeRenderingOn();
    }

    //Operators
    /*!
        \fn milxQtModel::operator=(const milxQtModel &operand)
        \brief Assignment operator for models. Only copies the data as given by GetOutput() of operand. EXPERIMENTAL.
    */
    milxQtModel& operator=(const milxQtModel &operand);

public slots:
    //Display
    /*!
        \fn milxQtModel::refresh()
        \brief Refresh the display of the model.
    */
    void refresh();
    /*!
        \fn milxQtModel::reset()
        \brief Reset the display of the model, includes camera and regeneration.
    */
    inline void reset()
    {
        generateModel(colourRed, colourGreen, colourBlue);
    }
    /*!
        \fn milxQtModel::toggleInterpolation(bool quiet = false)
        \brief Toggles the interpolation of the mesh between Phong and Gouraud. Default is Gouraud.
    */
    void toggleInterpolation(bool quiet = false);
    inline void interpolateDisplay(bool quiet = false)
    {   toggleInterpolation(quiet);    }
    inline void disableInterpolateDisplay(bool quiet = false)
    {   interpAct->setChecked(false);   interpolateDisplay(quiet);   }
    inline void enableInterpolateDisplay(bool quiet = false)
    {   interpAct->setChecked(true);   interpolateDisplay(quiet);   }
    /*!
        \fn milxQtModel::toggleSpecular(bool quiet = false)
        \brief Toggles the specular or shininess of the mesh between Flat and Shiny. Default is Shiny.
    */
    void toggleSpecular(bool quiet = false);
    inline void specularDisplay(bool quiet = false)
    {   toggleSpecular(quiet);    }
    inline void disableSpecularDisplay(bool quiet = false)
    {   specularAct->setChecked(false);   specularDisplay(quiet);    }
    inline void enableSpecularDisplay(bool quiet = false)
    {   specularAct->setChecked(true);   specularDisplay(quiet);    }
    /*!
        \fn milxQtModel::copyToContextMenu(QMenu *copyMenu)
        \brief Copies the menu, by duplicating the entries, to the context menu. Connections are assumed to be made before hand.

        This virtual implementation accounts for menus with custom connections.
    */
    virtual void copyToContextMenu(QMenu *copyMenu);

    //Computation
    /*!
        \fn milxQtModel::centroid()
        \brief Computes and returns the centroid of the model currently held.

        Does not re-compute unless the model has been changed.
    */
    coordinate& centroid();
    /*!
        \fn milxQtModel::centroidSize(bool average = false)
        \brief Computes and returns the centroid size of the model currently held which can be used for rescaling.
    */
    double centroidSize(bool average = false);
    /*!
        \fn milxQtModel::covarianceMatrix()
        \brief Computes and returns the covariance matrix of the model currently held.
    */
    vnl_matrix<double>& covarianceMatrix();

    //Filters
    /*!
        \fn milxQtModel::undoProcessing()
        \brief Undo/Redo the processing last done. Only applies to processing such as clean() etc.
    */
    void undoProcessing();
    /*!
        \fn milxQtModel::clean()
        \brief Clean the model of duplicate points by merging appropriately.
    */
    void clean();
    /*!
        \fn milxQtModel::triangulate()
        \brief Triangulate the model cells.
    */
    void triangulate();
    /*!
        \fn milxQtModel::decimate(double factor = 0.0)
        \brief Reduce the polygon count (decimate) the model. Factor to reduce by is acquired by input a dialog.
    */
    void decimate(double factor = 0.0);
    /*!
        \fn milxQtModel::quadricDecimate(double factor = 0.0)
        \brief Reduce the polygon count (decimate) the model using Quadric Decimation. The final result is the best simplification of the mesh, preserving topology.
    */
    void quadricDecimate(double factor = 0.0);
    /*!
        \fn milxQtModel::clusterDecimate()
        \brief Reduce the polygon count (decimate) the model using Quadric Clustering Decimation. The final result is the best simplification of the mesh, preserving topology.
    */
    void clusterDecimate();
    /*!
        \fn milxQtModel::renameScalars(QString newName = "")
        \brief Rename the scalars on the mesh.
    */
    void renameScalars(QString newName = "");
    /*!
        \fn milxQtModel::threshold(double lowValue = 0.0, double upValue = 0.0)
        \brief Clip the mesh based on scalar values on the mesh.
    */
    void threshold(double lowValue = 0.0, double upValue = 0.0);
    /*!
        \fn milxQtModel::thresholdScalars(double lowValue = 0.0, double upValue = 0.0, double outsideVal = 0.0)
        \brief Threshold scalar values on the mesh. Geometry and topology is unchanged.
    */
    void thresholdScalars(double lowValue = 0.0, double upValue = 0.0, double outsideVal = 0.0);
    /*!
        \fn milxQtModel::thresholdScalarsBinary(double lowValue = 0.0, double upValue = 0.0, double insideVal = 1.0, double outsideVal = 0.0)
        \brief Binary threshold scalar values on the mesh. Geometry and topology is unchanged, and scalars become two-valued, inside and outside.
    */
    void thresholdScalarsBinary(double lowValue = 0.0, double upValue = 0.0, double insideVal = 1.0, double outsideVal = 0.0);
    /*!
        \fn milxQtModel::maskScalars(QString filename = "")
        \brief Mask scalar values on the mesh using scalars on another mesh. Assumes point correspondence between meshes.
    */
    void maskScalars(QString filename = "");
    /*!
        \fn milxQtModel::gradient()
        \brief Compute the gradient of the scalar field on the model.
    */
    void gradient();
    /*!
        \fn milxQtModel::smooth(int iterations = 0)
        \brief Smooth the model using the standard Laplacian approach. Iterations to smooth by is acquired by input a dialog.
    */
    void smooth(int iterations = 0);
    /*!
        \fn milxQtModel::smoothSinc(int iterations = 0)
        \brief Smooth the model using the Windowed Sinc approach of Taubin et al. Iterations to smooth by is acquired by input a dialog.

        A pass band of 0.1 is used internally with all settings at default for the vtkWindowedSincPolyDataFilter filter.
        Paper: "Optimal Surface Smoothing as Filter Design" G. Taubin, T. Zhang and G. Golub
    */
    void smoothSinc(int iterations = 0);
    /*!
        \fn milxQtModel::curvature()
        \brief Compute the mean curvature of mesh and place the curvature as scalars on mesh.
    */
    void curvature();
    /*!
        \fn milxQtModel::transform(QString filename = "", bool inverse = false)
        \brief Transform the model using a loaded affine transform file. File is asked for within this member.
    */
    void transform(QString filename = "", bool inverse = false);
    /*!
        \fn milxQtModel::loadScalars(QString filename = "")
        \brief Loads the scalars from another model, where both models must have the same number of points.
    */
    void loadScalars(QString filename = "");
    /*!
        \fn milxQtModel::removeScalars()
        \brief Removes all the scalars from the model.
    */
    void removeScalars();
    /*!
        \fn milxQtModel::append(milxQtModel *mdl)
        \brief Appends the provided model to current model.
    */
    inline void append(milxQtModel *mdl)
    {
        AddInput(mdl->GetOutput());
    }
    /*!
        \fn milxQtModel::modelInfo()
        \brief Computes and prints the centroid, centroid size and covariance matrix to standard out.
    */
    void modelInfo();
    /*!
        \fn milxQtModel::matchInfo(QString filename = "", bool rescale = true, double factor = 1.0)
        \brief Loads another model and matches its centroid and (optionally) scales this model to match.

        The (optional) factor is how much (or the factor) of the rescaling to apply, e.g. 0.9 would be 90% of the rescaling.
    */
    void matchInfo(QString filename = "", bool rescale = true, double factor = 1.0);
    /*!
        \fn milxQtModel::registerICP(QString filename = "", bool similarity = false)
        \brief Register or align current mesh to another mesh using Iterative Closest Point algorithm, whose filename is provided. If not provided, then Open File dialog will apear.
    */
    void registerICP(QString filename = "", bool similarity = false);
    /*!
        \fn milxQtModel::registerLandmarks(QString filename = "", bool similarity = false)
        \brief Register or align current mesh to another mesh using the landmark transform, whose filename is provided. If not provided, then Open File dialog will apear.

        Warning: This registration method assumes meshes have the same number of points and point to point correspondence.
    */
    void registerLandmarks(QString filename = "", bool similarity = false);
    /*!
        \fn milxQtModel::voxelise()
        \brief Convert surface to an image volume.
    */
    void voxelise();
    /*!
        \fn milxQtModel::contour()
        \brief Draw contour interactively
    */
    virtual void contour();
    /*!
        \fn milxQtModel::selectPointsInContour()
        \brief Select model points inside contour and compute distance map.
    */
    virtual void selectPointsInContour();
    /*!
        \fn milxQtModel::weightGaussianFromContour(float stddev = -1.0)
        \brief Select model points inside contour, compute distance map and Gaussian weight the result.
    */
    virtual void weightGaussianFromContour(float stddev = -1.0, float clampValue = -1.0);
    /*!
        \fn milxQtModel::texture()
        \brief Emit texture of the model.
    */
    void texture();
    /*!
        \fn milxQtModel::changeColour(float red = 0, float green = 0, float blue = 0)
        \brief Changes the colour of the model if no scalars are present.
    */
    void changeColour(float red = -1, float green = -1, float blue = -1);
    /*!
        \fn milxQtModel::changeOpacity(float opacity = 1.0)
        \brief Changes the opacity/transparency of the model. You should use the changeColour() member and set its alpha channel.

        Action for this member has been removed because changeColour() provides this functionality.
    */
    void changeOpacity(float opacity = -1);
    /*!
        \fn milxQtModel::normals(const bool turnOn = false)
        \brief Shows the point normals of the model

        Argument provided to bypass GUI/checked action
    */
    void normals(const bool turnOn = false);
    /*!
        \fn milxQtModel::rotate(bool xAxis = false, bool yAxis = false, bool zAxis = false, float angle = 90)
        \brief Rotate model along any axis, useful for matching orientation with ITK etc.
    */
    void rotate(bool xAxis = false, bool yAxis = false, bool zAxis = false, float angle = 90);
    /*!
        \fn milxQtModel::flip(bool xAxis = false, bool yAxis = false, bool zAxis = false)
        \brief Flip model along any axis, useful for matching coordinate systems with ITK or flipping left to right.
    */
    void flip(bool xAxis = false, bool yAxis = false, bool zAxis = false);
    /*!
        \fn milxQtModel::background(bool white = false)
        \brief Changes background to white.
    */
    void background(bool white = false);
    /*!
        \fn milxQtModel::enableOutline(vtkDataObject *dataOutline = NULL)
        \brief Display outline box of data.
    */
    void enableOutline(vtkDataObject *dataOutline = NULL);
    /*!
        \fn milxQtModel::disableOutline()
        \brief Disables the outline box display.
    */
    void disableOutline();
    /*!
        \fn milxQtModel::outlineDisplay()
        \brief Toggles the outline box display.
    */
    void outlineDisplay();
    /*!
        \fn milxQtModel::enableCubeAxes(double *range = NULL, double *bounds = NULL)
        \brief Display cube (plot) axes of data.
    */
    void enableCubeAxes(double *range = NULL, double *bounds = NULL);
    /*!
        \fn milxQtModel::disableCubeAxes()
        \brief Disables the cube axes display.
    */
    void disableCubeAxes();
    /*!
        \fn milxQtModel::cubeAxesDisplay(double *range = NULL)
        \brief Toggles the cube (plot) axes display.
    */
    void cubeAxesDisplay(double *range = NULL);
    /*!
        \fn milxQtModel::removeOverlays()
        \brief Removes all actors that are overlaid with the model.
    */
    void removeOverlays();
    /*!
        \fn milxQtModel::enableScale(QString title = "", const bool quiet = false, double minRange = 0.0, double maxRange = 0.0, int noOfLabels = 3)
        \brief Enable scale bar display with the title provided.

        Quiet Boolean is to prevent possible popups to ask user parameters.
    */
    virtual void enableScale(QString title = "", const bool quiet = false, double minRange = 0.0, double maxRange = 0.0, int noOfLabels = 3);
    /*!
        \fn milxQtModel::scaleDisplay(const bool forceDisplay = false)
        \brief Toggles the scale bar display.

        forceDisplay Boolean is to overide possible previous settings and display bar.
    */
    virtual void scaleDisplay(const bool forceDisplay = false);
    /*!
        \fn milxQtModel::outputScalars(QString filename)
        \brief Output the scalars of model to file.
    */
    void outputScalars(QString filename = "");
    /*!
        \fn milxQtModel::viewToXYPlane()
        \brief Change view to xy-plane.
    */
    virtual void viewToXYPlane();
    /*!
        \fn milxQtModel::viewToZXPlane()
        \brief Change view to zx-plane.
    */
    virtual void viewToZXPlane();
    /*!
        \fn milxQtModel::viewToZYPlane()
        \brief Change view to zy-plane.
    */
    virtual void viewToZYPlane();
    /*!
        \fn milxQtModel::showArray(QAction * action)
        \brief Change scalars on model to one given by string or action.
    */
    void showArray(QAction * action);
    void showArray(const QString arrayName);
    /*!
        \fn milxQtModel::updateLookupTable()
        \brief Sets the necessary LUTs to model view.
    */
    virtual void updateLookupTable();

    /*!
        \fn milxQtModel::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);

    //Generators
    /*!
        \fn milxQtModel::generateDelaunayGraph(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the Delaunay 3D graph for the dataset. Arguments provided are for the colours of the edges of the graph.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the graph.
    */
    void generateDelaunayGraph(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateDelaunay2DTriangulation()
        \brief Generates the Delaunay 2D Triangulation (resulting in triangulated cells) for the model.
    */
    void generateDelaunay2DTriangulation();
    /*!
        \fn milxQtModel::generateDelaunayTriangulation(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the Delaunay 3D Triangulation (resulting in tetrahdral cells) for the dataset. Arguments provided are for the colours of the cells of the model.
    */
    void generateDelaunayTriangulation(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateVertices(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the vertices for the dataset. Arguments provided are for the colours of the vertices.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the vertices.
    */
    void generateVertices(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateVerticesAs(const GlyphType glyphType = Circle, float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the vertices as 2D circles for the dataset. Arguments provided are for the colours of the vertices.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the vertices.
    */
    void generateVerticesAs(const GlyphType glyphType = Circle, float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateTubes(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the tubes for the dataset. Arguments provided are for the colours of the tubes.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the graph.
    */
    void generateTubes(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateEdges(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the tubes for the dataset. Arguments provided are for the colours of the tubes. Same as generateTubes, alias for convenience.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the graph.
    */
    inline void generateEdges(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
    {
        generateTubes(red, blue, green);
    }
    /*!
        \fn milxQtModel::generatePointModel(double newScale = 1.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates a point model, i.e. glyphs at each point for the dataset. Arguments provided are for the colours of the edges of the graph.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the glyphs. Use generatePoints() to avoid changing the pipeline.

        You should have set the points (using SetPoints()) before call this member.
    */
    void generatePointModel(double newScale = 1.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateSampledPoints(float distance = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates points sampled "distance" apart on the model with glyphs at each point.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the glyphs.

        You should have set the points (using SetPoints()) before call this member.
    */
    void generateSampledPoints(float distance = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateVectorField(double newScale = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates a vector field model, i.e. arrows at each point for the dataset. Arguments provided are for the colours of the vectors.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the vectors.

        You should have set the points and vectors (using SetPoints() and SetVectors()) before call this member.
    */
    void generateVectorField(double newScale = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateTensorField(double newScale = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates a tensor field model, i.e. ellipses at each point for the dataset. Arguments provided are for the colours of the ellipses.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the ellipses.

        You should have set the points and tensors (using SetPoints() and SetTensors()) before call this member.
    */
    void generateTensorField(double newScale = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateStreamLines(vtkSmartPointer<vtkImageData> seeds, float timestep = 0.0)
        \brief Generates stream lines i.e. integrates a vector field starting at seed points.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the lines.

        You should have set the seeds as the model/surface (via SetInput()) before call this member.
    */
    void generateStreamLines(vtkSmartPointer<vtkImageData> seeds, float timestep = 0.0);
    /*!
        \fn milxQtModel::generateHedgehog(double newScale = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates a hedgehog model, i.e. lines at each point for the dataset whose lengths are given by the vector magnitudes. Arguments provided are for the scale and colours of the lines.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the vectors.

        You should have set the points and vectors (using SetPoints() and SetVectors()) before call this member.
    */
    void generateHedgehog(double newScale = 0.0, float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateNormals(int pointNormals = 0)
        \brief Generates normals of the model, with point normals being default.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the normals. Use normals() member to see these normals.
    */
    void generateNormals(int pointNormals = 0);
    /*!
        \fn milxQtModel::generateModel(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the model so that its ready for display. It requires that data has been set or assigned already.

        By default, the display is a surface. Use generatePoints() or generateWireframe() before this call to generate those displays of the model.
    */
    void generateModel(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateLabels()
        \brief Generates the model vertices with labels. It requires that data has been set or assigned already.
    */
    void generateLabels();
    /*!
        \fn milxQtModel::generatePointIDsScalars()
        \brief Generates scalars for the vertices in the model based on their point IDs. It requires that model has been set or assigned already.
    */
    void generatePointIDsScalars();
    /*!
        \fn milxQtModel::generatePoints(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates points display of the model, i.e. points at each point for the dataset. Arguments provided are for the colours of the points.

        You should have set the points (using SetPoints()) before call this member.
    */
    void generatePoints(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateWireframe(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates the 3D graph for the dataset. Arguments provided are for the colours of the edges of the graph.
        \warning Modifies the PolyData pipeline so that GetOutput() will return the wireframe.
    */
    void generateWireframe(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateSurface(float red = defaultColour, float green = defaultColour, float blue = defaultColour)
        \brief Generates a surface display of the model, i.e. points at each point for the dataset. Arguments provided are for the colours of the points.

        You should have set the points (using SetPoints()) before call this member.
    */
    void generateSurface(float red = defaultColour, float green = defaultColour, float blue = defaultColour);
    /*!
        \fn milxQtModel::generateIsoSurface(vtkSmartPointer<vtkImageData> img, int contourNumber = -1, double value = 0.0)
        \brief Generates a iso surface or contour from image data, i.e. the result of Marching Cubes algorithm at value.
    */
    void generateIsoSurface(vtkSmartPointer<vtkImageData> img, int contourNumber = -1, double value = 0.0);
    /*!
        \fn milxQtModel::generatePolyDataFromImage(vtkSmartPointer<vtkImageData> img)
        \brief Generates polygonal data from image data, i.e. the result is a surface.
    */
    void generatePolyDataFromImage(vtkSmartPointer<vtkImageData> img);
    /*!
        \fn milxQtModel::generateKMeansClustering(int numberOfClusters = 0)
        \brief Generates K-means clustering of the points in the model with the number of clusters provided.
    */
    void generateKMeansClustering(int numberOfClusters = 0);
    /*!
        \fn milxQtModel::generateQuantisedPoints(float quantiseFactor = 0.0)
        \brief Generates points that are snapped to a grid given by factor.
    */
    void generateQuantisedPoints(float quantiseFactor = 0.0);
    /*!
        \fn milxQtModel::generateCappedBoundaries()
        \brief Generates capped boundaries for meshes with boundaries (such as large holes or voids).

      This is useful to cap holes in an open mesh so it can be voxelised among other things.
    */
    void generateCappedBoundaries();
    /*!
        \fn milxQtModel::generateRegionLabels()
        \brief Generates scalars for unconnected regions of the mesh
    */
    void generateRegionLabels();
    /*!
        \fn milxQtModel::generateElevation()
        \brief Generates scalars based on the elevation from the bounding box extremes of the mesh.
    */
    void generateElevation();
    /*!
        \fn milxQtModel::generateReebGraphs()
        \brief Generates Reeb graph from the scalars of the mesh. Uses elevation scalars if no scalars found.
    */
    //~ void generateReebGraphs();

    /*!
        \fn milxQtModel::updateCoords(vtkObject *obj)
        \brief Picks the coordinates and pixel value from the current mouse position in the window.
    */
    virtual void updateCoords(vtkObject *obj);

signals:
    /*!
        \fn milxQtModel::surfaceToImage(vtkSmartPointer<vtkPolyData> )
        \brief Emit signal to compute surface to image process.
    */
    void surfaceToImage(vtkSmartPointer<vtkPolyData> );
    /*!
        \fn milxQtModel::resultAvailable(milxQtRenderWindow*)
        \brief Send signal that Resultant render window is available for showing.
    */
    void resultAvailable(milxQtRenderWindow*);
    /*!
        \fn milxQtModel::resultAvailable(milxQtModel*)
        \brief Send signal that Resultant model is available for showing.
    */
    void resultAvailable(milxQtModel*);

protected:
    //Flags
    bool appended; //!< Appended data present?
    bool modelled; //!< Model has been generated?
    bool labelled; //!< Labels generated?
    bool scalarsSet; //!< Scalars set to points of mesh?
    bool transformed; //!< Transformed?
    bool computedCentroid; //!< Computed centroid already?
    bool computedCovariance; //!< Computed covariance matrix already?
    bool largeMode; //!< Large Data set mode enabled?
    bool outlineBefore; //!< scale displayed?
    bool cubeAxesBefore; //!< scale displayed?
    bool undidProcess; //!< Undid a process?

    float colourRed; //!< Saved colour, used for refresh etc.
    float colourGreen; //!< Saved colour, used for refresh etc.
    float colourBlue; //!< Saved colour, used for refresh etc.

    milx::Model model; //!< Actual model and its operations (from SMILI)

    //Data attributes
    coordinate modelCentroid; //!< Centroid coordinate of the model
    vnl_matrix<double> modelCovarianceMatrix; //!< Covariance matrix of the model

    //Display Attributes
    vtkSmartPointer<vtkPolyDataMapper> modelMapper; //!< Model mapper
    vtkSmartPointer<vtkLODActor> modelActor; //!< Model actor
    vtkSmartPointer<vtkLODActor> delaunayActor; //!< Delaunay triangulator actor
    vtkSmartPointer<vtkActor2D> modelLabelsActor; //!< Label Actor
    vtkSmartPointer<vtkActor> outlineActor; //!< Outline box actor
    vtkSmartPointer<vtkOutlineFilter> outlineMesh; //!< outline box
    vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor; //!< outline box

    //------------------
    //Generate Menu
    QMenu* generateMenu; //!< Menu for generating viewing options
    QAction* genModelAct; //!< Action for generating model
    QAction* genVerticesAct; //!< Action for generating vertices
    QAction* genNormalsAct; //!< Action for generating normals
    QAction* genVectorsAct; //!< Action for generating vector fields
    QAction* genTensorsAct; //!< Action for generating vector fields
    QAction* genHedgehogAct; //!< Action for generating Hedgehog fields
    QAction* genDelaunayAct; //!< Action for generating Delaunay graphs
    QAction* genDelaunayTri2DAct; //!< Action for generating Delaunay triangulation
    QAction* genDelaunayTriAct; //!< Action for generating Delaunay triangulation
    QAction* genLabelsAct; //!< Action for generating labels
    QAction* genPointIDsAct; //!< Action for generating labels
    QAction* genPointModelAct; //!< Action for generating labels
    QAction* genSamplesAct; //!< Action for generating sampled points on mesh
    QAction* genKmeansAct; //!< Action for calculating k-means clustering of points in a model
    QAction* genQuantiseAct; //!< Action for quantising points in a model
    QAction* genCapBoundaries; //!< Action for capping the boundaries of a open mesh
    QAction* genRegionLabels; //!< Action for labelling unconnected regions of mesh
    QAction* genElevation; //!< Action for elevation scalars of mesh
    QAction* genReebGraph; //!< Action for Reeb graphs of mesh
    //Operations Menu
    QMenu* operateMenu; //!< Operate Menu
    QAction* cleanAct; //!< Action for cleaning model
    QAction* triAct; //!< Action for triangulate model
    QAction* decimateAct; //!< Action for decimating model
    QAction* quadricDecimateAct; //!< Action for Quadric decimating model
    QAction* clusterDecimateAct; //!< Action for Cluster decimating model
    QAction* smoothAct; //!< Action for smooth model
    QAction* smoothSincAct; //!< Action for smooth model
    QAction* curvatureAct; //!< Action for calculating mean curvation of a model
    //Transform Menu
    QMenu* transformMenu; //!< Transform Menu
    QAction* transformAct; //!< Action for transforming model
    QAction* matchAct; //!< Action for matching scale and centroids of two models
    QAction* registerAct; //!< Action for registering two models
    QAction* registerLandmarkAct; //!< Action for registering two models via landmark transform
    QAction* voxeliseAct; //!< Action for voxelising surface
    //Scalars Menu
    QMenu* scalarMenu; //!< Context Menu
    QActionGroup* scalarsGroup; //!< Grouping for scalars
    QAction* renameScalarsAct; //!< Action for renaming scalars array
    QAction* thresholdAct; //!< Action for thresholding model
    QAction* thresholdScalarsAct; //!< Action for thresholding scalar values on model
    QAction* thresholdScalarsBinaryAct; //!< Action for binary thresholding scalar values on model
    QAction* maskScalarsAct; //!< Action for masking scalar values on model
    QAction* gradientAct; //!< Action for calculating gradient of the scalars of a model
    QAction* removeScalarsAct; //!< Action for removing scalars from model
    QAction* scalarsAct; //!< Action for loading scalars from another model
    QAction* outScalarsAct; //!< Action for output of scalars
    //Show Menu
    QMenu* showMenu; //!< Camera Menu
    QAction* normalsAct; //!< Action for showing the normals of a model
    QAction* centroidAct; //!< show centroid action
    QAction* outlineAct; //!< show outline action
    QAction* cubeAxesAct; //!< show cube axes action
    QAction* overlaysAct; //!< removes actors/overlays action
    //Contour menu
    QAction* selectPointsAct; //!< select points within contour
    QAction* weightAct; //!< produce Gaussian weights from contour

    //Context Menu
    //menu defined in milxQtWindow
    //------------------
    QAction* doAct; //!< Action for undoing and redoing actions
    //------------------
    QAction* infoAct; //!< Action for showing the model info
    QAction* textureAct; //!< Action for extracting texture
    QAction* flipAct; //!< Action for flipping surface
    QAction* rotateAct; //!< Action for rotating surface
    //------------------
    QAction* colourAct; //!< Action for changing colours of a model
    QAction* interpAct; //!< Action for changing the interpolation of a model
    QAction* specularAct; //!< Action for changing the specular of a model to flat

    //View Menu
    //------------------
    QActionGroup* displayGroup; //!< Grouping for check boxes
    QAction* pointsAct; //!< Show points
    QAction* wireframeAct; //!< Show wireframe
    QAction* surfaceAct; //!< Show surface
    QActionGroup* arrayGroup; //!< Grouping for check boxes
    QMenu* arrayMenu; //!< arrays Menu
    //------------------

    //Helper Members
    /*!
        \fn milxQtModel::resetFlags()
        \brief Reset the flags for pointers and computations. Useful when the data changes and all previous results need updating.
    */
    void resetFlags();
    /*!
        \fn milxQtModel::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
        \fn milxQtModel::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
        \fn milxQtModel::setupTooltips()
        \brief Assign tooltips to each action for users benefit. The tooltips explains the function of each action
    */
    void setupTooltips();
    /*!
    	\fn milxQtModel::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);
    /*!
    	\fn milxQtModel::basicContextMenu()
    	\brief Return the basic context menu with the milxQtModel class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* basicContextMenu();
    /*!
    	\fn milxQtModel::generationMenu()
    	\brief Return the generate menu with the milxQtModel class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* generationMenu();
    /*!
    	\fn milxQtModel::operationsMenu()
    	\brief Return the operations menu with the milxQtModel class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* operationsMenu();
    /*!
    	\fn milxQtModel::transformsMenu()
    	\brief Return the transform menu with the milxQtModel class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* transformsMenu();
    /*!
    	\fn milxQtModel::scalarsMenu()
    	\brief Return the scalars menu with the milxQtModel class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* scalarsMenu();
    /*!
    	\fn milxQtModel::arraysMenu()
    	\brief Return the arrays menu with the milxQtModel class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* arraysMenu();

private:
    QPointer<milxQtModel> normalsModel;
    QPointer<milxQtModel> centroidModel;
    QPointer<milxQtModel> outlineModel;
};

#endif // MILXQTMODEL_H
