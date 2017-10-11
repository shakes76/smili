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
#ifndef MILXQTRENDERWINDOW_H
#define MILXQTRENDERWINDOW_H

#include "milxQtWindow.h"

#include <QWidgetAction>
#include <QLabel>
#include <QHBoxLayout>
#include <QStatusBar>

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkActorCollection.h>
#include <vtkAxesActor.h>
#include <vtkLookupTable.h>
#include <vtkImageData.h>
#include <vtkImageActor.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkTable.h>
#include <vtkScalarBarActor.h>
//Widgets
#include <vtkScalarBarWidget.h>
#include <vtkContourWidget.h>
#include <vtkLineWidget2.h>
#include <vtkDistanceWidget.h>
#include <vtkBiDimensionalWidget.h>
#include <vtkAngleWidget.h>
#include <vtkPlaneWidget.h>
#include <vtkBoxWidget2.h>
#include <vtkSphereRepresentation.h>
#include <vtkSphereWidget2.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkTextWidget.h>

//Enum for views
enum ViewType { AXIAL = 0, CORONAL = 1, SAGITTAL = 2 };
enum ViewerType { SINGLE = 0, SCANNER = 1};
enum OrientationType { RADIOLOGICAL = 0, NEUROLOGICAL = 1};

//Struct for image actors, since actors per view is only possible need to link to parent (otherwise you get white actor)
struct ImageActorItem
{
    vtkSmartPointer<vtkImageActor> parentActor; //The parent actor it's connected to
    vtkSmartPointer<vtkImageActor> imageActor; //image actor to display in 3D view
};

class LabelledAction : public QWidgetAction {
public:
    LabelledAction (const QString& title, QPixmap pix, QWidget *theParent = 0) :
      QWidgetAction(theParent) {
        pWidget = new QWidget(theParent);
        pLayout = new QHBoxLayout();
        pLabel = new QLabel(title);
        pPixmapLabel = new QLabel(title);
        pPixmapLabel->setPixmap(pix);
        pPixmapLabel->setScaledContents(true);
        pLayout->addWidget (pPixmapLabel);
        pLayout->addWidget (pLabel);
        pWidget->setLayout (pLayout);
        setDefaultWidget(pWidget);
    }

    QLabel* label() {
        return pLabel;
    }
    QLabel* pixmapLabel() {
        return pPixmapLabel;
    }

private:
    QLabel *pLabel;
    QLabel *pPixmapLabel;
    QWidget* pWidget;
    QHBoxLayout* pLayout;
};

/*!
    \class milxQtRenderWindow
    \brief This class represents the MILX Qt Render Window Display object using QVTK.
    \author Shekhar S. Chandra, 2013

    This is a base class for VTK windows. It encapsulates all the basic properties of VTK displays. You should
    use the milxQtModel milxQtImage etc. classes instead of this class.

    Usage Examples:
    Displaying a sphere
    \code
    ///Define Test Polygon
    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetThetaResolution(32);
        sphere->SetPhiResolution(32);
        sphere->SetRadius(10.0);
        sphere->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphere->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

    renderWin->setName("Sphere");
        renderWin->AddActor(actor);
        renderWin->generateRender();
        renderWin->show();
    \endcode
    Notice how in the above example requires a lot of VTK code. This is further reduced by using milxQtModel, when displaying models.

    The rendering is encapsulated within a milxQtWindow widget.
*/
class MILXQT_EXPORT milxQtRenderWindow : public milxQtWindow
{
    Q_OBJECT

public:
    /*!
        \fn milxQtRenderWindow::milxQtRenderWindow(QWidget *parent = 0, bool contextSystem = true)
        \brief The standard constructor
    */
    milxQtRenderWindow(QWidget *theParent = 0, bool contextSystem = true);
    /*!
        \fn milxQtRenderWindow::~milxQtRenderWindow()
        \brief The standard destructor
    */
    virtual ~milxQtRenderWindow();

    /*!
        \fn milxQtRenderWindow::contextMenuSystem(bool context)
        \brief Disable/enable context menu system, should be on by default otherwise VTK will steal all right click events.

        Also sets up render window object and widgets.

        WARNING: Use this caarefully and straight after creation.
    */
    void contextMenuSystem(bool context);

    /*!
        \fn milxQtRenderWindow::AddActor(vtkSmartPointer<vtkProp> actor)
        \brief Add a VTK actor to this window.
    */
    inline void AddActor(vtkSmartPointer<vtkProp> actor)
    {
        renderer->AddActor(actor);
    }
    /*!
        \fn milxQtRenderWindow::AddVolume(vtkSmartPointer<vtkVolume> actor)
        \brief Add a VTK volume to this window.
    */
    inline void AddVolume(vtkSmartPointer<vtkVolume> actor)
    {
        renderer->AddVolume(actor);
    }
    /*!
        \fn milxQtRenderWindow::AddActor2D(vtkSmartPointer<vtkActor2D> actor)
        \brief Add a VTK actor to this window.
    */
    inline void AddActor2D(vtkSmartPointer<vtkProp> actor)
    {
        renderer->AddActor2D(actor);
    }

    /*!
        \fn milxQtRenderWindow::SetActor(vtkSmartPointer<vtkActor> actor)
        \brief Add a VTK actor to this window. Alias for AddActor.
    */
    inline void SetActor(vtkSmartPointer<vtkProp> actor)
    {
        renderer->AddActor(actor);
    }
    /*!
        \fn milxQtRenderWindow::SetActor2D(vtkSmartPointer<vtkActor2D> actor)
        \brief Add a VTK actor to this window. Alias for AddActor2D.
    */
    inline void SetActor2D(vtkSmartPointer<vtkProp> actor)
    {
        renderer->AddActor2D(actor);
    }
    /*!
        \fn milxQtRenderWindow::SetBackground(double red, double green, double blue)
        \brief Set the background of the display window.
    */
    inline void SetBackground(double red, double green, double blue)
    {
        renderer->GradientBackgroundOff();
        renderer->SetBackground(red, green, blue);
    }
    /*!
        \fn milxQtRenderWindow::SetSize(int height, int width)
        \brief Set the size of the display window.
    */
    inline void SetSize(int height, int width)
    {
        QVTKWidget::GetRenderWindow()->SetSize(width, height);
        int *winSize = QVTKWidget::GetRenderWindow()->GetSize();
        QVTKWidget::resize(winSize[0],winSize[1]);
    }
    /*!
        \fn milxQtRenderWindow::SetLookupTable(vtkSmartPointer<vtkLookupTable> lut)
        \brief Set the lookup table of the colours in display window.
    */
    inline void SetLookupTable(vtkSmartPointer<vtkLookupTable> lut)
    {
        lookupTable = lut;
    }

    /*!
        \fn milxQtRenderWindow::RemoveAllActors()
        \brief Remove all VTK actors from this window.
    */
    inline void RemoveAllActors()
    {
        renderer->RemoveAllViewProps();
    }
    /*!
        \fn milxQtRenderWindow::RemoveActor(vtkSmartPointer<vtkProp> actor)
        \brief Remove the VTK actor from this window.
    */
    inline void RemoveActor(vtkSmartPointer<vtkProp> actor)
    {
        renderer->RemoveActor(actor);
    }
    /*!
        \fn milxQtRenderWindow::RemoveActor2D(vtkSmartPointer<vtkProp> actor)
        \brief Remove the VTK 2D actor from this window.
    */
    inline void RemoveActor2D(vtkSmartPointer<vtkProp> actor)
    {
        renderer->RemoveActor(actor);
    }

    /**
      \brief Get all the actors in the render window
    */
    inline vtkActorCollection* GetActors()
    {
        return renderer->GetActors();
    }
    /**
      \brief Get all the actors in the render window
    */
    inline vtkVolumeCollection* GetVolumes()
    {
        return renderer->GetVolumes();
    }
    /**
      \brief Get the image actors, Implement in derived class that uses images
    */
    inline virtual vtkImageActor* GetImageActor()
    {
        if(!imageActors.empty())
            return imageActors.last().imageActor;       //Implement in derived class
        else
            return NULL;
    }
    /*!
        \fn milxQtRenderWindow::GetLookupTable()
        \brief Get the lookup table of the colours in display window.
    */
    inline vtkLookupTable* GetLookupTable()
    {
        return lookupTable;
    }
    /*!
        \fn milxQtRenderWindow::SetScalarBar(vtkScalarBarActor* bar)
        \brief Get the scalar bar of the display window.
    */
    inline void SetScalarBar(vtkScalarBarActor* bar)
    {
        scalarBar->SetScalarBarActor(bar);
    }
    /*!
        \fn milxQtRenderWindow::GetScalarBar()
        \brief Get the scalar bar of the display window.
    */
    inline vtkScalarBarActor* GetScalarBar()
    {
      return scale;
    }
    /*!
        \fn milxQtRenderWindow::GetBackground(double &red, double &green, double &blue)
        \brief Get the background of the display window.
    */
    inline void GetBackground(double &red, double &green, double &blue)
    {
        renderer->GetBackground(red, green, blue);
    }
    /**
      \brief Get the data, Implement in derived class that uses data derived from vtkDataSet, like vtkPolyData or vtkImageData.

      Useful for getting scalar range etc.
    */
    inline virtual vtkDataSet* GetDataSet()
    {
        return NULL;       //Implement in derived class that uses images
    }

    /*!
        \fn milxQtRenderWindow::SetRenderer(vtkRenderer* rnder)
        \brief Assigns a VTK Renderer object.
    */
    inline void SetRenderer(vtkRenderer* rnder)
    {
        renderer = rnder;
    }
    /*!
        \fn milxQtRenderWindow::GetRenderer()
        \brief Returns the VTK Renderer object.
    */
    inline vtkRenderer* GetRenderer()
    {
        return renderer;
    }
    /*!
      \brief Get the interactor associated with the view rendering
    */
    inline virtual vtkRenderWindowInteractor* GetVTKInteractor()
    {
        return QVTKWidget::GetRenderWindow()->GetInteractor();
    }

    /*!
        \fn milxQtRenderWindow::OffScreenRenderingOn()
        \brief Enable off-screen rendering, which is faster for large datasets in some instances.
    */
    inline void OffScreenRenderingOn()
    {
        QVTKWidget::GetRenderWindow()->OffScreenRenderingOn();
    }
    /*!
        \fn milxQtRenderWindow::OffScreenRenderingOff()
        \brief Disable off-screen rendering, which is faster for large datasets in some instances.
    */
    inline void OffScreenRenderingOff()
    {
        QVTKWidget::GetRenderWindow()->OffScreenRenderingOff();
    }

    /*!
        \fn milxQtRenderWindow::Render()
        \brief Force Render or Update of the display.
    */
    inline void Render()
    {
        QVTKWidget::GetRenderWindow()->Render();
    }

    /*!
        \fn milxQtRenderWindow::enableUpdates(QStatusBar *bar)
        \brief Enables the update of coordinates to the status bar provided directly. Update is done on mouse movement.
    */
    void enableUpdates(QStatusBar *bar);

public slots:
    /*!
        \fn milxQtRenderWindow::addModelActor(vtkSmartPointer<vtkActor> mdlActor)
        \brief Directly add model actor to generic view.
    */
    void addModelActor(vtkSmartPointer<vtkActor> mdlActor);
    /*!
        \fn milxQtRenderWindow::removeModelActor(vtkSmartPointer<vtkActor> mdlActor)
        \brief Directly remove model actor from generic view.
    */
    void removeModelActor(vtkSmartPointer<vtkActor> mdlActor);
    /*!
        \fn milxQtRenderWindow::addActor(vtkSmartPointer<vtkActor> imgActor, vtkSmartPointer<vtkMatrix4x4> transformMatrix = NULL)
        \brief Directly add actor to generic view with transform matrix.
    */
    void addActor(vtkSmartPointer<vtkActor> imgActor, vtkMatrix4x4 *transformMatrix = NULL);
    /*!
    \fn milxQtRenderWindow::addImageActor(vtkSmartPointer<vtkImageActor> imgActor, vtkSmartPointer<vtkMatrix4x4> transformMatrix = NULL)
    \brief Directly add image actor to generic view.
    */
    void addImageActor(vtkSmartPointer<vtkImageActor> imgActor, vtkMatrix4x4 *transformMatrix = NULL);
    /*!
        \fn milxQtRenderWindow::removeImageActor(vtkSmartPointer<vtkImageActor> imgActor)
        \brief Directly remove image actor from generic view.
    */
    void removeImageActor(vtkSmartPointer<vtkImageActor> imgActor);
    /*!
        \fn milxQtRenderWindow::importFrom(milxQtRenderWindow *windowToImportFrom)
        \brief Import actors/scene directly from another renderwindow
    */
    void importFrom(milxQtRenderWindow *windowToImportFrom);
    /*!
        \fn milxQtRenderWindow::importViewFrom(milxQtRenderWindow *windowToImportFrom)
        \brief Same as importFrom()
    */
    inline void importViewFrom(milxQtRenderWindow *windowToImportFrom)
    {   importFrom(windowToImportFrom);   }
    /*!
      \fn milxQtRenderWindow::updateTextActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
      \brief Update any text actors in display to current slice.
    */
    virtual void updateTextActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command);
    /*!
      \fn milxQtRenderWindow::updateImageActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
      \brief Update any image actors in display to current slice when slice number changes.
    */
    virtual void updateImageActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command);
    virtual void updateImageActor(vtkSmartPointer<vtkImageActor> actor);

    void userEvent(QMouseEvent *event = NULL);
    void userEvent(QKeyEvent *event);
    void userEvent(QWheelEvent *event);

    /*!
        \fn milxQtRenderWindow::refresh()
        \brief Refresh the display of the model including widgets. Camera remains as is.
    */
    void refresh();
    /*!
        \fn milxQtRenderWindow::getTransformMatrix()
        \brief Internal transform matrix used to transform actors appropriately for direction and coordinate systems
    */
    inline vtkSmartPointer<vtkMatrix4x4> getTransformMatrix()
    {   return transformMatrix;   }
    /*!
        \fn reset()
        \brief Reset the rendering, camera and windowing.
    */
    void reset();
    /*!
        \fn milxQtRenderWindow::isLoaded()
        \brief Is the data successfully loaded? Use this to check if all is ok with data display.
    */
    inline bool isLoaded()
    {
        return loaded;
    }
    /*!
        \fn milxQtRenderWindow::linkProgressEventOf(vtkObject * obj)
        \brief Link the progress of filters etc to keep the UI responsive.
    */
    void linkProgressEventOf(vtkObject * obj);
    /*!
        \fn milxQtRenderWindow::updateCoords(vtkObject *obj)
        \brief Picks the coordinates and pixel value from the current mouse position in the window.
    */
    virtual void updateCoords(vtkObject *obj);
    /*!
        \fn milxQtRenderWindow::contour()
        \brief Draw contour interactively
    */
    virtual void contour();
    /*!
        \fn milxQtRenderWindow::contourAsPolyData()
        \brief Save contour as polydata/model. This saves all intermediate points, not just what the user places.
    */
    void contourAsPolyData();
    /*!
        \fn milxQtRenderWindow::contourAsNodePolyData()
        \brief Save contour as polydata/model. This only saves the points the user has placed.
    */
    void contourAsNodePolyData();
    /*!
        \fn milxQtRenderWindow::contourInitFromPolyData(QString filename = "")
        \brief Load contour from polydata/model
    */
    void contourInitFromPolyData(QString filename = "");
    /*!
        \fn milxQtRenderWindow::enableAxes(float xScale = 1.0, float yScale = 1.0, float zScale = 1.0)
        \brief Enable axes display with each arrow for dimension scaled as provided.
    */
    void enableAxes(float xScale = 1.0, float yScale = 1.0, float zScale = 1.0);
    /*!
        \fn milxQtRenderWindow::disableAxes()
        \brief Disables the axes display.
    */
    void disableAxes();
    /*!
        \fn milxQtRenderWindow::axesDisplay()
        \brief Toggles the axes display.
    */
    void axesDisplay();
    /*!
        \fn milxQtRenderWindow::disableOrient()
        \brief Disables orientation marker.
    */
    void disableOrient();
    /*!
        \fn milxQtRenderWindow::orientDisplay()
        \brief Toggles orientation marker.
    */
    void orientDisplay();
    /*!
        \fn milxQtRenderWindow::background(bool white = false)
        \brief Changes background to white.
    */
    virtual void background(bool white = false);
    /*!
        \fn milxQtRenderWindow::lighting()
        \brief Toggles two sided lighting in display.
    */
    void lighting();
    /*!
        \fn milxQtRenderWindow::textDisplay()
        \brief Toggles text in display.
    */
    void textDisplay();
    /*!
        \fn milxQtRenderWindow::crosshair()
        \brief Toggles corsshair.
    */
    void crosshair();

    //Viewing
    /*!
        \fn milxQtRenderWindow::viewToXYPlane()
        \brief Change view to xy-plane.
    */
    virtual void viewToXYPlane();
    inline void viewToAxial()
    {
        viewToXYPlane();
    }
    /*!
        \fn milxQtRenderWindow::viewToZXPlane()
        \brief Change view to zx-plane.
    */
    virtual void viewToZXPlane();
    inline void viewToCoronal()
    {
        viewToZXPlane();
    }
    /*!
        \fn milxQtRenderWindow::viewToZYPlane()
        \brief Change view to zy-plane.
    */
    virtual void viewToZYPlane();
    inline void viewToSagittal()
    {
        viewToZYPlane();
    }
    /*!
        \fn milxQtRenderWindow::getDefaultOrientation()
        \brief Get the orientation mode to one of the supported standards.

        0-Radiological: Feet first view
        1-Neurological: Head first view
    */
    inline int getDefaultOrientation()
    {   return orientationView;  }
    /*!
        \fn milxQtRenderWindow::setDefaultOrientation(int orientMode)
        \brief Change orientation mode to one of the supported standards. Default: Radiological.

        0-Radiological: Feet first view
        1-Neurological: Head first view
    */
    inline void setDefaultOrientation(int orientMode)
    {   orientationView = orientMode;  }
    /*!
        \fn milxQtRenderWindow::getView()
        \brief Return current view mode identified by number.
        0-axial, 1-coronal, 2-sagittal
    */
    inline int getView()
    {   return currentView;   }
    /*!
        \fn milxQtRenderWindow::setView(int viewMode)
        \brief Change view to view mode identified by number.
        0-axial, 1-coronal, 2-sagittal
    */
    void setView(int viewMode);
    /*!
        \fn milxQtRenderWindow::enableActionBasedOnView()
        \brief Enables the view actions corresponding to current view set.
    */
    void enableActionBasedOnView();
    /*!
		\fn milxQtRenderWindow::disableDefaultView()
		\brief Disables the default view whenever data is process or displayed. Enabled by default.

		Used in the derived classes in general.
	*/
    inline void disableDefaultView()
    {	useDefaultView = false;	}
    /*!
		\fn milxQtRenderWindow::disableDefaultView()
		\brief Enables the default view whenever data is process or displayed. Enabled by default.

		Used in the derived classes in general.
	*/
    inline void enableDefaultView()
    {	useDefaultView = true;	}
    /*!
        \fn milxQtRenderWindow::setDefaultView(int viewMode)
        \brief Change default view to view mode identified by number.
        0-axial, 1-coronal, 2-sagittal
    */
    inline void setDefaultView(int viewMode)
    {   defaultView = viewMode;   }
    /*!
        \fn milxQtRenderWindow::setView(QString filename = "")
        \brief Saves the camera details to internal variables that can be restore anytime with loadView().
    */
    void saveView(QString filename = "");
    void saveViewFile();
    /*!
        \fn milxQtRenderWindow::loadView(QString filename = "")
        \brief Saves the camera details to internal variables that can be restore anytime with loadView().
    */
    void loadView(QString filename = "");
    void loadViewFile();
    /*!
        \fn milxQtRenderWindow::enableScale(QString title = "", const bool quiet = false, double minRange = 0.0, double maxRange = 0.0, int noOfLabels = 3)
        \brief Enable scale bar display with the title provided.

        Quiet Boolean is to prevent possible popups to ask user parameters.
    */
    virtual void enableScale(QString title = "", const bool quiet = false, double minRange = 0.0, double maxRange = 0.0, int noOfLabels = 3) {}
    /*!
        \fn milxQtRenderWindow::disableScale()
        \brief Disables the scale bar display.
    */
    virtual void disableScale();
    /*!
        \fn milxQtRenderWindow::enableCrosshair()
        \brief Enables the mouse pointer as a crosshair instead. Scene must be rendered before calling.
    */
    virtual inline void enableCrosshair()
    {
        if(rendered)
        {
            renderWindow->SetCurrentCursor(VTK_CURSOR_CROSSHAIR);
            crosshairAct->setChecked(true);
        }
    }
    /*!
        \fn milxQtRenderWindow::disableCrosshair()
        \brief Restores the mouse pointer to default.
    */
    virtual inline void disableCrosshair()
    {
        renderWindow->SetCurrentCursor(0);
        crosshairAct->setChecked(false);
    }

    /*!
        \fn milxQtRenderWindow::scaleDisplay(const bool forceDisplay = false)
        \brief Toggles the scale bar display.

        forceDisplay Boolean is to overide possible previous settings and display bar.
    */
    virtual void scaleDisplay(const bool forceDisplay = false) {}
    /*!
        \fn milxQtRenderWindow::colourMapToRainbow(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to default VTK which is rainbow
    */
    virtual void colourMapToRainbow(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToVTK(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to inverse of default VTK
    */
    virtual void colourMapToVTK(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToGray(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Gray
    */
    virtual void colourMapToGray(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToSeismic(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Seismic
    */
    virtual void colourMapToSeismic(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToLogGray(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Logarithmic (base 10) Gray
    */
    virtual void colourMapToLogGray(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToJet(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Jet
    */
    virtual void colourMapToJet(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToNIH(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to NIH
    */
    virtual void colourMapToNIH(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToNIH_Fire(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to NIH Fire
    */
    virtual void colourMapToNIH_Fire(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToAAL(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to AAL
    */
    virtual void colourMapToAAL(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToFS(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to FS (FreeSurfer)
    */
    virtual void colourMapToFS(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToHOT(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to HOT
    */
    virtual void colourMapToHOT(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToCOOL(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to COOL
    */
    virtual void colourMapToCOOL(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToCOOLWARM(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to COOL
    */
    virtual void colourMapToCOOLWARM(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToKnee(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Knee
    */
    virtual void colourMapToKnee(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToBone(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Bone
    */
    virtual void colourMapToBone(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToSpectral(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to Spectral
    */
    virtual void colourMapToSpectral(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToGNUPlot(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to GNUPlot
    */
    virtual void colourMapToGNUPlot(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToCubeHelix(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to CubeHelix
    */
    virtual void colourMapToCubeHelix(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::colourMapToHSV(double minRange = 0.0, double maxRange = 0.0)
        \brief Change the colour map to HSV
    */
    virtual void colourMapToHSV(double minRange = 0.0, double maxRange = 0.0);
    /*!
        \fn milxQtRenderWindow::updateLookupTable()
        \brief Implement this into your derived class to ensure the new colourmap is passed on to your viewing data (image, model etc.)
    */
    virtual void updateLookupTable();

    /*!
        \fn milxQtRenderWindow::generateRender()
        \brief Generate the render so it is ready for display. Should be called before showing the window.

        Consider using refresh() or reset() to regenerate the view rather than this member for performance reasons.
        You should only need to call this member once.
    */
    void generateRender();

    //Custom
    /*!
        \fn milxQtRenderWindow::setCustomActionGroup(QActionGroup *actionGrp)
        \brief Set the action group for the menu with custom connections.
    */
    inline void setCustomActionGroup(QActionGroup *actionGrp)
    {
        windowActionGroup = actionGrp;
    }
    /*!
        \fn milxQtRenderWindow::createCustomConnections(QActionGroup *actGroup)
        \brief Create a series of connections on each action in group which will be placed into the context menu
    */
    void createCustomConnections(QList<QAction *> actionsInMenu);
    /*!
        \fn milxQtRenderWindow::createCustomConnections(QMenu *fromMenu)
        \brief Create a series of connections on each action in group which will be placed into the context menu
    */
    void createCustomConnections(QMenu *fromMenu);
    /*!
        \fn milxQtRenderWindow::customOperation()
        \brief Custom operation, data dependent for viewing in unified environment. By default, this allows importing of actors into the current views.

        Default Implementation involves handling generic actors and image actor only. Implement further in child class.
    */
    virtual void customOperation();

    /*!
        \fn milxQtRenderWindow::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);

    /*!
        \fn milxQtRenderWindow::openModelUsingQt(const QString filename, vtkSmartPointer<vtkPolyData> &data)
        \brief Opens a model file using Qt file objects, which can be a Wavefront Object file (*.obj) only. This member is intended to be used with the Qt resource system.

        Only supports triangulated meshes atm.
        This is useful for support for 'qrc:/' and ':/' paths, hence allowing meshes to be stored within the executable binary. Returns true if successful.
    */
    bool openModelUsingQt(const QString filename, vtkSmartPointer<vtkPolyData> &data);

signals:
    /*!
        \fn milxQtRenderWindow::imageAvailable(vtkImageData*, QString )
        \brief Send signal that an image is available for showing.
    */
    void imageAvailable(vtkImageData*, QString );
    /*!
        \fn milxQtRenderWindow::modelAvailable(vtkPolyData*, QString )
        \brief Send signal that an surface etc. is available for showing.
    */
    void modelAvailable(vtkPolyData*, QString );
    /*!
        \fn milxQtRenderWindow::tableToPlot(vtkSmartPointer<vtkTable>, QString)
        \brief Emit signal to show the table as a plot.
    */
    void tableToPlot(vtkSmartPointer<vtkTable>, QString);
    /*!
        \fn milxQtRenderWindow::modified(vtkSmartPointer<vtkImageActor> )
        \brief Emit signal to allow updating of image actors if present.
    */
    void modified(vtkSmartPointer<vtkImageActor> );

protected:
    //Flags
    bool loaded; //!< Loaded Image from file?
    bool rendered; //!< Scene as been setup for rendering? (by generateRender())
    bool axesBefore; //!< Axes displayed?
    bool scaleBefore; //!< scale displayed?
    bool customScalarBar; //! Using custom scalar bar?
    bool logScale; //! Using log scalar map?
    bool useDefaultView; //!< Use default view whenever possible?
    bool orientationAxes; //!< Display orientation (posterior etc.) axes?
    bool contextMenuEnabled; //!< Display orientation (posterior etc.) axes?

    //Variables
    int defaultView; //!< Default view for data (default is axial)
    int currentView; //!< Current view for data
    int orientationView; //!< view orientation standard

    //Window Menu
    QMenu* windowPropertiesMenu; //!< Context Menu
    QAction* backgroundAct; //!< Action for axes of the display
    QAction* axesAct; //!< Action for axes of the display
    QAction* lightingAct; //!< Action for two-sided lighting of the display
    QAction* lineAct; //!< Action for distance measuring display
    QAction* distanceAct; //!< Action for distance measuring display
    QAction* biDirectionAct; //!< Action for cross distance measuring display
    QAction* angleAct; //!< Action for angle measuring display
    QAction* planeAct; //!< Action for drawing planes
    QAction* boxAct; //!< Action for box drawing display
    QAction* sphereAct; //!< Action for sphere annotate display
    QAction* humanAct; //!< Show human view orientation glyph?
    QAction* textAct; //!< Action for angle measuring display
    QAction* crosshairAct; //!< Action for crosshair

    //Contouring Menu
    QMenu* contourMenu; //!< Contour Menu
    QAction* contourAct; //!< Action for contouring surface points using Dijkstras algorithm
    QAction* contourPolyDataAct; //!< Save contour surface points as polydata.
    QAction* contourNodePolyDataAct; //!< Save contour surface nodes as polydata.
    QAction* contourInitAct; //!< Load contour surface points from polydata.
//    QAction* contourImageAct; //!< Save contour surface points as binary image.

    //View Menu
    QMenu* viewMenu; //!< Context Menu
    QActionGroup* viewGroup; //!< Grouping for check boxes
    QAction* viewXY; //!< Change view to xy-plane (Axial)
    QAction* viewZY; //!< Change view to zy-plane (Coronal)
    QAction* viewZX; //!< Change view to zx-plane (Sagittal)
    QAction* saveViewAct; //!< Save camera view
    QAction* loadViewAct; //!< Load camera view
    QAction* saveViewFileAct; //!< Save camera view to file
    QAction* loadViewFileAct; //!< Load camera view to file
    QAction* scaleAct; //!< Action for the scale bar of the display
    QMenu* colourMapMenu; //!< Colour map menu
    QActionGroup* mapGroup; //!< Grouping for check boxes
    QAction* actionDefault; //!< Default colours
    QAction* actionJet; //!< Jet colours
    QAction* actionRainbow; //!< Default colours
    QAction* actionInvRainbow; //!< Inv Default colours
    LabelledAction* actionGray; //!< Gray colours
    LabelledAction* actionSeismic; //!< Gray2 colours
    QAction* actionLogGray; //!< Logarithmic Gray colours
    QAction* actionNIH; //!< NIH colours
    QAction* actionNIH_FIRE; //!< NIH_FIRE colours
    QAction* actionAAL; //!< AAL colours
    QAction* actionFS; //!< FS colours
    LabelledAction* actionHOT; //!< HOT colours
    QAction* actionCOOL; //!< COOL colours
    LabelledAction* actionCOOLWARM; //!< COOL colours
    QAction* actionKnee; //!< Knee colours
    LabelledAction* actionBone; //!< Spectral colours
    LabelledAction* actionSpectral; //!< Spectral colours
    LabelledAction* actionGNUPlot; //!< GNUPlot colours
    LabelledAction* actionCubeHelix; //!< Cube Helix colours
    LabelledAction* actionHSV; //!< Cube Helix colours

    //Actions
    QAction* resetAct; //!< Action for refreshing the display
    QAction* refreshAct; //!< Action for refreshing the display

    QPoint dragStartPosition; //!< Dragging helper variable

    vtkSmartPointer<vtkEventQtSlotConnect> Connector; //!< VTK Events to slots convertor
    vtkSmartPointer<vtkMatrix4x4> transformMatrix; //!< Transform matrix for actor(s)

    //VTK Rendering Internals
    vtkSmartPointer<vtkRenderer> renderer; //!< Renderer for the data
    vtkSmartPointer<vtkRenderWindow> renderWindow; //!< Render Window used if no other set, only for deletion, access through QVTKWidget members
    vtkSmartPointer<vtkLookupTable> lookupTable; //!< Lookup table for the shapes/images, base class is used to allow references to different look up table types
    vtkSmartPointer<vtkCamera> camera; //!< camera for the view

    vtkSmartPointer<vtkScalarBarActor> scale; //!< Scale for the display
    vtkSmartPointer<vtkScalarBarWidget> scalarBar; //!< Scalar Bar Widget for the display
    vtkSmartPointer<vtkAxesActor> axes; //!< Axes for the model
    vtkSmartPointer<vtkAxesActor> orientAxes; //!< Orientation (posterior etc.) axes for the orientation marker
    vtkSmartPointer<vtkPointPicker> dataPicker; //!< For determining coordinates and points from the window

    //Widgets
    vtkSmartPointer<vtkContourWidget> contourWidget; //!< contour interaction
    vtkSmartPointer<vtkLineWidget2> lineWidget; //!< Used for drawing lines
    vtkSmartPointer<vtkDistanceWidget> distanceWidget; //!< Used for measuring distances
    vtkSmartPointer<vtkBiDimensionalWidget> biDirectionWidget; //!< Used for measuring cross distances
    vtkSmartPointer<vtkAngleWidget> angleWidget; //!< Used for measuring angles
    vtkSmartPointer<vtkPlaneWidget> planeWidget; //!< Used for drawing planes
    vtkSmartPointer<vtkBoxWidget2> boxWidget; //!< Used for measuring angles
    vtkSmartPointer<vtkSphereRepresentation> sphereRep; //!< Sphere for widgets
    vtkSmartPointer<vtkSphereWidget2> sphereWidget; //!< Used for measuring angles
    vtkSmartPointer<vtkOrientationMarkerWidget> humanGlyph; //!< Glyph for showing equivalent view on human
    QList< vtkSmartPointer<vtkTextWidget> > textWidgets; //!< Used for displaying movable texts in render window
    QList< vtkSmartPointer<vtkTextActor> > textActors; //!< text actors for text widgets
    QList<ImageActorItem> imageActors; //!< Images actors being displayed in model view

    QStatusBar *updateBar; //!< Pointer to bar, not allocated or deleted. To be passed to only.

    QActionGroup* windowActionGroup; //!< used for the custom menu

    /*!
        \fn milxQtRenderWindow::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
        \fn milxQtRenderWindow::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
    	\fn milxQtRenderWindow::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);
    /*!
    	\fn milxQtRenderWindow::colourMapsMenu()
    	\brief Return the colourmaps menu with the milxQtRenderWindow class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* colourMapsMenu();
    /*!
        \fn milxQtRenderWindow::dragMoveEvent(QDragMoveEvent *event)
        \brief Part of the Drag and Drop feature members. Accepts drags.
    */
    void dragMoveEvent(QDragMoveEvent *event);
    /*!
        \fn milxQtRenderWindow::dragLeaveEvent(QDragLeaveEvent *event)
        \brief Part of the Drag and Drop feature members. Accepts drags.
    */
    void dragLeaveEvent(QDragLeaveEvent *event);
    /*!
        \fn milxQtRenderWindow::dragEnterEvent(QDragEnterEvent *event)
        \brief Part of the Drag and Drop feature members. Tells what drags to accept.
    */
    void dragEnterEvent(QDragEnterEvent *event);
    /*!
        \fn milxQtRenderWindow::mouseDoubleClickEvent(QMouseEvent *event)
        \brief Part of the Drag and Drop feature members. Triggers the dragging on double click.
    */
    void mouseDoubleClickEvent(QMouseEvent *event);
    /*!
        \fn milxQtRenderWindow::dropEvent(QDropEvent *event)
        \brief Part of the Drag and Drop feature members. Opens the dropped files.
    */
    void dropEvent(QDropEvent *event);

    /*!
        \fn milxQtRenderWindow::SetupWidgets(vtkRenderWindowInteractor *interactor)
        \brief Update the interactor so that the widgets are usable. Call this if you change the render window or interator manually.
    */
    virtual void SetupWidgets(vtkRenderWindowInteractor *interactor);

    /*!
        \fn milxQtRenderWindow::setupHumanGlyph(vtkSmartPointer<vtkMatrix4x4> mat = NULL)
        \brief Setup the human orientation glyph even if some transform if present between glyph and view.
    */
    void setupHumanGlyph(vtkSmartPointer<vtkMatrix4x4> mat = NULL);

private:

};

#endif // MILXQTRENDERWINDOW_H
