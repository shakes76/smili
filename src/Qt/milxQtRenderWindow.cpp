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
#include "milxQtRenderWindow.h"

#include <vtkRendererCollection.h>
#include <vtkPropAssembly.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkDistanceRepresentation.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkBoxRepresentation.h>
#include <vtkWidgetEvent.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCellArray.h> //needed for loading orientation glyph model
#include <vtkFloatArray.h> //needed for loading orientation glyph model
#include <vtkPointData.h>//needed for loading orientation glyph model
#if VTK_MAJOR_VERSION > 5
    #include <vtkImageMapper3D.h>
#endif

#include "milxColourMap.h"
#include "milxFile.h"

milxQtRenderWindow::milxQtRenderWindow(QWidget *theParent, bool contextSystem) : milxQtWindow(theParent)
{
    ///Init
    loaded = false;
    rendered = false;
    axesBefore = false;
    scaleBefore = false;
    customScalarBar = false;
    logScale = false;
    useDefaultView = true;
    orientationAxes = true;
    contextMenuEnabled = contextSystem;
    defaultView = AXIAL; //axial
    currentView = AXIAL; //axial
    orientationView = RADIOLOGICAL;

    ///Allocate critical aspects
    renderer = vtkSmartPointer<vtkRenderer>::New();
    lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    camera = vtkSmartPointer<vtkCamera>::New();

    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();

    contextMenuSystem(contextMenuEnabled); //also sets up render window object and widgets

    //orientation glyph (showing equivalent view on human) setup
    setupHumanGlyph();

    setAcceptDrops(true);

    createActions();

    Connector = vtkSmartPointer<vtkEventQtSlotConnect>::New();

    createConnections();

    reset();
}

milxQtRenderWindow::~milxQtRenderWindow()
{
    textWidgets.clear();
    textActors.clear();
}

void milxQtRenderWindow::contextMenuSystem(bool context)
{
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    QVTKWidget::SetRenderWindow(renderWindow);
    QVTKWidget::GetRenderWindow()->StereoRenderOff(); ///Force Stereo Render off
  #if(VTK_MAJOR_VERSION < 6)
    QVTKWidget::GetRenderWindow()->ReportGraphicErrorsOn(); ///Force Error reporting
  #endif

    //widget setups
    SetupWidgets(QVTKWidget::GetInteractor());

    ///Performance
    QVTKWidget::GetRenderWindow()->GetInteractor()->SetDesiredUpdateRate(20); //for LOD actors
    QVTKWidget::GetRenderWindow()->GetInteractor()->SetStillUpdateRate(0.01); //for LOD actors

    if(context)
    {
        ///Unbind the right mouse button events as Qt uses context menu.
        QVTKWidget::GetRenderWindow()->GetInteractor()->RemoveObservers(vtkCommand::RightButtonPressEvent);
        QVTKWidget::GetRenderWindow()->GetInteractor()->RemoveObservers(vtkCommand::RightButtonReleaseEvent);
    }
}

void milxQtRenderWindow::linkProgressEventOf(vtkObject * obj)
{
    Connector->Connect(obj,
                       vtkCommand::ProgressEvent,
                       this,
                       SLOT( updateQtEvents() ),
                       NULL, 1.0); //High Priority
}

void milxQtRenderWindow::SetupWidgets(vtkRenderWindowInteractor *interactor)
{
    //distance widget setup
    lineWidget = vtkSmartPointer<vtkLineWidget2>::New();
        lineWidget->SetInteractor(interactor);
        lineWidget->CreateDefaultRepresentation();

    //distance widget setup
    distanceWidget = vtkSmartPointer<vtkDistanceWidget>::New();
        distanceWidget->SetInteractor(interactor);
        distanceWidget->CreateDefaultRepresentation();
        vtkDistanceRepresentation::SafeDownCast(distanceWidget->GetRepresentation())->SetLabelFormat("%-#6.3g mm");

    //cross distance widget setup
    biDirectionWidget = vtkSmartPointer<vtkBiDimensionalWidget>::New();
        biDirectionWidget->SetInteractor(interactor);
        biDirectionWidget->CreateDefaultRepresentation();

    //angle widget setup
    angleWidget = vtkSmartPointer<vtkAngleWidget>::New();
        angleWidget->SetInteractor(interactor);
        angleWidget->CreateDefaultRepresentation();

    //plane widget setup
    planeWidget = vtkSmartPointer<vtkPlaneWidget>::New();
          planeWidget->SetInteractor(interactor);

    //box widget
    boxWidget = vtkSmartPointer<vtkBoxWidget2>::New();
      boxWidget->SetInteractor(interactor);
      //boxWidget->CreateDefaultRepresentation();

    vtkSmartPointer<vtkBoxRepresentation> boxRepresentation = vtkSmartPointer<vtkBoxRepresentation>::New();
    boxWidget->SetRepresentation(boxRepresentation);

    //sphere widget setup
    sphereRep = vtkSmartPointer<vtkSphereRepresentation>::New();
        sphereRep->SetRepresentationToWireframe();
        sphereRep->SetThetaResolution(36);
        sphereRep->SetPhiResolution(36);
        sphereRep->SetRadius(10);
        sphereRep->HandleTextOn();
        sphereRep->RadialLineOn();
        sphereRep->HandleVisibilityOn();
    sphereWidget = vtkSmartPointer<vtkSphereWidget2>::New();
        sphereWidget->SetInteractor(interactor);
//        sphereWidget->CreateDefaultRepresentation();
        sphereWidget->SetRepresentation(sphereRep);
        sphereWidget->ScalingEnabledOn();
    // Change bindings.
    vtkWidgetEventTranslator *eventTranslator = sphereWidget->GetEventTranslator();
    //unbind
//    eventTranslator->RemoveTranslation( vtkCommand::LeftButtonPressEvent );
//    eventTranslator->RemoveTranslation( vtkCommand::LeftButtonReleaseEvent );
    eventTranslator->RemoveTranslation( vtkCommand::MiddleButtonPressEvent );
    eventTranslator->RemoveTranslation( vtkCommand::MiddleButtonReleaseEvent );
    eventTranslator->RemoveTranslation( vtkCommand::RightButtonPressEvent );
    eventTranslator->RemoveTranslation( vtkCommand::RightButtonReleaseEvent );
    //bind
//    eventTranslator->SetTranslation(vtkCommand::LeftButtonPressEvent, vtkWidgetEvent::Scale);
//    eventTranslator->SetTranslation(vtkCommand::LeftButtonReleaseEvent, vtkWidgetEvent::EndScale);
//    eventTranslator->SetTranslation(vtkCommand::MiddleButtonPressEvent, vtkWidgetEvent::Translate);
//    eventTranslator->SetTranslation(vtkCommand::MiddleButtonReleaseEvent, vtkWidgetEvent::EndTranslate);
    eventTranslator->SetTranslation(vtkCommand::MiddleButtonPressEvent, vtkWidgetEvent::Scale);
    eventTranslator->SetTranslation(vtkCommand::MiddleButtonReleaseEvent, vtkWidgetEvent::EndScale);
//    eventTranslator->SetTranslation(vtkCommand::MouseMoveEvent, vtkWidgetEvent::Scale);
}

void milxQtRenderWindow::contour()
{
    if(!contourWidget)
        contourWidget = vtkSmartPointer<vtkContourWidget>::New();
        contourWidget->SetInteractor(QVTKWidget::GetInteractor());

    if(contourAct->isChecked())
    {
        contourWidget->EnabledOn();
        printInfo("Enabled Contouring");
    }
    else
    {
        contourWidget->EnabledOff();
        printInfo("Disabled Contouring");
    }

//    vtkSmartPointer<vtkOrientedGlyphContourRepresentation> rep = vtkOrientedGlyphContourRepresentation::SafeDownCast( contourWidget->GetRepresentation() );
//    rep->GetLinesProperty()->SetColor(1, 0.2, 0);
//    rep->GetLinesProperty()->SetLineWidth(3.0);
//
//    vtkSmartPointer<vtkImageActorPointPlacer> pointPlacer = vtkSmartPointer<vtkImageActorPointPlacer>::New();
//        pointPlacer->AddProp(modelActor);
//        pointPlacer->GetPolys()->AddItem(currentMesh);
//    rep->SetPointPlacer(pointPlacer);
//
//    vtkSmartPointer<vtkDijkstraImageContourLineInterpolator> interpolator = vtkSmartPointer<vtkDijkstraImageContourLineInterpolator>::New();
//        interpolator->GetPolys()->AddItem(currentMesh);
//    rep->SetLineInterpolator(interpolator);

    refresh();
}

void milxQtRenderWindow::contourAsPolyData()
{
    if(!contourAct->isChecked())
        return;

    vtkSmartPointer<vtkContourRepresentation> rep = contourWidget->GetContourRepresentation();

    //rep->GetContourRepresentationAsPolyData() gives intermediate points also
    //GetNodePolyData() gives only the nodes placed by user
    emit modelAvailable(rep->GetContourRepresentationAsPolyData(), "Contour Model");
}

void milxQtRenderWindow::contourAsNodePolyData()
{
  if(!contourAct->isChecked())
    return;

  vtkSmartPointer<vtkContourRepresentation> rep = contourWidget->GetContourRepresentation();

  vtkSmartPointer<vtkPolyData> nodes = vtkSmartPointer<vtkPolyData>::New();
  rep->GetNodePolyData(nodes);

  //rep->GetContourRepresentationAsPolyData() gives intermediate points also
  //GetNodePolyData() gives only the nodes placed by user
  emit modelAvailable(nodes, "Contour");
}

void milxQtRenderWindow::contourInitFromPolyData(QString filename)
{
  if(!contourAct->isChecked())
    return;

  QSettings settings("Shekhar Chandra", "milxQt");
  //Add supported file entry
  QString exts = openModelExts.c_str();
  QString path = settings.value("recentPath").toString();

  if(filename.isEmpty())
  {
      QFileDialog *fileOpener = new QFileDialog(this);
      filename = fileOpener->getOpenFileName(this,
                              tr("Select Model(s) to Open"),
                              path,
                              tr(exts.toStdString().c_str()) ); //!< \todo Check and validate extensions support at Open in Main class
  }

  printDebug("Opening contour model");
  vtkSmartPointer<vtkPolyData> mesh;
  milx::File::OpenModel(filename.toStdString(), mesh);

  printDebug("Initialising Contouring");
  contourWidget->On();
  contourWidget->Initialize(mesh);
  contourWidget->Render();
}

void milxQtRenderWindow::enableAxes(float xScale, float yScale, float zScale)
{
    if(!axes)
        axes = vtkSmartPointer<vtkAxesActor>::New();

    axes->SetTotalLength(xScale, yScale, zScale);
    axes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToViewport();
    axes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToViewport();
    axes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToViewport();

    AddActor(axes);
    axesBefore = true;
}

void milxQtRenderWindow::disableAxes()
{
    if(!axesBefore)
        return; //nothing to do then

    RemoveActor(axes);
    axesBefore = false;
}

void milxQtRenderWindow::background(bool white)
{
    if(backgroundAct->isChecked() || white)
    {
        renderer->GradientBackgroundOff();
        renderer->SetBackground(1 , 1, 1);
        backgroundAct->setChecked(true);
        if(orientAxes)
        {
            orientAxes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(0, 0, 0);
            orientAxes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(0, 0, 0);
            orientAxes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(0, 0, 0);
        }
    }
    else
    {
        renderer->GradientBackgroundOn();
        renderer->SetBackground(1 , 1, 1);
        //renderer->SetBackground2(0.6, 0.5, 0.4);
        renderer->SetBackground2(0.4, 0.5, 0.6);
        if(orientAxes)
        {
            orientAxes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
            orientAxes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
            orientAxes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
        }
    }

    Render();
}

void milxQtRenderWindow::axesDisplay()
{
    if(axesAct->isChecked())
        enableAxes();
    else
        disableAxes();

    refresh();
}

void milxQtRenderWindow::disableOrient()
{
    humanAct->setChecked(false);
    refresh();
}

void milxQtRenderWindow::orientDisplay()
{
    humanAct->setChecked(!humanAct->isChecked());
    refresh();
}

void milxQtRenderWindow::lighting()
{
    if(lightingAct->isChecked())
        renderer->TwoSidedLightingOn();
    else
    {
        renderer->TwoSidedLightingOff();
        printInfo("Disabled Two-Sided Lighting.");
    }

    Render();
}

void milxQtRenderWindow::textDisplay()
{
    bool ok, ok1, ok2, ok3;
    QString newText = QInputDialog::getText(this, tr("Please Enter text"),
                                      tr("Text:"), QLineEdit::Normal, "A", &ok);
    float rValue = QInputDialog::getDouble(this, tr("Please Provide red channel (0.0-1.0)"),
                     tr("Red:"), 1.0, -2147483647, 2147483647, 3, &ok1);
    float gValue = QInputDialog::getDouble(this, tr("Please Provide green channel (0.0-1.0)"),
                     tr("Green:"), 0.0, -2147483647, 2147483647, 3, &ok2);
    float bValue = QInputDialog::getDouble(this, tr("Please Provide blue channel (0.0-1.0)"),
                     tr("Blue:"), 0.0, -2147483647, 2147483647, 3, &ok3);

    if(!ok || !ok1 || !ok2 || !ok3)
        return;

    vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
        textActor->SetInput(newText.toStdString().c_str());
        textActor->SetTextScaleModeToProp();
        textActor->GetTextProperty()->SetColor( rValue, gValue, bValue );
        textActors.append(textActor);

    vtkSmartPointer<vtkTextWidget> textWidget = vtkSmartPointer<vtkTextWidget>::New();
        textWidget->SetInteractor(QVTKWidget::GetInteractor());
        textWidget->SetTextActor(textActor);
        textWidget->On();
        textWidgets.append(textWidget);

    Connector->Connect(textWidget,
                       vtkCommand::WidgetActivateEvent,
                       this,
                       SLOT( updateTextActor(vtkObject*, unsigned long, void*, void*, vtkCommand*) )); //High Priority

    Render();
}

void milxQtRenderWindow::crosshair()
{
    printDebug("Toggling Crosshair.");
    Render();

    if(crosshairAct->isChecked())
        enableCrosshair();
    else
        disableCrosshair();
}

void milxQtRenderWindow::viewToXYPlane()
{
    if(rendered)
    {
        //viewToAxial
        vtkCamera *cam = milxQtRenderWindow::GetRenderer()->GetActiveCamera();
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            cam->SetFocalPoint(0, 0, 0);
            if(orientationView == RADIOLOGICAL)
            {
                cam->SetPosition(0,0,1);
                cam->SetViewUp(0,-1,0);
            }
            else //Neurological
            {
                cam->SetPosition(0,0,-1);
                cam->SetViewUp(0,1,0);
            }
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            milxQtRenderWindow::Render();

            currentView = AXIAL;
    }
}

void milxQtRenderWindow::viewToZXPlane()
{
    if(rendered)
    {
        //viewToSagittal
        vtkCamera *cam = milxQtRenderWindow::GetRenderer()->GetActiveCamera();
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            cam->SetFocalPoint(0, 0, 0);
            if(orientationView == RADIOLOGICAL)
            {
                cam->SetPosition(0,-1,0);
                cam->SetViewUp(0,0,1);
            }
            else //Neurological
            {
                cam->SetPosition(0,1,0);
                cam->SetViewUp(0,0,1);
            }
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            milxQtRenderWindow::Render();

            currentView = CORONAL;
    }
}

void milxQtRenderWindow::viewToZYPlane()
{
    if(rendered)
    {
        //viewToCoronal
        vtkCamera *cam = milxQtRenderWindow::GetRenderer()->GetActiveCamera();
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            cam->SetFocalPoint(0, 0, 0);
            cam->SetPosition(-1,0,0); //no orientation standard dependence
            cam->SetViewUp(0,0,1);
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            milxQtRenderWindow::Render();

            currentView = SAGITTAL;
    }
}

void milxQtRenderWindow::setView(int viewMode)
{
    printDebug("Changing view to " + QString::number(viewMode));
    if(viewMode == SAGITTAL) //sagittal
        viewToSagittal();
    else if(viewMode == CORONAL) //coronal
        viewToCoronal();
    else
        viewToAxial();
}

void milxQtRenderWindow::enableActionBasedOnView()
{
    if(currentView == CORONAL) //sagittal
        viewZX->setChecked(true);
    else if(currentView == SAGITTAL) //coronal
        viewZY->setChecked(true);
    else
        viewXY->setChecked(true);
}

void milxQtRenderWindow::saveView(QString filename)
{
    printInfo("Saving camera settings.");

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");

        settings.setValue("cameraSettings.ClippingRange.x", renderer->GetActiveCamera()->GetClippingRange()[0]);
        settings.setValue("cameraSettings.ClippingRange.y", renderer->GetActiveCamera()->GetClippingRange()[1]);
    //    file << renderer->GetActiveCamera()->GetClippingRange()[0] << " "
    //    << renderer->GetActiveCamera()->GetClippingRange()[1] << " "

        settings.setValue("cameraSettings.FocalPoint.x", renderer->GetActiveCamera()->GetFocalPoint()[0]);
        settings.setValue("cameraSettings.FocalPoint.y", renderer->GetActiveCamera()->GetFocalPoint()[1]);
        settings.setValue("cameraSettings.FocalPoint.z", renderer->GetActiveCamera()->GetFocalPoint()[2]);
    //    << renderer->GetActiveCamera()->GetFocalPoint()[0] << " "
    //    << renderer->GetActiveCamera()->GetFocalPoint()[1] << " "
    //    << renderer->GetActiveCamera()->GetFocalPoint()[2] << " "

        settings.setValue("cameraSettings.Position.x", renderer->GetActiveCamera()->GetPosition()[0]);
        settings.setValue("cameraSettings.Position.y", renderer->GetActiveCamera()->GetPosition()[1]);
        settings.setValue("cameraSettings.Position.z", renderer->GetActiveCamera()->GetPosition()[2]);
    //    << renderer->GetActiveCamera()->GetPosition()[0] << " "
    //    << renderer->GetActiveCamera()->GetPosition()[1] << " "
    //    << renderer->GetActiveCamera()->GetPosition()[2] << " "

        settings.setValue("cameraSettings.ViewAngle", renderer->GetActiveCamera()->GetViewAngle());
    //    << renderer->GetActiveCamera()->GetViewAngle() << " "

        settings.setValue("cameraSettings.ViewUp.x", renderer->GetActiveCamera()->GetViewUp()[0]);
        settings.setValue("cameraSettings.ViewUp.y", renderer->GetActiveCamera()->GetViewUp()[1]);
        settings.setValue("cameraSettings.ViewUp.z", renderer->GetActiveCamera()->GetViewUp()[2]);
    //    << renderer->GetActiveCamera()->GetViewUp()[0] << " "
    //    << renderer->GetActiveCamera()->GetViewUp()[1] << " "
    //    << renderer->GetActiveCamera()->GetViewUp()[2] << " "
    }
    else
    {
        milx::File::SaveCamera(filename.toStdString(), renderer->GetActiveCamera());
    }
}

void milxQtRenderWindow::saveViewFile()
{
    QFileDialog *fileSaver = new QFileDialog(this);
    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString();

    QString filename = fileSaver->getSaveFileName(this,
                                tr("Select File Name to Save"),
                                path,
                                tr("Camera Files (*.cam *.txt)"));

    if(filename.isEmpty())
        return;

    saveView(filename);
}

void milxQtRenderWindow::loadView(QString filename)
{
    printInfo("Laoding camera settings.");

    QFile file;
    file.setFileName(filename);
    if(!file.exists())
    {
        printWarning("File provided for loading did not exist. Using internally saved view instead");
        filename = "";
    }

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");

        double clippingRange[2];
        clippingRange[0] = settings.value("cameraSettings.ClippingRange.x").toDouble();
        clippingRange[1] = settings.value("cameraSettings.ClippingRange.y").toDouble();

        double focalPoint[3];
        focalPoint[0] = settings.value("cameraSettings.FocalPoint.x").toDouble();
        focalPoint[1] = settings.value("cameraSettings.FocalPoint.y").toDouble();
        focalPoint[2] = settings.value("cameraSettings.FocalPoint.z").toDouble();

        double position[3];
        position[0] = settings.value("cameraSettings.Position.x").toDouble();
        position[1] = settings.value("cameraSettings.Position.y").toDouble();
        position[2] = settings.value("cameraSettings.Position.z").toDouble();

        double angle = settings.value("cameraSettings.ViewAngle").toDouble();

        double viewUp[3];
        viewUp[0] = settings.value("cameraSettings.ViewUp.x").toDouble();
        viewUp[1] = settings.value("cameraSettings.ViewUp.y").toDouble();
        viewUp[2] = settings.value("cameraSettings.ViewUp.z").toDouble();

        camera->SetClippingRange(clippingRange);
        camera->SetFocalPoint(focalPoint);
        camera->SetPosition(position);
        camera->SetViewAngle(angle);
        camera->SetViewUp(viewUp);
    }
    else
    {
        milx::File::LoadCamera(filename.toStdString(), camera);
    }

    renderer->SetActiveCamera(camera);
    QVTKWidget::GetRenderWindow()->Render();
}

void milxQtRenderWindow::loadViewFile()
{
    QFileDialog *fileOpener = new QFileDialog(this);
    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString();

    QString filename = fileOpener->getOpenFileName(this,
                            tr("Select File to Open"),
                            path,
                            tr("Camera Files (*.cam *.txt)") );

    if(filename.isEmpty())
        return;

    loadView(filename);
}

void milxQtRenderWindow::addModelActor(vtkSmartPointer<vtkActor> mdlActor)
{
    AddActor(mdlActor); //!< transfer model actor to current render
    if(!rendered)
        generateRender();
}

void milxQtRenderWindow::removeModelActor(vtkSmartPointer<vtkActor> mdlActor)
{
    RemoveActor(mdlActor); //!< remove model actor to current render
    if(!rendered)
        generateRender();
}

void milxQtRenderWindow::addActor(vtkSmartPointer<vtkActor> actor, vtkMatrix4x4 *transformMatrix)
{
  vtkSmartPointer<vtkActor> rndActor = vtkSmartPointer<vtkActor>::New();
  rndActor->SetMapper(actor->GetMapper());
  if (transformMatrix)
    rndActor->SetUserMatrix(transformMatrix);
  AddActor(rndActor); //!< transfer actor to current render
}

void milxQtRenderWindow::addImageActor(vtkSmartPointer<vtkImageActor> imgActor, vtkMatrix4x4 *transformMatrix)
{
    //Check if already in view
    foreach(ImageActorItem item, imageActors)
    {
        if(imgActor == item.parentActor)
            return; //Dont need to do anything.
    }

    int extents[6];
    imgActor->GetDisplayExtent(extents);

    //Copy window level and actor to avoid leveling bugs like actor going blank or white when multiple actors in 3D view
    ImageActorItem item;
    item.parentActor = imgActor;
    item.imageActor = vtkSmartPointer<vtkImageActor>::New();
    item.imageActor->SetDisplayExtent(extents);
#if(VTK_MAJOR_VERSION > 5)
    item.imageActor->GetMapper()->SetInputData(imgActor->GetInput());
#else
    item.imageActor->SetInput(imgActor->GetInput());
#endif
    if(transformMatrix)
    {
//        printDebug("Transforming Actor to Image Space");
        item.imageActor->SetUserMatrix( transformMatrix );
        //Fix the actor properties, as they get reset to default ones. \todo fix this
    }
    ///Here we use the Qt signals and slots directly as it was found that the VTK-Qt connector caused problems
//    Connector->Connect(imgActor,
//                       vtkCommand::ModifiedEvent,
//                       this,
//                       SLOT( updateImageActor(vtkObject*, unsigned long, void*, void*, vtkCommand*) )); //High Priority

    //Added to image actor list
    imageActors.append(item);

    AddActor( item.imageActor );
    if(!rendered)
        generateRender();
}

void milxQtRenderWindow::updateTextActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
{
    vtkSmartPointer<vtkTextWidget> textWidget = vtkTextWidget::SafeDownCast(obj);
    vtkSmartPointer<vtkTextActor> textActor = textWidget->GetTextActor();

    QMessageBox msgBox;
    msgBox.setText("Would you like to remove this text?");
    msgBox.setInformativeText("Remove this text?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();

    if(ret == QMessageBox::No)
    {
        bool ok, ok1, ok2, ok3;
        QString newText = QInputDialog::getText(this, tr("Please Enter text"),
                                          tr("Text:"), QLineEdit::Normal, textActor->GetInput(), &ok);
        float rValue = QInputDialog::getDouble(this, tr("Please Provide red channel (0.0-1.0)"),
                         tr("Red:"), 1.0, -2147483647, 2147483647, 3, &ok1);
        float gValue = QInputDialog::getDouble(this, tr("Please Provide green channel (0.0-1.0)"),
                         tr("Green:"), 0.0, -2147483647, 2147483647, 3, &ok2);
        float bValue = QInputDialog::getDouble(this, tr("Please Provide blue channel (0.0-1.0)"),
                         tr("Blue:"), 0.0, -2147483647, 2147483647, 3, &ok3);

        if(!ok || !ok1 || !ok2 || !ok3)
            return;

        textActor->SetInput(newText.toStdString().c_str());
        textActor->SetTextScaleModeToProp();
        textActor->GetTextProperty()->SetColor( rValue, gValue, bValue );

        textWidget->SetTextActor(textActor);
    }
    else
        textWidget->Off();

    Render();
}

void milxQtRenderWindow::updateImageActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
{
  vtkSmartPointer<vtkImageActor> actor = vtkImageActor::SafeDownCast(obj);

  //Update the actor which is somehere in list
  foreach(ImageActorItem item, imageActors)
    {
      if(actor == item.parentActor)
        {
          int extents[6];

          actor->GetDisplayExtent(extents);

          item.imageActor->SetDisplayExtent(extents);
          item.imageActor->SetInterpolate(actor->GetInterpolate());
          //            item.imageActor->SetInput(actor->GetInput());

          break; //Only one of each is allowed
        }
    }

  Render();
}

void milxQtRenderWindow::updateImageActor(vtkSmartPointer<vtkImageActor> actor)
{
    //Update the actor which is somehere in list
    foreach(ImageActorItem item, imageActors)
    {
        if(actor == item.parentActor)
        {
            int extents[6];

            actor->GetDisplayExtent(extents);

            item.imageActor->SetDisplayExtent(extents);
            item.imageActor->SetInterpolate(actor->GetInterpolate());
//            item.imageActor->SetInput(actor->GetInput());

            break; //Only one of each is allowed
        }
    }

    Render();
}

void milxQtRenderWindow::userEvent(QMouseEvent *event)
{

}

void milxQtRenderWindow::userEvent(QKeyEvent *event)
{

}

void milxQtRenderWindow::userEvent(QWheelEvent *event)
{

}

void milxQtRenderWindow::refresh()
{
    //Check widgets
    if(lineAct->isChecked())
        lineWidget->On(); //enable lines
    else
        lineWidget->Off();
    if(distanceAct->isChecked())
      distanceWidget->On(); //enable distance measuring
    else
      distanceWidget->Off();
    if(biDirectionAct->isChecked())
        biDirectionWidget->On(); //enable cross distance measuring
    else
        biDirectionWidget->Off();
    if(angleAct->isChecked())
        angleWidget->On(); //enable angle measuring
    else
        angleWidget->Off();
    if(planeAct->isChecked())
        planeWidget->On(); //enable plane drawing
    else
        planeWidget->Off();
    if(boxAct->isChecked())
      boxWidget->On(); //enable box drawing
    else
      boxWidget->Off();
    if(sphereAct->isChecked())
        sphereWidget->On(); //enable sphere drawing
    else
        sphereWidget->Off();
    if(humanAct->isChecked())
        humanGlyph->On(); //enable orientation glyph
    else
        humanGlyph->Off();

    Render();
}

void milxQtRenderWindow::removeImageActor(vtkSmartPointer<vtkImageActor> imgActor)
{
    ImageActorItem foundItem;

    //Check if already in view
    size_t count = 0;
    foreach(ImageActorItem item, imageActors)
    {
        if(imgActor == item.imageActor || imgActor == item.parentActor)
        {
            foundItem.imageActor = item.imageActor; //Dont need to do anything.
            foundItem.parentActor = item.parentActor; //Dont need to do anything.
            imageActors.removeAt(count);
            break;
        }
        count ++;
    }

    RemoveActor( foundItem.imageActor ); //!< remove model actor to current render
    if(!rendered)
        generateRender();
}

void milxQtRenderWindow::importFrom(milxQtRenderWindow *windowToImportFrom)
{
    if(windowToImportFrom)
    {
        if(windowToImportFrom->GetImageActor()) //Is image so...
        {
            //Transform actor to correct orientation
            vtkSmartPointer<vtkMatrix4x4> transformModel = vtkSmartPointer<vtkMatrix4x4>::New();
            transformModel->DeepCopy(windowToImportFrom->getTransformMatrix());
            transformModel->Invert();

            addImageActor( windowToImportFrom->GetImageActor(), transformModel );

            ///Here we use the Qt signals and slots directly as it was found that the VTK-Qt connector caused problems
            ///with the image actors.
            connect(windowToImportFrom, SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), this, SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        }
        else
        {
            vtkActorCollection *actors = windowToImportFrom->GetActors();

            size_t n = actors->GetNumberOfItems();

            actors->InitTraversal();
            for(size_t j = 0; j < n; j ++)
            {
                printDebug("Adding Generic Actor into scene.");

                addModelActor(actors->GetNextItem()); //!< transfer model actor to current model

                qApp->processEvents(); ///Keep UI responsive
            }
        }
        Render();
    }
}

void milxQtRenderWindow::disableScale()
{
    if(!scaleBefore)
        return; //nothing to do then

    scalarBar->EnabledOff();
    scaleAct->setChecked(false);
    scaleBefore = false;
}

void milxQtRenderWindow::colourMapToRainbow(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toRainbow();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionRainbow->isChecked())
        actionRainbow->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToVTK(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toVTK();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionInvRainbow->isChecked())
        actionInvRainbow->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToJet(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toJet();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionJet->isChecked())
        actionJet->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToGray(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toGray();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionGray->isChecked())
        actionGray->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToSeismic(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
    colours->toSeismic();
    colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionSeismic->isChecked())
        actionSeismic->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToLogGray(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toLogGray();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionLogGray->isChecked())
        actionLogGray->setChecked(true);

    logScale = true;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToNIH(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toNIH();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionNIH->isChecked())
        actionNIH->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToNIH_Fire(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toNIH_FIRE();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionNIH_FIRE->isChecked())
        actionNIH_FIRE->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToAAL(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toAAL();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionAAL->isChecked())
        actionAAL->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToFS(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toFS(); //freesurfer
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionFS->isChecked())
        actionFS->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToHOT(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toHOT();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionHOT->isChecked())
        actionHOT->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToCOOL(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toCOOL();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionCOOL->isChecked())
        actionCOOL->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToCOOLWARM(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toCOOLWARM();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionCOOLWARM->isChecked())
        actionCOOLWARM->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToKnee(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toKnee();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionKnee->isChecked())
        actionKnee->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToBone(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toBone();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionBone->isChecked())
        actionBone->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToSpectral(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toSpectral();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionSpectral->isChecked())
        actionSpectral->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToGNUPlot(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toGNUPlot();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionGNUPlot->isChecked())
        actionGNUPlot->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToCubeHelix(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
        colours->toCubeHelix();
        colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionCubeHelix->isChecked())
        actionCubeHelix->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::colourMapToHSV(double minRange, double maxRange)
{
    double range[2];
    if(minRange == 0.0 && maxRange == 0.0)
        GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes
    else
    {
        range[0] = minRange;
        range[1] = maxRange;
    }

    milx::ColourMap *colours = new milx::ColourMap;
    colours->toHSV();
    colours->SetRange(range);

    lookupTable = colours->GetOutput();

    if(!actionHSV->isChecked())
        actionHSV->setChecked(true);

    logScale = false;
    updateLookupTable();
}

void milxQtRenderWindow::updateLookupTable()
{

}

void milxQtRenderWindow::reset()
{
    if(!rendered)
        return;

    renderer->ResetCamera();
    refresh();
}

void milxQtRenderWindow::generateRender()
{
    printDebug("Generating Render");

    ///Create the renderer for the actor
    if(!GetImageActor()) //not image then
        background();

    ///Create the Render Window with renderer
    if(QVTKWidget::GetRenderWindow()->GetRenderers()->GetNumberOfItems() == 0 || renderer != QVTKWidget::GetRenderWindow()->GetRenderers()->GetFirstRenderer()) //support only single renderer atm
        QVTKWidget::GetRenderWindow()->AddRenderer(renderer);

    ///Resize Window
    printDebug("Size of window: " + QString::number(QVTKWidget::size().height()) + "x" + QString::number(QVTKWidget::size().width()));
    printDebug("Actual Size of window: " + QString::number(QVTKWidget::GetRenderWindow()->GetSize()[0]) + "x" + QString::number(QVTKWidget::GetRenderWindow()->GetSize()[1]));
    if(QVTKWidget::size().height() < minWindowSize || QVTKWidget::size().width() < minWindowSize)
    {
        QVTKWidget::GetRenderWindow()->SetSize(minWindowSize, minWindowSize);
        printDebug("Resized to minimum size");
    }

    int *winSize = QVTKWidget::GetRenderWindow()->GetSize();
    QVTKWidget::resize(winSize[0], winSize[1]);

    reset();

    rendered = true;
    printDebug("Completed Generating Render");
}

void milxQtRenderWindow::createCustomConnections(QList<QAction *> actionsInMenu)
{
    foreach(QAction *action, actionsInMenu)
    {
        printDebug("Adding action for custom operations.");
        connect(action, SIGNAL(triggered()), this, SLOT(customOperation()));
    }
}

void milxQtRenderWindow::createCustomConnections(QMenu *fromMenu)
{
    QList<QAction *> actionsInMenu = fromMenu->actions();

    foreach(QAction *action, actionsInMenu)
    {
        printDebug("Adding action for custom operations.");
        connect(action, SIGNAL(triggered()), this, SLOT(customOperation()));
    }
}

void milxQtRenderWindow::customOperation()
{
    milxQtRenderWindow *currentWindow = qobject_cast<milxQtRenderWindow *>( windowActionGroup->checkedAction()->parent() );

    if(currentWindow == 0)
    {
        printError("Error in determining parent of action");
        return;
    }

    printDebug("Checked Window is: " + currentWindow->strippedNamePrefix());

    //Check for image actor
    if(currentWindow->GetImageActor()) //Is image so...
    {
        //Transform actor to correct orientation
        vtkSmartPointer<vtkMatrix4x4> transformModel = vtkSmartPointer<vtkMatrix4x4>::New();
        transformModel->DeepCopy(currentWindow->getTransformMatrix());
        transformModel->Invert();

        addImageActor( currentWindow->GetImageActor(), transformModel );

        ///Here we use the Qt signals and slots directly as it was found that the VTK-Qt connector caused problems
        ///with the image actors.
        connect(currentWindow, SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), this, SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
    }
    else
    {
        vtkActorCollection *actors = currentWindow->GetActors();

        size_t n = actors->GetNumberOfItems();

        actors->InitTraversal();
        for(size_t j = 0; j < n; j ++)
        {
            printDebug("Adding Generic Actor into scene.");

            addModelActor(actors->GetNextItem()); //!< transfer model actor to current model

            qApp->processEvents(); ///Keep UI responsive
        }
    }

    Render();
}

void milxQtRenderWindow::createMenu(QMenu *menu)
{
    if(!menu)
        return;

    menu->clear();

    foreach(QAction *currAct, milxQtWindow::actionsToAdd)
    {
      menu->addAction(currAct);
    }
    foreach(QMenu *currMenu, milxQtWindow::menusToAdd)
    {
      menu->addMenu(currMenu);
    }
    contourPolyDataAct->setDisabled(!contourAct->isChecked());
    contourInitAct->setDisabled(!contourAct->isChecked());
  //    contourImageAct->setDisabled(!contourAct->isChecked());

    menu->addMenu(contourMenu);
    menu->addAction(backgroundAct);
    menu->addAction(axesAct);
    menu->addAction(lightingAct);
    menu->addAction(lineAct);
    menu->addAction(distanceAct);
    menu->addAction(biDirectionAct);
    menu->addAction(angleAct);
    menu->addAction(planeAct);
    menu->addAction(boxAct);
    menu->addAction(sphereAct);
    menu->addAction(humanAct);
    menu->addAction(textAct);

    ///Change View of Volume
    menu->addMenu(viewMenu);
    viewMenu->addAction(viewXY);
    viewMenu->addAction(viewZX);
    viewMenu->addAction(viewZY);
    viewMenu->addAction(saveViewAct);
    viewMenu->addAction(saveViewFileAct);
    viewMenu->addAction(loadViewAct);
    viewMenu->addAction(loadViewFileAct);
    enableActionBasedOnView();
    colourMapsMenu();

    menu->addSeparator();
    foreach(QAction *currAct, milxQtWindow::actionsToAppend)
    {
      menu->addAction(currAct);
    }
    foreach(QMenu *currMenu, milxQtWindow::menusToAppend)
    {
      menu->addMenu(currMenu);
    }

    menu->addAction(refreshAct);
    menu->addAction(resetAct);
}

bool milxQtRenderWindow::openModelUsingQt(const QString filename, vtkSmartPointer<vtkPolyData> &data)
{
    QFile qtFile(filename);
    if(!qtFile.open(QIODevice::ReadOnly))
    {
        printError("Unable to open model using Qt for reading.");
        return false;
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
        normals->SetNumberOfComponents(3);

    QTextStream inFile(&qtFile);
    while (!inFile.atEnd())
    {
        QString lineType;
        inFile >> lineType;

        if(lineType == "" || lineType == "#" || lineType == "mtllib" || lineType == "usemtl" || lineType == "g" || lineType == "s" || lineType == "vt") //no material, texture or scalars support yet
            continue;
        else if(lineType == "v") //vertices v
        {
            double position[3];
            inFile >> position[0] >> position[1] >> position[2];
//            qDebug() << position[0] << "," << position[1] << "," << position[2];
            points->InsertNextPoint(position);
        }
        else if(lineType == "vn") //normals vn
        {
            float normal[3];
            inFile >> normal[0] >> normal[1] >> normal[2];
//            qDebug() << normal[0] << "," << normal[1] << "," << normal[2];
		#if VTK_MAJOR_VERSION > 7
			normals->InsertNextTypedTuple(normal); //InsertNextTypedTuple
		#else
			normals->InsertNextTupleValue(normal); //InsertNextTypedTuple
		#endif
        }
        else if(lineType == "f") //faces/cells
        {
            QString face[3]; ///tri's only \todo generalise
            inFile >> face[0] >> face[1] >> face[2];
            // v/vt/vn format read

            cells->InsertNextCell(3);
            for(size_t j = 0; j < 3; j ++)
            {
                QStringList faceParts = face[j].split('/');
//                qDebug() << faceParts;
                cells->InsertCellPoint(faceParts[0].toInt()-1); //obj indices indexed from 1 not 0
            }
        }
    }

    if(!data)
        data = vtkSmartPointer<vtkPolyData>::New();
        data->SetPoints(points);
        data->SetPolys(cells);
        data->GetPointData()->SetNormals(normals);

    return true;
}

void milxQtRenderWindow::enableUpdates(QStatusBar *bar)
{
    if(!dataPicker)
        dataPicker = vtkSmartPointer<vtkPointPicker>::New();

    Connector->Connect( QVTKWidget::GetRenderWindow()->GetInteractor(),
                        vtkCommand::MouseMoveEvent,
                        this,
                        SLOT( updateCoords(vtkObject*) ) );

    updateBar = bar;
}

void milxQtRenderWindow::updateCoords(vtkObject *obj)
{
    ///Get interactor
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    QString message = "";

    ///Get event position
    ///Code initial by Mark Wyszomierski 2003-2007 @ devsample
    ///Modified by Shekhar Chandra
    ///Do the pick. It will return a non-zero value if we intersected the image.
    if (dataPicker->Pick(iren->GetEventPosition()[0],
                         iren->GetEventPosition()[1],
                         0,  // always zero.
                         renderer))
    {
        // Get the mapped position of the mouse using the picker.
        double point[3];
        dataPicker->GetPickPosition(point);

        // Get the volume index within the entire volume now.
        vtkIdType nVolIdx = dataPicker->GetPointId();

        message = "Point " + QString::number(nVolIdx) + ": (" + QString::number(point[0]) + ", " + QString::number(point[1]) + ", " + QString::number(point[2]) + ")";
    }

    ///Write message to status bar
    updateBar->showMessage(message);
}

void milxQtRenderWindow::createActions()
{
    contextMenu = new QMenu(this);
    contextMenu->setTitle(QApplication::translate("MainWindow", "Rendering", 0, QApplication::UnicodeUTF8));

    contourAct = new QAction(this);
    contourAct->setText(QApplication::translate("Render", "Contour Mode", 0, QApplication::UnicodeUTF8));
    contourAct->setCheckable(true);
    contourAct->setChecked(false);
    contourPolyDataAct = new QAction(this);
    contourPolyDataAct->setText(QApplication::translate("Render", "Create Model From Contour", 0, QApplication::UnicodeUTF8));
    contourNodePolyDataAct = new QAction(this);
    contourNodePolyDataAct->setText(QApplication::translate("Render", "Save Contour", 0, QApplication::UnicodeUTF8));
    contourInitAct = new QAction(this);
    contourInitAct->setText(QApplication::translate("Render", "Load Contour", 0, QApplication::UnicodeUTF8));
//    contourImageAct = new QAction(this);
//    contourImageAct->setText(QApplication::translate("Render", "Create Binary Image From Contour", 0, QApplication::UnicodeUTF8));

    contourMenu = new QMenu(this);
    contourMenu->setTitle(QApplication::translate("Render", "Contouring", 0, QApplication::UnicodeUTF8));
    contourMenu->addAction(contourAct);
    contourMenu->addAction(contourPolyDataAct);
    contourMenu->addAction(contourNodePolyDataAct);
    contourMenu->addAction(contourInitAct);
//    contourMenu->addAction(contourImageAct);

    backgroundAct = new QAction(this);
    backgroundAct->setText(QApplication::translate("Render", "&White Background", 0, QApplication::UnicodeUTF8));
    backgroundAct->setCheckable(true);

    axesAct = new QAction(this);
    axesAct->setText(QApplication::translate("Render", "&Axes", 0, QApplication::UnicodeUTF8));
    axesAct->setCheckable(true);

    lightingAct = new QAction(this);
    lightingAct->setText(QApplication::translate("Render", "&Two-Sided Lighting", 0, QApplication::UnicodeUTF8));
    lightingAct->setCheckable(true);
    lightingAct->setChecked(true);

    lineAct = new QAction(this);
    lineAct->setText(QApplication::translate("Render", "&Draw Line", 0, QApplication::UnicodeUTF8));
    lineAct->setCheckable(true);
    lineAct->setChecked(false);

    distanceAct = new QAction(this);
    distanceAct->setText(QApplication::translate("Render", "&Measure Distance", 0, QApplication::UnicodeUTF8));
    distanceAct->setCheckable(true);
    distanceAct->setChecked(false);

    biDirectionAct = new QAction(this);
    biDirectionAct->setText(QApplication::translate("Render", "&Measure Cross Distance", 0, QApplication::UnicodeUTF8));
    biDirectionAct->setCheckable(true);
    biDirectionAct->setChecked(false);

    angleAct = new QAction(this);
    angleAct->setText(QApplication::translate("Render", "&Measure Angle", 0, QApplication::UnicodeUTF8));
    angleAct->setCheckable(true);
    angleAct->setChecked(false);

    planeAct = new QAction(this);
    planeAct->setText(QApplication::translate("Render", "&Draw Plane", 0, QApplication::UnicodeUTF8));
    planeAct->setCheckable(true);
    planeAct->setChecked(false);

    boxAct = new QAction(this);
    boxAct->setText(QApplication::translate("Render", "&Draw Box/Cuboid", 0, QApplication::UnicodeUTF8));
    boxAct->setCheckable(true);
    boxAct->setChecked(false);

    sphereAct = new QAction(this);
    sphereAct->setText(QApplication::translate("Render", "&Draw Sphere/Circle", 0, QApplication::UnicodeUTF8));
    sphereAct->setCheckable(true);
    sphereAct->setChecked(false);

    humanAct = new QAction(this);
    humanAct->setText(QApplication::translate("Render", "Show Human Orientation", 0, QApplication::UnicodeUTF8));
    humanAct->setCheckable(true);
    humanAct->setChecked(true);

    textAct = new QAction(this);
    textAct->setText(QApplication::translate("Render", "&Insert Text", 0, QApplication::UnicodeUTF8));

    crosshairAct = new QAction(this);
    crosshairAct->setText(QApplication::translate("Render", "&Crosshair", 0, QApplication::UnicodeUTF8));
    crosshairAct->setCheckable(true);
    crosshairAct->setChecked(false);

    windowPropertiesMenu = new QMenu(this);
    windowPropertiesMenu->setTitle(QApplication::translate("Render", "Window Properties", 0, QApplication::UnicodeUTF8));
    windowPropertiesMenu->addAction(backgroundAct);
    windowPropertiesMenu->addAction(axesAct);
    windowPropertiesMenu->addAction(lightingAct);
    windowPropertiesMenu->addAction(lineAct);
    windowPropertiesMenu->addAction(distanceAct);
    windowPropertiesMenu->addAction(biDirectionAct);
    windowPropertiesMenu->addAction(angleAct);
    windowPropertiesMenu->addAction(planeAct);
    windowPropertiesMenu->addAction(boxAct);
    windowPropertiesMenu->addAction(sphereAct);
    windowPropertiesMenu->addAction(textAct);
    windowPropertiesMenu->addAction(crosshairAct);

    viewMenu = new QMenu(this);
    viewMenu->setTitle("View");
    viewXY = new QAction(this);
    viewXY->setText(QApplication::translate("Render", "Axial (xy-plane)", 0, QApplication::UnicodeUTF8));
    viewXY->setShortcut(tr("Alt+a"));
    viewXY->setCheckable(true);
    viewXY->setChecked(true);
    viewZX = new QAction(this);
    viewZX->setText(QApplication::translate("Render", "Coronal (zx-plane)", 0, QApplication::UnicodeUTF8));
    viewZX->setShortcut(tr("Alt+c"));
    viewZX->setCheckable(true);
    viewZY = new QAction(this);
    viewZY->setText(QApplication::translate("Render", "Sagittal (zy-plane)", 0, QApplication::UnicodeUTF8));
    viewZY->setShortcut(tr("Alt+s"));
    viewZY->setCheckable(true);
    saveViewAct = new QAction(this);
    saveViewAct->setText(QApplication::translate("Render", "Save View", 0, QApplication::UnicodeUTF8));
    saveViewAct->setShortcut(tr("Ctrl+Alt+s"));
    loadViewAct = new QAction(this);
    loadViewAct->setText(QApplication::translate("Render", "Load View", 0, QApplication::UnicodeUTF8));
    loadViewAct->setShortcut(tr("Ctrl+Alt+l"));
    saveViewFileAct = new QAction(this);
    saveViewFileAct->setText(QApplication::translate("Render", "Save View to File", 0, QApplication::UnicodeUTF8));
    saveViewFileAct->setShortcut(tr("Shift+Alt+s"));
    loadViewFileAct = new QAction(this);
    loadViewFileAct->setText(QApplication::translate("Render", "Load View from File", 0, QApplication::UnicodeUTF8));
    loadViewFileAct->setShortcut(tr("Shift+Alt+l"));
    scaleAct = new QAction(this);
    scaleAct->setText(QApplication::translate("Render", "&Scalar Bar", 0, QApplication::UnicodeUTF8));
    scaleAct->setCheckable(true);
    viewGroup = new QActionGroup(this);
    viewGroup->addAction(viewXY);
    viewGroup->addAction(viewZX);
    viewGroup->addAction(viewZY);

    //Colour map
    colourMapMenu = new QMenu(this);
    colourMapMenu->setTitle("Colour Maps");
    actionDefault = new QAction(this);
    actionDefault->setText(QApplication::translate("Render", "None", 0, QApplication::UnicodeUTF8));
//    actionDefault->setDisabled(false);
    actionDefault->setCheckable(true);
    actionDefault->setChecked(true);
    actionJet = new LabelledAction("Jet", QPixmap(":resources/cmaps/cm_jet.png"));
//    actionJet->setText(QApplication::translate("Render", "Jet", 0, QApplication::UnicodeUTF8));
    actionJet->setCheckable(true);
    actionJet->setChecked(false);
    actionRainbow = new LabelledAction("Rainbow", QPixmap(":resources/cmaps/cm_rainbow.png"));
//    actionRainbow->setText(QApplication::translate("Render", "Rainbow", 0, QApplication::UnicodeUTF8));
    actionRainbow->setCheckable(true);
    actionRainbow->setChecked(false);
    actionInvRainbow = new LabelledAction("VTK", QPixmap(":resources/cmaps/cm_vtk.png"), this);
//    actionInvRainbow->setText(QApplication::translate("Render", "VTK", 0, QApplication::UnicodeUTF8));
    actionInvRainbow->setCheckable(true);
    actionInvRainbow->setChecked(false);
    actionGray = new LabelledAction("Gray", QPixmap(":resources/cmaps/cm_gray.png"));
//    actionGray->setText(QApplication::translate("Render", "Gray", 0, QApplication::UnicodeUTF8));
    actionGray->setCheckable(true);
    actionGray->setChecked(false);
    actionLogGray = new LabelledAction("Gray (Log)", QPixmap(":resources/cmaps/cm_gray.png"));
//    actionLogGray->setText(QApplication::translate("Render", "Gray (Log 10)", 0, QApplication::UnicodeUTF8));
    actionLogGray->setCheckable(true);
    actionLogGray->setChecked(false);
    actionSeismic = new LabelledAction("Seismic", QPixmap(":resources/cmaps/cm_seismic.png"));
    actionSeismic->setCheckable(true);
    actionSeismic->setChecked(false);
    actionNIH = new LabelledAction("NIH", QPixmap(":resources/cmaps/cm_nih.png"), this);
//    actionNIH->setText(QApplication::translate("Render", "NIH", 0, QApplication::UnicodeUTF8));
    actionNIH->setCheckable(true);
    actionNIH->setChecked(false);
    actionNIH_FIRE = new LabelledAction("NIH Fire", QPixmap(":resources/cmaps/cm_nih_fire.png"), this);
//    actionNIH_FIRE->setText(QApplication::translate("Render", "NIH Fire", 0, QApplication::UnicodeUTF8));
    actionNIH_FIRE->setCheckable(true);
    actionNIH_FIRE->setChecked(false);
    actionAAL = new LabelledAction("AAL", QPixmap(":resources/cmaps/cm_aal.png"), this);
//    actionAAL->setText(QApplication::translate("Render", "AAL", 0, QApplication::UnicodeUTF8));
    actionAAL->setCheckable(true);
    actionAAL->setChecked(false);
    actionFS = new LabelledAction("FreeSurfer", QPixmap(":resources/cmaps/cm_fs.png"), this);
//    actionFS->setText(QApplication::translate("Render", "FreeSurfer", 0, QApplication::UnicodeUTF8));
    actionFS->setCheckable(true);
    actionFS->setChecked(false);
    actionHOT = new LabelledAction("HOT", QPixmap(":resources/cmaps/cm_hot.png"), this);
//    actionHOT->setText(QApplication::translate("Render", "HOT", 0, QApplication::UnicodeUTF8));
//    actionHOT->setIcon(QIcon(":resources/cmaps/cm_hot.png"));
    actionHOT->setCheckable(true);
    actionHOT->setChecked(false);
    actionCOOL = new LabelledAction("COOL", QPixmap(":resources/cmaps/cm_cool.png"), this);
//    actionCOOL->setText(QApplication::translate("Render", "COOL", 0, QApplication::UnicodeUTF8));
    actionCOOL->setCheckable(true);
    actionCOOL->setChecked(false);
    actionCOOLWARM = new LabelledAction("COOL-WARM", QPixmap(":resources/cmaps/cm_coolwarm.png"), this);
//    actionCOOLWARM->setText(QApplication::translate("Render", "COOL-WARM", 0, QApplication::UnicodeUTF8));
//    actionCOOLWARM->setIcon(QIcon(":resources/cmaps/cm_coolwarm.png"));
    actionCOOLWARM->setCheckable(true);
    actionCOOLWARM->setChecked(false);
    actionKnee = new LabelledAction("Knee", QPixmap(":resources/cmaps/cm_knee.png"), this);
//    actionKnee->setText(QApplication::translate("Render", "Knee", 0, QApplication::UnicodeUTF8));
    actionKnee->setCheckable(true);
    actionKnee->setChecked(false);
    actionBone = new LabelledAction("Bone", QPixmap(":resources/cmaps/cm_bone.png"), this);
//    actionBone->setText(QApplication::translate("Render", "Bone", 0, QApplication::UnicodeUTF8));
//    actionBone->setIcon(QIcon(":resources/cmaps/cm_bone.png"));
    actionBone->setCheckable(true);
    actionBone->setChecked(false);
    actionSpectral = new LabelledAction("Spectral", QPixmap(":resources/cmaps/cm_Spectral.png"), this);
//    actionSpectral->setText(QApplication::translate("Render", "Spectral", 0, QApplication::UnicodeUTF8));
//    actionSpectral->setIcon(QIcon(":resources/cmaps/cm_Spectral.png"));
    actionSpectral->setCheckable(true);
    actionSpectral->setChecked(false);
    actionGNUPlot = new LabelledAction("GNUPlot", QPixmap(":resources/cmaps/cm_gnuplot.png"), this);
//    actionGNUPlot->setText(QApplication::translate("Render", "GNUPlot", 0, QApplication::UnicodeUTF8));
//    actionGNUPlot->setIcon(QIcon(":resources/cmaps/cm_gnuplot.png"));
    actionGNUPlot->setCheckable(true);
    actionGNUPlot->setChecked(false);
    actionCubeHelix = new LabelledAction("CubeHelix", QPixmap(":resources/cmaps/cm_cubehelix.png"), this);
//    actionCubeHelix->setText(QApplication::translate("Render", "CubeHelix", 0, QApplication::UnicodeUTF8));
//    actionCubeHelix->setIcon(QIcon(":resources/cmaps/cm_gnuplot.png"));
    actionCubeHelix->setCheckable(true);
    actionCubeHelix->setChecked(false);
    actionHSV = new LabelledAction("HSV", QPixmap(":resources/cmaps/cm_hsv.png"), this);
    actionHSV->setCheckable(true);
    actionHSV->setChecked(false);
    mapGroup = new QActionGroup(this);
    mapGroup->addAction(actionDefault);
    mapGroup->addAction(actionJet);
    mapGroup->addAction(actionRainbow);
    mapGroup->addAction(actionInvRainbow);
    mapGroup->addAction(actionGray);
    mapGroup->addAction(actionLogGray);
    mapGroup->addAction(actionSeismic);
    mapGroup->addAction(actionNIH);
    mapGroup->addAction(actionNIH_FIRE);
    mapGroup->addAction(actionAAL);
    mapGroup->addAction(actionFS);
    mapGroup->addAction(actionHOT);
    mapGroup->addAction(actionCOOL);
    mapGroup->addAction(actionCOOLWARM);
    mapGroup->addAction(actionKnee);
    mapGroup->addAction(actionBone);
    mapGroup->addAction(actionSpectral);
    mapGroup->addAction(actionGNUPlot);
    colourMapsMenu();

    refreshAct = new QAction(this);
    refreshAct->setText(QApplication::translate("Render", "&Refresh", 0, QApplication::UnicodeUTF8));
    refreshAct->setShortcut(tr("F5"));
    resetAct = new QAction(this);
    resetAct->setText(QApplication::translate("Render", "Reset", 0, QApplication::UnicodeUTF8));
    resetAct->setIcon(QIcon(":/resources/toolbar/refresh.png"));
    resetAct->setShortcut(tr("F7"));
}

void milxQtRenderWindow::createConnections()
{
    //Re-direct mouse events
    connect(this, SIGNAL(mouseEvent(QMouseEvent*)), this, SLOT(userEvent(QMouseEvent*)));

    //Stop VTK intercepting right click events
    /*Connector->Connect(QVTKWidget::GetInteractor(),
                       vtkCommand::RightButtonPressEvent,
                       this,
                       SLOT( consumeVTKEvent(vtkObject*, unsigned long, void*, void*, vtkCommand*) ),
                       NULL, 1.0); //High Priority
    Connector->Connect(QVTKWidget::GetInteractor(),
                       vtkCommand::RightButtonReleaseEvent,
                       this,
                       SLOT( consumeVTKEvent(vtkObject*, unsigned long, void*, void*, vtkCommand*) ),
                       NULL, 1.0); //High Priority*/

    connect(contourAct, SIGNAL(triggered()), this, SLOT(contour()));
    connect(contourPolyDataAct, SIGNAL(triggered()), this, SLOT(contourAsPolyData()));
    connect(contourNodePolyDataAct, SIGNAL(triggered()), this, SLOT(contourAsNodePolyData()));
    connect(contourInitAct, SIGNAL(triggered()), this, SLOT(contourInitFromPolyData()));
//    connect(contourImageAct, SIGNAL(triggered()), this, SLOT(contourAsBinaryImage()));
    connect(backgroundAct, SIGNAL(triggered()), this, SLOT(background()));
    connect(axesAct, SIGNAL(triggered()), this, SLOT(axesDisplay()));
    connect(lightingAct, SIGNAL(triggered()), this, SLOT(lighting()));
    connect(lineAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(distanceAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(biDirectionAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(angleAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(planeAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(boxAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(sphereAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(humanAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(textAct, SIGNAL(triggered()), this, SLOT(textDisplay()));
    connect(crosshairAct, SIGNAL(triggered()), this, SLOT(crosshair()));
    connect(viewXY, SIGNAL(triggered()), this, SLOT(viewToXYPlane()));
    connect(viewZX, SIGNAL(triggered()), this, SLOT(viewToZXPlane()));
    connect(viewZY, SIGNAL(triggered()), this, SLOT(viewToZYPlane()));
    connect(actionJet, SIGNAL(triggered()), this, SLOT(colourMapToJet()));
    connect(actionRainbow, SIGNAL(triggered()), this, SLOT(colourMapToRainbow()));
    connect(actionInvRainbow, SIGNAL(triggered()), this, SLOT(colourMapToVTK()));
    connect(actionGray, SIGNAL(triggered()), this, SLOT(colourMapToGray()));
    connect(actionLogGray, SIGNAL(triggered()), this, SLOT(colourMapToLogGray()));
    connect(actionSeismic, SIGNAL(triggered()), this, SLOT(colourMapToSeismic()));
    connect(actionNIH, SIGNAL(triggered()), this, SLOT(colourMapToNIH()));
    connect(actionNIH_FIRE, SIGNAL(triggered()), this, SLOT(colourMapToNIH_Fire()));
    connect(actionAAL, SIGNAL(triggered()), this, SLOT(colourMapToAAL()));
    connect(actionFS, SIGNAL(triggered()), this, SLOT(colourMapToFS()));
    connect(actionHOT, SIGNAL(triggered()), this, SLOT(colourMapToHOT()));
    connect(actionCOOL, SIGNAL(triggered()), this, SLOT(colourMapToCOOL()));
    connect(actionCOOLWARM, SIGNAL(triggered()), this, SLOT(colourMapToCOOLWARM()));
    connect(actionKnee, SIGNAL(triggered()), this, SLOT(colourMapToKnee()));
    connect(actionBone, SIGNAL(triggered()), this, SLOT(colourMapToBone()));
    connect(actionSpectral, SIGNAL(triggered()), this, SLOT(colourMapToSpectral()));
    connect(actionGNUPlot, SIGNAL(triggered()), this, SLOT(colourMapToGNUPlot()));
    connect(actionCubeHelix, SIGNAL(triggered()), this, SLOT(colourMapToCubeHelix()));
    connect(actionHSV, SIGNAL(triggered()), this, SLOT(colourMapToHSV()));
    connect(saveViewAct, SIGNAL(triggered()), this, SLOT(saveView()));
    connect(loadViewAct, SIGNAL(triggered()), this, SLOT(loadView()));
    connect(saveViewFileAct, SIGNAL(triggered()), this, SLOT(saveViewFile()));
    connect(loadViewFileAct, SIGNAL(triggered()), this, SLOT(loadViewFile()));
    connect(scaleAct, SIGNAL(triggered()), this, SLOT(scaleDisplay()));
    connect(refreshAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(resetAct, SIGNAL(triggered()), this, SLOT(reset()));
}

void milxQtRenderWindow::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    createMenu(contextMenu); //!< Only exists for the duration of the context selection

    contextMenu->exec(currentEvent->globalPos());
}

QMenu* milxQtRenderWindow::colourMapsMenu()
{
    colourMapMenu->addAction(actionDefault);
    colourMapMenu->addAction(actionJet);
    colourMapMenu->addAction(actionRainbow);
    colourMapMenu->addAction(actionGray);
    colourMapMenu->addAction(actionLogGray);
    colourMapMenu->addAction(actionSeismic);
    colourMapMenu->addAction(actionSpectral);
    colourMapMenu->addAction(actionGNUPlot);
    colourMapMenu->addAction(actionBone);
    colourMapMenu->addAction(actionHOT);
    colourMapMenu->addAction(actionCOOL);
    colourMapMenu->addAction(actionCOOLWARM);
    colourMapMenu->addAction(actionCubeHelix);
    colourMapMenu->addAction(actionHSV);
    colourMapMenu->addAction(actionInvRainbow);
    colourMapMenu->addAction(actionNIH);
    colourMapMenu->addAction(actionNIH_FIRE);
    colourMapMenu->addAction(actionAAL);
    colourMapMenu->addAction(actionFS);
    colourMapMenu->addAction(actionKnee);

    return colourMapMenu;
}

//Drag members
//Drop members
void milxQtRenderWindow::dragLeaveEvent(QDragLeaveEvent *event)
{
    event->accept();
}

void milxQtRenderWindow::dragMoveEvent(QDragMoveEvent *event)
{
    event->acceptProposedAction();
}

void milxQtRenderWindow::dragEnterEvent(QDragEnterEvent *currentEvent)
{
    if(currentEvent->mimeData()->hasFormat("text/uri-list") || currentEvent->mimeData()->hasFormat("text/plain"))
        currentEvent->acceptProposedAction();
}

void milxQtRenderWindow::mouseDoubleClickEvent(QMouseEvent *event)
{
    QDrag *drag = new QDrag(this);
    QMimeData *mimeData = new QMimeData;

    QList<QUrl> urls;
    urls.append( QUrl("file://localhost/" + getName()) );

    mimeData->setUrls(urls);
    drag->setMimeData(mimeData);

    drag->exec(Qt::CopyAction | Qt::MoveAction | Qt::LinkAction, Qt::CopyAction);
}

void milxQtRenderWindow::dropEvent(QDropEvent *currentEvent)
{
    if(currentEvent->source() == this) //disallow self drops
        return;

    milxQtRenderWindow *sourceWindow = qobject_cast<milxQtRenderWindow *>(currentEvent->source());

    QList<QUrl> urlsList = currentEvent->mimeData()->urls();
    QString tmp, typeString = currentEvent->mimeData()->text();

    for(int j = 0; j < urlsList.size(); j ++)
    {
        if(urlsList[j].isValid())
        {
            printInfo("Dropped Path into Render Window: " + urlsList[j].path());

            if(sourceWindow)
            {
                importFrom(sourceWindow);
            }
            else
            {
                printError("Data type unsupported. Passing on drop... ");
                currentEvent->ignore(); //dont accept drop
            }
        }
    }

    currentEvent->acceptProposedAction();
}

void milxQtRenderWindow::setupHumanGlyph(vtkSmartPointer<vtkMatrix4x4> mat)
{
    vtkSmartPointer<vtkPolyData> human;
    openModelUsingQt(":/resources/human.obj", human); //needs to be done usingQt method to access model as a resource

    vtkSmartPointer<vtkPolyData> humanTransformed;
    if(mat)
    {
        vtkSmartPointer<vtkTransform> transf = vtkSmartPointer<vtkTransform>::New();
            transf->SetMatrix(mat);

        vtkSmartPointer<vtkTransformPolyDataFilter> transformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        #if VTK_MAJOR_VERSION <= 5
            transformer->SetInput(human);
        #else
            transformer->SetInputData(human);
        #endif
            transformer->SetTransform(transf);
            transformer->Update();

        vtkSmartPointer<vtkPolyDataNormals> normalsHuman = vtkSmartPointer<vtkPolyDataNormals>::New();
        #if VTK_MAJOR_VERSION <= 5
            normalsHuman->SetInput(transformer->GetOutput());
        #else
            normalsHuman->SetInputData(transformer->GetOutput());
        #endif
            normalsHuman->Update();

        humanTransformed = normalsHuman->GetOutput();
    }
    else
        humanTransformed = human;

    vtkSmartPointer<vtkPropAssembly> propAssembly = vtkSmartPointer<vtkPropAssembly>::New();
    orientAxes = vtkSmartPointer<vtkAxesActor>::New();

    if(orientationAxes)
    {
        orientAxes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->ShadowOff();
        orientAxes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(14);
        orientAxes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
        orientAxes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->ShadowOff();
        orientAxes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(14);
        orientAxes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
        orientAxes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->ShadowOff();
        orientAxes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(14);
        orientAxes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
        orientAxes->SetXAxisLabelText("Left");
        orientAxes->SetYAxisLabelText("P");
        orientAxes->SetZAxisLabelText("S");

        if(mat)
            orientAxes->SetUserMatrix(mat);

        propAssembly->AddPart(orientAxes);
    }

    vtkSmartPointer<vtkPolyDataMapper> humanMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
        humanMapper->SetInput(humanTransformed);
    #else
        humanMapper->SetInputData(humanTransformed);
    #endif
    vtkSmartPointer<vtkActor> humanActor = vtkSmartPointer<vtkActor>::New();
        humanActor->SetMapper(humanMapper);
        humanActor->GetProperty()->SetInterpolationToPhong();
//        humanActor->GetProperty()->SetColor(235/255.0, 180/255.0, 173/255.0); //skin
        humanActor->GetProperty()->SetColor(60/255.0, 232/255.0, 30/255.0); //Varian Eclipse green
    propAssembly->AddPart(humanActor);
    humanGlyph = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        humanGlyph->SetOutlineColor(0.0, 0, 0.0);
        humanGlyph->SetOrientationMarker(propAssembly);
        humanGlyph->SetInteractor(QVTKWidget::GetInteractor());
        humanGlyph->SetViewport(0.0, 0.0, 0.3, 0.3);
        humanGlyph->SetDefaultRenderer(renderer);
}
