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
#include <limits>

#include "milxQtUnifiedWindow.h"
#include <QDebug>
//ITK
#include <itkLinearInterpolateImageFunction.h>
//VTK
#include <vtkSmartPointer.h>
//#include <vtkPolyDataPointSampler.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <vtkLineSource.h> //Testing
#include <vtkDoubleArray.h>
//Checkboarding
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkCheckerboardRepresentation.h>

milxQtUnifiedWindow::milxQtUnifiedWindow(QWidget *theParent) : milxQtRenderWindow(theParent)
{
    ///Set strings
    milxQtWindow::prefix = "Uni: ";

    createActions();

    createConnections();
}

milxQtUnifiedWindow::~milxQtUnifiedWindow()
{
    //dtor
}

void milxQtUnifiedWindow::addToWindow(milxQtModel *model)
{
    bool alreadyPresent = false;

    if(!model)
        return;

    setCommonProperties(model);

    if(!unifyModels.isEmpty())
        alreadyPresent = unifyModels.contains(model);

    if(!alreadyPresent)
        unifyModels.append(model);

    refresh();
}

void milxQtUnifiedWindow::addToWindow(milxQtImage *image)
{
    bool alreadyPresent = false;

    if(!image)
        return;

    setCommonProperties(image);

    if(!unifyImages.isEmpty())
        alreadyPresent = unifyImages.contains(image);

    if(!alreadyPresent)
        unifyImages.append(image);

    refresh();
}

void milxQtUnifiedWindow::removeFromWindow(QWidget *passedWindow)
{
    bool present = false;

    if(!passedWindow)
        return;

    ///Up cast
    milxQtModel *model = qobject_cast<milxQtModel *>(passedWindow);
//    if(model == 0)
//        milxQtImage *image = qobject_cast<milxQtImage *>(passedWindow);

    if(!unifyModels.isEmpty())
        present = unifyModels.contains(model);

    if(present)
    {
        unifyModels.removeAll(model);
        RemoveActor(model->GetActor()); //!< Remove data from general display
        unionAct->setChecked(true); //reset the unified window display type just in case
    }
}

void milxQtUnifiedWindow::refresh()
{
    if(unionAct->isChecked())
        generateUnion();
    if(geoDifferenceAct->isChecked())
        generateDifference();
}

void milxQtUnifiedWindow::generateUnion()
{
    if(unifyModels.isEmpty() && unifyImages.isEmpty())
        return;

    printInfo("Generating Union");
    RemoveAllActors();
    foreach(milxQtModel *model, unifyModels) //Add models
    {
        if(unionAct->isChecked())
            addModelActor(model->GetActor()); //!< Add data to general display
        else
            removeModelActor(model->GetActor()); //!< Remove data from general display

        qApp->processEvents(); ///Keep UI responsive
    }
    foreach(milxQtImage *img, unifyImages) //Add images
    {
        if(unionAct->isChecked())
        {
            //Transform actor to correct orientation
            vtkSmartPointer<vtkMatrix4x4> transformModel = vtkSmartPointer<vtkMatrix4x4>::New();
                transformModel->DeepCopy(img->getTransformMatrix());
                transformModel->Invert();

            addImageActor(img->GetImageActor(), transformModel); //!< Add data to general display

            ///Here we use the Qt signals and slots directly as it was found that the VTK-Qt connector caused problems
            ///with the image actors.
            connect(img, SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), this, SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        }
        else
            removeImageActor(img->GetImageActor()); //!< Remove data from general display

        qApp->processEvents(); ///Keep UI responsive
    }

    milxQtRenderWindow::generateRender();
}

void milxQtUnifiedWindow::generateDifference(double pseudoInfinityFactor)
{
    if(unifyModels.isEmpty())
        return;

    emit working(-1);
    double bounds[6];
    unifyModels.first()->GetOutput()->GetBounds(bounds);
    printInfo("Original Bounds: " + QString::number(bounds[0]) + ", " + QString::number(bounds[1]) + ", " + QString::number(bounds[2]) + ", " + QString::number(bounds[3]) + ", " + QString::number(bounds[4]) + ", " + QString::number(bounds[5]));

    ///Find max bounds
    double bigBounds[6];
    foreach(QPointer<milxQtModel> mdl, unifyModels)
    {
        if(mdl == unifyModels.first()) ///Compute first model to every other
            continue;

        double mdlBounds[6];
        mdl->GetOutput()->GetBounds(mdlBounds);
        printInfo("Bounds: " + QString::number(mdlBounds[0]) + ", " + QString::number(mdlBounds[1]) + ", " + QString::number(mdlBounds[2]) + ", " + QString::number(mdlBounds[3]) + ", " + QString::number(mdlBounds[4]) + ", " + QString::number(mdlBounds[5]));

        //Get largest bounds
        bigBounds[0] = milx::Minimum<double>(mdlBounds[0], bounds[0]);
        bigBounds[1] = milx::Maximum<double>(mdlBounds[1], bounds[1]);
        bigBounds[2] = milx::Minimum<double>(mdlBounds[2], bounds[2]);
        bigBounds[3] = milx::Maximum<double>(mdlBounds[3], bounds[3]);
        bigBounds[4] = milx::Minimum<double>(mdlBounds[4], bounds[4]);
        bigBounds[5] = milx::Maximum<double>(mdlBounds[5], bounds[5]);
    }

    //pad bounds by 5% more
//    bigBounds[0] -= 0.05*bigBounds[0];
//    bigBounds[1] += 0.05*bigBounds[1];
//    bigBounds[2] -= 0.05*bigBounds[2];
//    bigBounds[3] += 0.05*bigBounds[3];
//    bigBounds[4] -= 0.05*bigBounds[4];
//    bigBounds[5] += 0.05*bigBounds[5];
    printInfo("Overall Bounds: " + QString::number(bigBounds[0]) + ", " + QString::number(bigBounds[1]) + ", " + QString::number(bigBounds[2]) + ", " + QString::number(bigBounds[3]) + ", " + QString::number(bigBounds[4]) + ", " + QString::number(bigBounds[5]));

    printInfo("Computing Distance map of first model ...");
    QPointer<milxQtImage> image = new milxQtImage; //List deletion
        image->setName("Voxelisation Source");
        image->generateVoxelisedSurface(unifyModels.first()->GetOutput(), bigBounds);
        image->distanceMap();
        image->generateImage();
        emit imageAvailable(image);

    ///create larger image with label in it
    // compute dimensions of image having all surfaces
    /*int dim[3];
    for (int i = 0; i < 3; i++)
    {
        dim[i] = static_cast<int>( ceil((bigBounds[i * 2 + 1] - bigBounds[i * 2]) / image->GetCharImage()->GetSpacing()[i]) );
    }
    printInfo("Dimensions of Distance map is " + QString::number(dim[0]) + "x" + QString::number(dim[1]) + "x" + QString::number(dim[2]));

    image->resize(dim[0], dim[1], dim[2]);
    emit imageAvailable(image);*/

    /*QPointer<milxQtImage> distImage = new milxQtImage; //List deletion
      distImage->set8BitImage();
      distImage->zeros(dim[0], dim[1], dim[2], image);
      distImage->setName("Distance Image Source");
      distImage->add(image);
      distImage->distanceMap();
      distImage->generateImage();
      emit imageAvailable(image);*/

    ///Get DT
    const bool absValues = true;
    floatImageType::Pointer distanceImage = image->GetFloatImage();

    ///Voxel based computation of absolute Euclidean surface distance
    printInfo("Computing Distance map for other models ...");
    foreach(QPointer<milxQtModel> mdl, unifyModels)
    {
        if(mdl == unifyModels.first()) ///Compute first model to every other
            continue;

        printInfo(mdl->strippedBaseName() + " info:");
        printInfo("Computing Distance map ...");
        QPointer<milxQtImage> image2 = new milxQtImage; //List deletion
            image2->generateVoxelisedSurface(mdl->GetOutput(), bigBounds);
            image2->generateImage();
            image2->distanceMap();
            emit imageAvailable(image2);

        floatImageType::Pointer distanceImage2 = image2->GetFloatImage();
        printInfo("Computing Forward Distances ...");
        vtkSmartPointer<vtkFloatArray> scalars1 = surfaceScalarsFromImage(mdl->GetOutput(), distanceImage, absValues);
        printInfo("Computing Backward Distances ...");
        vtkSmartPointer<vtkFloatArray> scalars2 = surfaceScalarsFromImage(unifyModels.first()->GetOutput(), distanceImage2, absValues);

        double range1[2];
        scalars1->GetRange(range1);
        double range2[2];
        scalars2->GetRange(range2);
        double hausdorff = milx::Maximum<double>(range1[1], range2[1]);
        printInfo("Hausdorff Distance: " + QString::number(hausdorff));

        mdl->GetOutput()->GetPointData()->SetScalars(scalars1);
        mdl->GetOutput()->Modified();
    #if VTK_MAJOR_VERSION <=5
        mdl->GetOutput()->Update();
    #endif // VTK_MAJOR_VERSION
        mdl->colourMapToHOT();
        unifyModels.first()->GetOutput()->GetPointData()->SetScalars(scalars2);
        unifyModels.first()->GetOutput()->Modified();
    #if VTK_MAJOR_VERSION <=5
        unifyModels.first()->GetOutput()->Update();
    #endif // VTK_MAJOR_VERSION
    }
    unifyModels.first()->colourMapToHOT();
    emit done(-1);
}

void milxQtUnifiedWindow::generateScalarDifference()
{
    printWarning("Scalar Difference not implemented yet!");

    if(unifyModels.isEmpty())
        return;
}

void milxQtUnifiedWindow::generateCheckerBoard()
{
    if(unifyImages.isEmpty() && unifyImages.size() < 2)
    {
        printError("Not enough images being compared. Add more.");
        return;
    }

    printInfo("Generating Checkboard to View Images");

    // Create a checker pipeline
    checker = vtkSmartPointer<vtkImageCheckerboard>::New();
        checker->SetNumberOfDivisions(5,5,1);
        size_t count = 0;
        foreach(milxQtImage *img, unifyImages)
        {
            printDebug("Adding " + img->getName() + " to checkerboard");
            int *displayExtents = img->GetDisplayExtent();
            cerr   << "Slice: " << displayExtents[0] << ", " << displayExtents[1] << ", " << displayExtents[2]
                   << ", " << displayExtents[3] << ", " << displayExtents[4] << ", " << displayExtents[5];// << endl;
            vtkSmartPointer<vtkImageReslice> slice = vtkSmartPointer<vtkImageReslice>::New();
            #if VTK_MAJOR_VERSION <=5
                slice->SetInput(img->GetOutput());
            #else
                slice->SetInputData(img->GetOutput());
            #endif // VTK_MAJOR_VERSION
//                slice->SetOutputExtent(0, 9, 0, 100, 0, 0);
                slice->SetOutputExtent(displayExtents);
                slice->Update();
        #if VTK_MAJOR_VERSION <=5
            checker->SetInput(count, slice->GetOutput()); //Needs to be 2D image here
        #else
            checker->SetInputData(count, slice->GetOutput()); //Needs to be 2D image here
        #endif // VTK_MAJOR_VERSION
            slices.append(slice);
            count ++;
        }
        printDebug("Done adding");
        checker->Update();

    printDebug("Setting up checkerboard actor");
    vtkSmartPointer<vtkImageMapToWindowLevelColors> windowLevel = vtkSmartPointer<vtkImageMapToWindowLevelColors>::New();
    #if VTK_MAJOR_VERSION <=5
        windowLevel->SetInput(checker->GetOutput());
    #else
        windowLevel->SetInputData(checker->GetOutput());
    #endif

    vtkSmartPointer<vtkImageActor> checkerActor = vtkSmartPointer<vtkImageActor>::New();
//        checkerActor->GetMapper()->SetInputConnection(windowLevel->GetOutputPort());
    #if VTK_MAJOR_VERSION <=5
        checkerActor->SetInput(windowLevel->GetOutput());
    #else
        checkerActor->SetInputData(windowLevel->GetOutput());
    #endif

    // VTK widgets consist of two parts: the widget part that handles
    // event processing; and the widget representation that defines how
    // the widget appears in the scene
    // (i.e., matters pertaining to geometry).
    printDebug("Setting up checkerboard widget");
    checkerWidget = vtkSmartPointer<vtkCheckerboardWidget>::New();
        checkerWidget->SetInteractor(milxQtRenderWindow::GetRenderWindow()->GetInteractor());

    vtkCheckerboardRepresentation *checkerWidgetRep = static_cast<vtkCheckerboardRepresentation *>(checkerWidget->GetRepresentation());
        checkerWidgetRep->SetImageActor(checkerActor);
        checkerWidgetRep->SetCheckerboard(checker);

    // Add the actors to the renderer, set the background and size
    //
    milxQtRenderWindow::RemoveAllActors();
    milxQtRenderWindow::AddActor(checkerActor);
    checkerWidget->On();
    milxQtRenderWindow::refresh();
}

vtkSmartPointer<vtkFloatArray> milxQtUnifiedWindow::surfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<floatImageType> img, const bool absoluteValues)
{
  const int numberOfPoints = surface->GetNumberOfPoints();

  typedef itk::Point<double, 3> InputImagePointType;
  typedef itk::ContinuousIndex<double, 3 > ContinuousIndexType;

  typedef itk::LinearInterpolateImageFunction<floatImageType, double> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(img);

  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
  scalars->SetName("Distance");
  scalars->SetNumberOfTuples(numberOfPoints);
  scalars->SetNumberOfComponents(1);

  InputImagePointType point;
  ContinuousIndexType index;
  for(int i = 0; i < numberOfPoints; i++)
  {
    double position[3];
    surface->GetPoint(i, position);

    point[0] = position[0];
    point[1] = position[1];
    point[2] = position[2];

    img->TransformPhysicalPointToContinuousIndex(point, index);

    // If inside, mark
    if(interpolator->IsInsideBuffer(index))
    {

      double valueFound = 0.0;
      if(absoluteValues)
        valueFound = fabs(interpolator->EvaluateAtContinuousIndex(index));
      else
        valueFound = interpolator->EvaluateAtContinuousIndex(index);
      scalars->SetTuple1(i, valueFound);
    }
  }

  return scalars;
}

void milxQtUnifiedWindow::createActions()
{
    //Modes
    unionAct = new QAction(this);
    unionAct->setText(QApplication::translate("Unified", "&Union of Data", 0));
    unionAct->setShortcut(tr("Alt+u"));
    unionAct->setCheckable(true);
    unionAct->setChecked(true);
    geoDifferenceAct = new QAction(this);
    geoDifferenceAct->setText(QApplication::translate("Unified", "&Surface (Hausdoff) Distance", 0));
    geoDifferenceAct->setShortcut(tr("Alt+d"));
    geoDifferenceAct->setCheckable(true);
    scalarDifferenceAct = new QAction(this);
    scalarDifferenceAct->setText(QApplication::translate("Unified", "&Scalar Difference", 0));
    scalarDifferenceAct->setShortcut(tr("Alt+s"));
    scalarDifferenceAct->setCheckable(true);
    scalarDifferenceAct->setDisabled(true);
    checkerBoardAct = new QAction(this);
    checkerBoardAct->setText(QApplication::translate("Unified", "&Checkerboard", 0));
    checkerBoardAct->setShortcut(tr("Alt+c"));
    checkerBoardAct->setCheckable(true);
    checkerBoardAct->setDisabled(true);
    modeGroup = new QActionGroup(this); //Exclusive by default
    modeGroup->addAction(unionAct);
    modeGroup->addAction(geoDifferenceAct);
    modeGroup->addAction(scalarDifferenceAct);
    modeGroup->addAction(checkerBoardAct);

    //Refresh part of render window class
}

void milxQtUnifiedWindow::createConnections()
{
    //Modes
    connect(unionAct, SIGNAL(triggered()), this, SLOT(generateUnion()));
    connect(geoDifferenceAct, SIGNAL(triggered()), this, SLOT(generateDifference()));
    connect(scalarDifferenceAct, SIGNAL(triggered()), this, SLOT(generateScalarDifference()));
    connect(checkerBoardAct, SIGNAL(triggered()), this, SLOT(generateCheckerBoard()));

    connect(milxQtRenderWindow::refreshAct, SIGNAL(triggered()), this, SLOT(refresh()));
}

void milxQtUnifiedWindow::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection

    if(unifyImages.isEmpty())
    {
        checkerBoardAct->setDisabled(true);
    }
    else
    {
        checkerBoardAct->setDisabled(false);
    }
    if(unifyModels.isEmpty())
    {
        geoDifferenceAct->setDisabled(true);
        scalarDifferenceAct->setDisabled(true);
    }
    else
    {
        geoDifferenceAct->setDisabled(false);
        scalarDifferenceAct->setDisabled(false);
    }

    contextMenu->addAction(unionAct);
    contextMenu->addAction(geoDifferenceAct);
    contextMenu->addAction(scalarDifferenceAct);
    contextMenu->addAction(checkerBoardAct);
    contextMenu->addSeparator();
    ///Change View of Volume
    contextMenu->addMenu(milxQtRenderWindow::viewMenu);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewXY);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZX);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZY);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewAct);
    contextMenu->addMenu(milxQtRenderWindow::colourMapMenu);
    contextMenu->addAction(milxQtRenderWindow::axesAct);
    contextMenu->addAction(milxQtRenderWindow::refreshAct);

    contextMenu->exec(currentEvent->globalPos());
}

void milxQtUnifiedWindow::setCommonProperties(QWidget *passedWindow)
{
    connect(passedWindow, SIGNAL(closing(QWidget *)), this, SLOT(removeFromWindow(QWidget *)));
    //passedWindow->setDeletableOnClose(true);
}
