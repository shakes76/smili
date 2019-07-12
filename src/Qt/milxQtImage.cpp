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
#include "milxQtImage.h"

//ITK
//#include <itkImageToHistogramFilter.h>
//VTK Libraries
#include <vtkCamera.h>
#include <vtkImageMagnify.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageActor.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkInteractorStyleImage.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkImageFFT.h>
#include <vtkImageButterworthHighPass.h>
#include <vtkImageRFFT.h>
#include <vtkImageMagnitude.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageShiftScale.h>
#include <vtkImageCast.h>
#include <vtkImageFlip.h>
#include <vtkImageReslice.h>
#include <vtkImageBlend.h>
#include <vtkDistanceRepresentation.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkBoxRepresentation.h>
#include <vtkImageActorPointPlacer.h>
#include <vtkDijkstraImageContourLineInterpolator.h>
#include <vtkDijkstraImageGeodesicPath.h>
#include <vtkTextProperty.h>
#include <vtkLogLookupTable.h>

#include "milxColourMap.h"
#include "milxQtFile.h"

//milxQtImage
milxQtImage::milxQtImage(QWidget *theParent, bool contextSystem) : milxQtRenderWindow(theParent, contextSystem)
{
    ///Initialise variables
    usingVTKImage = false;
    imported = false;
    appendedData = false;
    eightbit = false;
    integer = false;
    rgb = false;
    vectorised = false;
    viewerSetup = false;
    volume = false;
    flipped = true;
    track = false;

    meanValue = 0;
    stddevValue = 0;
    minValue = 0;
    maxValue = 0;

    actualNumberOfDimensions = 3;

    ///Set strings
    milxQtWindow::prefix = "Img: ";

    ///Connect Image Object for progress updates
    linkProgressEventOf(milx::ProgressUpdates->GetUpdateObject()); //link itk observer propagator

    ///Allocate critical aspects
    imageChar = charImageType::New();
    imageInt = intImageType::New();
    imageRGB = rgbImageType::New();
    imageFloat = floatImageType::New();
    imageVector = NULL; //rare so allocate as needed

    imageData = vtkSmartPointer<vtkImageData>::New();
    viewer = vtkSmartPointer<vtkImageViewer3>::New();
    observeProgress = itkEventQtObserver::New();
    transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();

    milxQtWindow::setDeletableOnClose(true);

    //View defaults
    milxQtRenderWindow::setDefaultView(AXIAL);
    setView(AXIAL);
    setDefaultOrientation(RADIOLOGICAL);

    createActions();
    createConnections();
}

milxQtImage::~milxQtImage()
{
    ///Smart pointers handle deletion
}

void milxQtImage::setData(QPointer<milxQtImage> newImg, const bool forceDeepCopy)
{
    if(newImg->is8BitImage())
        setData(newImg->GetCharImage(), newImg->isDisplayFlipped());
    else if(newImg->is32BitImage())
        setData(newImg->GetIntImage(), newImg->isDisplayFlipped());
    else if(newImg->isVectorImage())
        setData(newImg->GetVectorImage(), newImg->isDisplayFlipped(), forceDeepCopy);
    else if(newImg->isRGBImage())
        setData(newImg->GetRGBImage(), newImg->isDisplayFlipped());
    else
        setData(newImg->GetFloatImage(), newImg->isDisplayFlipped());
}

void milxQtImage::setData(charImageType::Pointer newImg, const bool flipY)
{
    imageChar = milx::Image<charImageType>::DuplicateImage(newImg);

    flipped = flipY;

    loaded = true;
    eightbit = true;
    integer = false;
    rgb = false;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setData(intImageType::Pointer newImg, const bool flipY)
{
    imageInt = milx::Image<intImageType>::DuplicateImage(newImg);

    flipped = flipY;

    loaded = true;
    eightbit = false;
    integer = true;
    rgb = false;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setData(rgbImageType::Pointer newImg, const bool flipY)
{
    imageRGB = milx::Image<rgbImageType>::DuplicateImage(newImg);

    flipped = flipY;

    loaded = true;
    eightbit = false;
    integer = false;
    rgb = true;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setData(floatImageType::Pointer newImg, const bool flipY)
{
    imageFloat = milx::Image<floatImageType>::DuplicateImage(newImg);

    flipped = flipY;

    loaded = true;
    eightbit = false;
    integer = false;
    rgb = false;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setData(vectorImageType::Pointer newImg, const bool flipY, const bool deepCopy)
{
    if(deepCopy)
        imageVector = milx::Image<vectorImageType>::DuplicateImage(newImg);
    else
        imageVector = newImg;

    flipped = flipY;

    loaded = true;
    eightbit = false;
    integer = false;
    rgb = false;
    vectorised = true;
    usingVTKImage = false;

    magnitude();
}

void milxQtImage::setData(vnl_matrix<double> &newData)
{
	typedef double doublePixelType;
	typedef itk::Image<doublePixelType, milx::imgDimension> doubleImageType;

	doubleImageType::Pointer imageDouble = milx::Image<doubleImageType>::ImportMatrixToImage<double>(newData);
	imageFloat = milx::Image<doubleImageType>::CastImage<floatImageType>(imageDouble);
	printWarning("Matrix of Double type has been converted to Float type for display");
    loaded = true;
    usingVTKImage = false;
    integer = false;
    eightbit = false;
    rgb = false;
    vectorised = false;
}

void milxQtImage::setData(const unsigned slice, vnl_matrix<double> &newData)
{
    setData(newData);
}

void milxQtImage::setData(vtkSmartPointer<vtkImageData> newImg)
{
    imageData->DeepCopy(newImg);

    usingVTKImage = true;
    eightbit = false;
    integer = false;
    if( (newImg->GetNumberOfScalarComponents() == 4 || newImg->GetNumberOfScalarComponents() == 3) && newImg->GetScalarType() == VTK_UNSIGNED_CHAR )
        rgb = true;
    else
        rgb = false;
    loaded = true;
    vectorised = false;
    flipped = false;

    //~ histogram(256, 0, 255, false);
}

void milxQtImage::setDisplayData(QPointer<milxQtImage> newImg)
{
    if(newImg->is8BitImage())
    	setDisplayData(newImg->GetCharImage(), newImg->isDisplayFlipped());
    else if(newImg->is32BitImage())
      setDisplayData(newImg->GetIntImage(), newImg->isDisplayFlipped());
    else if(newImg->isRGBImage())
    	setDisplayData(newImg->GetRGBImage(), newImg->isDisplayFlipped());
    else if(newImg->isVectorImage())
    	setDisplayData(newImg->GetVectorImage(), newImg->isDisplayFlipped());
    else
    	setDisplayData(newImg->GetFloatImage(), newImg->isDisplayFlipped());
}

void milxQtImage::setDisplayData(charImageType::Pointer newImg, const bool flipY)
{
    imageChar = newImg;

    flipped = flipY;

    loaded = true;
    eightbit = true;
    integer = false;
    rgb = false;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setDisplayData(intImageType::Pointer newImg, const bool flipY)
{
  imageInt = newImg;

  flipped = flipY;

  loaded = true;
  eightbit = false;
  integer = true;
  rgb = false;
  vectorised = false;
  usingVTKImage = false;
}

void milxQtImage::setDisplayData(rgbImageType::Pointer newImg, const bool flipY)
{
    imageRGB = newImg;

    flipped = flipY;

    loaded = true;
    eightbit = false;
    integer = false;
    rgb = true;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setDisplayData(floatImageType::Pointer newImg, const bool flipY)
{
    imageFloat = newImg;

    flipped = flipY;

    loaded = true;
    eightbit = false;
    integer = false;
    rgb = false;
    vectorised = false;
    usingVTKImage = false;
}

void milxQtImage::setDisplayData(vectorImageType::Pointer newImg, const bool flipY)
{
    setData(newImg, flipY, false); //due to the large memory requirements not wise to duplicate
}

void milxQtImage::setSlice(int slice)
{
    if(!loaded || !volume)
        return;

    if(slice >= viewer->GetSliceMin() && slice <= viewer->GetSliceMax())
        viewer->SetSlice(slice);
    else
        printError("Slice out of range! Ignoring.");
}

void milxQtImage::SetTransform(vtkSmartPointer<vtkTransform> newTransform)
{
    if(loaded)
    {
        emit working(-1);
        ///Transform this model
        vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
        #if VTK_MAJOR_VERSION <= 5
            reslice->SetInput(imageData);
        #else
            reslice->SetInputData(imageData);
        #endif
            linkProgressEventOf(reslice);
            reslice->AutoCropOutputOn();
            reslice->TransformInputSamplingOn();
            reslice->SetOutputDimensionality(milx::imgDimension);
            reslice->SetResliceTransform(newTransform);
            reslice->Update();

        imageData = reslice->GetOutput();
        imageData->Modified();
    #if VTK_MAJOR_VERSION <= 5
        imageData->Update();
    #endif
        emit done(-1);

        usingVTKImage = true;
    }
}

void milxQtImage::generateImage(const bool quietly)
{
    if (loaded)
    {
        printDebug("Generating Image");
        int bounds[6];

        emit working(-1);
        updateData(orientAct->isChecked());

        imageData->GetExtent(bounds);

        imageInformation();

        ///Ensure small images are rendered correctly using Magnify class
        if(bounds[5] > 1)
            volume = true;

        ///Setup Viewer
    #if VTK_MAJOR_VERSION <= 5
        viewer->SetInput(imageData);
    #else
        viewer->SetInputData(imageData);
    #endif
        if(!viewerSetup)
        {
            printDebug("Setting up viewer");
            linkProgressEventOf(viewer);
            milxQtRenderWindow::SetRenderer(viewer->GetRenderer());
            printDebug("Size of Image window: " + QString::number(milxQtRenderWindow::GetRenderWindow()->GetSize()[0]) + "x" + QString::number(milxQtRenderWindow::GetRenderWindow()->GetSize()[1]));
            QVTKWidget::SetRenderWindow(viewer->GetRenderWindow());
            viewer->SetupInteractor(QVTKWidget::GetInteractor());
            SetupWidgets(viewer->GetRenderWindow()->GetInteractor());
            if(volume)
                viewer->SetSlice(bounds[5]/2); //show middle of volume

            //Remove VTK events for the right mouse button for Qt context menu
            QVTKWidget::GetInteractor()->RemoveObservers(vtkCommand::RightButtonPressEvent);
            QVTKWidget::GetInteractor()->RemoveObservers(vtkCommand::RightButtonReleaseEvent);
            Connector->Connect(QVTKWidget::GetInteractor(),
                       vtkCommand::EndWindowLevelEvent,
                       this,
                       SLOT( userEvent() ));
            Connector->Connect(QVTKWidget::GetInteractor(),
                       vtkCommand::MouseWheelForwardEvent,
                       this,
                       SLOT( userEvent() ));
            Connector->Connect(QVTKWidget::GetInteractor(),
                       vtkCommand::MouseWheelBackwardEvent,
                       this,
                       SLOT( userEvent() ));
            Connector->Connect(QVTKWidget::GetInteractor(),
                       vtkCommand::KeyPressEvent,
                       this,
                       SLOT( userEvent() ));

            viewer->GetRenderer()->ResetCamera();

            //Human Glyph setup
            setupHumanGlyph(transformMatrix);
            humanGlyph->SetDefaultRenderer(viewer->GetRenderer());
            humanGlyph->SetInteractor(viewer->GetRenderWindow()->GetInteractor());
            if(actualNumberOfDimensions > 2)
            {
                milxQtRenderWindow::humanAct->setEnabled(true);
                humanGlyph->On();
            }
            else
            {
                milxQtRenderWindow::humanAct->setEnabled(false);
                humanGlyph->Off();
            }
//            humanGlyph->InteractiveOn();

            //Sphere annotate
            sphereWidget->SetDefaultRenderer(viewer->GetRenderer());
            sphereWidget->SetInteractor(viewer->GetRenderWindow()->GetInteractor());

            printDebug("Viewer Setup Complete");
        }
        viewer->GetInteractorStyle()->InvokeEvent(vtkCommand::ResetWindowLevelEvent); //Reset window level as if pressing 'r'
        viewer->GetRenderer()->ResetCamera(); //Reset window view as if pressing 'Shift+r'
        viewer->UpdateCursor();
        viewer->Render();

        ///Check for magnification
        if(!viewerSetup)
        {
            /*viewer->GetInput()->GetExtent(magBounds);
            if(magBounds[1] > bounds[1] || magBounds[3] > bounds[3] || magBounds[5] > bounds[5])
            {
                printInfo( "Magnified Generate Size: " + QString::number(magBounds[1]+1) + "x" + QString::number(magBounds[3]+1) + "x" + QString::number(magBounds[5]+1) );
                printWarning("Image magnified by " + QString::number( static_cast<double>(magBounds[1])/bounds[1] ) + " for better display (likely because it was small).");
            }*/

            setupEvents();

            ///Resize Window
            /*if (bounds[1] > minWindowSize && bounds[3] > minWindowSize)
            {
                viewer->SetSize(bounds[1],bounds[3]);
                QVTKWidget::resize(bounds[1],bounds[3]);
            }
            else
            {
                viewer->SetSize(minWindowSize,minWindowSize);
                QVTKWidget::resize(minWindowSize,minWindowSize);
            }*/
            viewerSetup = true;
        }

        //Pass through data with no LUT initially if not RGB or RGBA image
        //According to VTK docs for vtkImageMapToWindowLevelColors: If the lookup table is not set, or is set to NULL, then the input data will be passed through if it is already of type UNSIGNED_CHAR.
        printInfo("Number of Image Components: " + QString::number(imageData->GetNumberOfScalarComponents()));
        if(imageData->GetNumberOfScalarComponents() > 2)
        {
            GetWindowLevel()->SetLookupTable(NULL);
            GetWindowLevel()->SetWindow(255);
            GetWindowLevel()->SetLevel(127.5);
        }
        lookupTable = NULL;

        if(milxQtRenderWindow::useDefaultView)
            setView(milxQtRenderWindow::defaultView); //Default view

        emit milxQtRenderWindow::modified(GetImageActor());
        if(!quietly)
            emit modified(this);
        emit done(-1);
    }
    printDebug("Completed Generating Image");
}

void milxQtImage::generateVoxelisedSurface(vtkSmartPointer<vtkPolyData> surfaceToVoxelise, double *bounds, double *spacing)
{
    bool ok1, ok2, ok3;
    double xSpacing = 1.0, ySpacing = 1.0, zSpacing = 1.0;

    if(spacing == NULL)
    {
    xSpacing = QInputDialog::getDouble(this, tr("Please enter x spacing value"),
                      tr("X Spacing Value:"), 0.5, -2147483647, 2147483647, 7, &ok1);
    ySpacing = QInputDialog::getDouble(this, tr("Please enter y spacing value"),
                      tr("Y SpacingValue:"), 0.5, -2147483647, 2147483647, 7, &ok2);
    zSpacing = QInputDialog::getDouble(this, tr("Please enter z spacing value"),
                      tr("Z Spacing Value:"), 0.5, -2147483647, 2147483647, 7, &ok3);
    }

    if(ok1 && ok2 && ok3)
    {
        emit working(-1);

        printDebug("Creating Filled Image");
        vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();

        double totalBounds[6];
        if(bounds == NULL)
        {
            bounds = &totalBounds[0];
            surfaceToVoxelise->GetBounds(bounds);
        }
        double spacing[3]; // desired volume spacing
        spacing[0] = xSpacing;
        spacing[1] = ySpacing;
        spacing[2] = zSpacing;
        whiteImage->SetSpacing(spacing);

        //Pad the bounds by a slice
        bounds[0] -= spacing[0];
        bounds[1] += spacing[0];
        bounds[2] -= spacing[1];
        bounds[3] += spacing[1];
        bounds[4] -= spacing[2];
        bounds[5] += spacing[2];

        // compute dimensions
        int dim[3];
        for (int i = 0; i < 3; i++)
        {
            dim[i] = static_cast<int>( ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]) );
        }
        printInfo("Dimensions of Voxelisation is " + QString::number(dim[0]) + "x" + QString::number(dim[1]) + "x" + QString::number(dim[2]));
        whiteImage->SetDimensions(dim);
        whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

        double origin[3];
        // NOTE: I am not sure whether or not we had to add some offset!
        origin[0] = bounds[0];// + spacing[0] / 2;
        origin[1] = bounds[2];// + spacing[1] / 2;
        origin[2] = bounds[4];// + spacing[2] / 2;
        whiteImage->SetOrigin(origin);
    #if VTK_MAJOR_VERSION <= 5
        whiteImage->SetScalarTypeToUnsignedChar();
        whiteImage->SetNumberOfScalarComponents(1);
        whiteImage->AllocateScalars();
    #else
        whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
    #endif

        // fill the image with foreground voxels:
        unsigned char inval = 255;
        unsigned char outval = 0;
        for (vtkIdType i = 0; i < whiteImage->GetNumberOfPoints(); ++i)
            whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);

        // polygonal data --> image stencil:
        printDebug("Polygonal data --> Image stencil");
        vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    #if VTK_MAJOR_VERSION <= 5
        pol2stenc->SetInput(surfaceToVoxelise);
    #else
        pol2stenc->SetInputData(surfaceToVoxelise);
    #endif
        pol2stenc->SetOutputOrigin(origin);
        pol2stenc->SetOutputSpacing(spacing);
        pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
        linkProgressEventOf(pol2stenc); //Keeps UI responsive
        pol2stenc->Update();

        // cut the corresponding white image and set the background:
        printDebug("Image stencil --> Image");
        vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
    #if VTK_MAJOR_VERSION <= 5
        imgstenc->SetInput(whiteImage);
        imgstenc->SetStencil(pol2stenc->GetOutput());
    #else
        imgstenc->SetInputData(whiteImage);
        imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
    #endif
        imgstenc->ReverseStencilOff();
        imgstenc->SetBackgroundValue(outval);
        linkProgressEventOf(imgstenc); //Keeps UI responsive
        imgstenc->Update();

        printDebug("Converting to ITK Image");
        setData( milx::Image<charImageType>::ConvertVTKImageToITKImage(imgstenc->GetOutput()) );
        usingVTKImage = false;
        setName("Voxelised Surface");
        imageInformation();
        emit done(-1);
    }
}

//VTK Filters
vtkSmartPointer<vtkImageData> milxQtImage::butterWorthHighPass(vtkSmartPointer<vtkImageData> img)
{
    printInfo("Computing the FFT of Image.");
    vtkSmartPointer<vtkImageFFT> fft = vtkSmartPointer<vtkImageFFT>::New();
//        fft->SetDimensionality(3);
#if VTK_MAJOR_VERSION <= 5
    fft->SetInput(img);
#else
    fft->SetInputData(img);
#endif
    linkProgressEventOf(fft);
    fft->Update();

    printInfo("Filtering FFT space of Image.");
    vtkSmartPointer<vtkImageButterworthHighPass> highPass = vtkSmartPointer<vtkImageButterworthHighPass>::New();
    highPass->SetInputConnection(fft->GetOutputPort());
    highPass->SetOrder(2);
    highPass->SetXCutOff(0.2);
    highPass->SetYCutOff(0.1);
    linkProgressEventOf(highPass);
    highPass->ReleaseDataFlagOff();

    printInfo("Computing the Inverse FFT of Image.");
    vtkSmartPointer<vtkImageRFFT> ifft = vtkSmartPointer<vtkImageRFFT>::New();
    ifft->SetInputConnection(highPass->GetOutputPort());
    linkProgressEventOf(ifft);
    ifft->Update();

    printInfo("Extracting the magnitude and casting to float image type.");
    vtkSmartPointer<vtkImageMagnitude> magnitudeFilter = vtkSmartPointer<vtkImageMagnitude>::New();
    magnitudeFilter->SetInputConnection(ifft->GetOutputPort());
    linkProgressEventOf(magnitudeFilter);

    vtkSmartPointer<vtkImageCast> rfftCastFilter = vtkSmartPointer<vtkImageCast>::New();
    rfftCastFilter->SetInputConnection(magnitudeFilter->GetOutputPort());
    rfftCastFilter->SetOutputScalarTypeToFloat();
    linkProgressEventOf(rfftCastFilter);
    rfftCastFilter->Update();

    return rfftCastFilter->GetOutput();
}

void milxQtImage::trackView(milxQtImage *windowToTrack, ViewType viewTo)
{
    printDebug("Tracking View");
    track = true;
    viewToTrack = viewTo;
    enableCrosshair();
    milxQtRenderWindow::Connector->Connect(windowToTrack->GetVTKInteractor(),
                                  vtkCommand::MiddleButtonPressEvent,
//                                  vtkCommand::LeftButtonPressEvent,
                                  this,
                                  SLOT( updateTrackedView(vtkObject*) ),
                                  NULL, 1.0); //High Priority
}

void milxQtImage::userEvent(QMouseEvent *event)
{
    emit milxQtRenderWindow::modified(GetImageActor());
}

void milxQtImage::userEvent(QKeyEvent *event)
{
    emit milxQtRenderWindow::modified(GetImageActor());
}

void milxQtImage::userEvent(QWheelEvent *event)
{
    emit milxQtRenderWindow::modified(GetImageActor());
}

void milxQtImage::updateCoords(vtkObject *obj)
{
    ///Get interactor
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    if(!iren)
        return;

    QString message = "";

    // Get a shortcut to the pixel data.
    vtkImageData* pImageData = viewer->GetInput();

    if(!pImageData)
        return;

    if( (pImageData->GetNumberOfScalarComponents() == 4 || pImageData->GetNumberOfScalarComponents() == 3)
        && pImageData->GetScalarType() == VTK_UNSIGNED_CHAR ) ///\todo Why is update of coordinates for RGB images become really slow?
    {
        message = "RGB Image Pixel Cooridinates Diabled for Performance";
        updateBar->showMessage(message);
        return;
    }

    ///Get event position
    ///Code initial by Mark Wyszomierski 2003-2007 @ devsample
    ///Modified and simplified by Shekhar Chandra
    ///Do the pick. It will return a non-zero value if we intersected the image.
    if (dataPicker->Pick(iren->GetEventPosition()[0],
                         iren->GetEventPosition()[1],
                         0,  // always zero.
                         renderer))
    {
        // Get the volume index within the entire volume now.
        //~ const vtkIdType nVolIdx = dataPicker->GetPointId();

        double ptMapped[3];
        dataPicker->GetMapperPosition(ptMapped);
        //~ printDebug("Mapper Position: " + QString::number(ptMapped[0]) + ", " + QString::number(ptMapped[1]) + ", " + QString::number(ptMapped[2]));

        //~ if(nVolIdx >= 0) //-1 means no point picked
        //~ {
            //~ coordinate current( pImageData->GetPoint(nVolIdx) ), cell(0.0); ///Get the pixel in real space
            coordinate current( ptMapped ), cell(0.0); ///Get the pixel in real space
            const int zAxis = viewer->GetSliceOrientation();
            int actual[3], relative[3];
            double spacing[3], origin[3];

            //~ printDebug("Prev. Current: " + QString::number(current[0]) + ", " + QString::number(current[1]) + ", " + QString::number(current[2]));
            pImageData->GetSpacing(spacing);
            pImageData->GetOrigin(origin);

            //pick is always only 2D, add ``z-axis''
            if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_XY) //z-axis is z-axis
            {
                current[2] = origin[2] + viewer->GetSlice()*spacing[2];
                //~ printDebug("Axial");
            }
            else if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_XZ) //y-axis is z-axis
            {
                //~ current[2] = current[1];
                current[1] = origin[1] + viewer->GetSlice()*spacing[1];
                //~ printDebug("Coronal");
            }
            else if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_YZ) //x-axis is z-axis
            {
                //~ current[2] = current[1];
                //~ current[1] = current[0];
                current[0] = origin[0] + viewer->GetSlice()*spacing[0];
                //~ printDebug("Sagittal");
            }

            ///Compute the actual coordinate for pixel access
            int nVolIdx = pImageData->FindPoint(current.data_block());
            const bool isInside = pImageData->ComputeStructuredCoordinates(current.data_block(), actual, cell.data_block());

            relative[0] = actual[0];
            relative[1] = actual[1];
            relative[2] = actual[2];

            //~ printDebug("Current: " + QString::number(current[0]) + ", " + QString::number(current[1]) + ", " + QString::number(current[2]));
            //~ printDebug("Actual: " + QString::number(actual[0]) + ", " + QString::number(actual[1]) + ", " + QString::number(actual[2]));
            if(isInside)
            {
                if(track)
                    enableCrosshair();

                //~ printDebug("Relative: " + QString::number(relative[0]) + ", " + QString::number(relative[1]) + ", " + QString::number(relative[2]));
                ///We have to handle different number of scalar components.
                switch (pImageData->GetNumberOfScalarComponents())
                {
                case 1:
                {
                    ///Extract first scalar component only and build message
                    const double scalar = pImageData->GetScalarComponentAsDouble(actual[0], actual[1], actual[2], 0); //Negative cause image flipped so use origin
                    message = "Point " + QString::number(nVolIdx) + ", Indices(" + QString::number(relative[0]) + ", " + QString::number(relative[1]) + ", " + QString::number(relative[2]) + "), "
                              + "Real(" + QString::number(current[0]) + ", " + QString::number(current[1]) + ", " + QString::number(current[2]) + ") = "
                              + QString::number(scalar);
                }
                break;

                case 2:
                {
                    ///Extract all three scalar components and build message
                    const double scalarR = pImageData->GetScalarComponentAsDouble(actual[0], actual[1], actual[2], 0); //Negative cause image flipped so use origin
                    const double scalarG = pImageData->GetScalarComponentAsDouble(actual[0], actual[1], actual[2], 1);
                    message = "Point " + QString::number(nVolIdx) + ", Indices(" + QString::number(relative[0]) + ", " + QString::number(relative[1]) + ", " + QString::number(relative[2]) + "), "
                              + "Real(" + QString::number(current[0]) + ", " + QString::number(current[1]) + ", " + QString::number(current[2]) + ") = "
                              + "[" + QString::number(scalarR) + "," + QString::number(scalarG) + "]";
                }
                break;

                case 3:
                {
                    ///Extract all three scalar components and build message
                    const double scalarR = pImageData->GetScalarComponentAsDouble(actual[0], actual[1], actual[2], 0); //Negative cause image flipped so use origin
                    const double scalarG = pImageData->GetScalarComponentAsDouble(actual[0], actual[1], actual[2], 1);
                    const double scalarB = pImageData->GetScalarComponentAsDouble(actual[0], actual[1], actual[2], 2);
                    message = "Point " + QString::number(nVolIdx) + ", Indices(" + QString::number(relative[0]) + ", " + QString::number(relative[1]) + ", " + QString::number(relative[2]) + "), "
                              + "Real(" + QString::number(current[0]) + ", " + QString::number(current[1]) + ", " + QString::number(current[2]) + ") = "
                              + "[" + QString::number(scalarR) + "," + QString::number(scalarG) + "," + QString::number(scalarB) + "]";
                }
                break;

                default:
                    message = "Unsupported number of scalar components";
                    break;
                }

//                if(track)
//                    emit coordinateChanged(actual[0], actual[1], actual[2]);
            }
        //~ }
    }

    ///Write message to status bar
    updateBar->showMessage(message);
}

void milxQtImage::updateSlice(vtkObject *obj)
{
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    string keyPressed = iren->GetKeySym();

    //printDebug("Updating Slice Display");
    if (keyPressed == "Up")
    {
        if (viewer->GetSlice()+1 <= viewer->GetSliceMax())
        {
            printDebug("Incrementing from Slice " + QString::number(viewer->GetSlice()) + " to Slice " + QString::number(viewer->GetSlice()+1));
            viewer->SetSlice(viewer->GetSlice()+1);
        }
    }
    else if (keyPressed == "Down")
    {
        if (viewer->GetSlice()-1 >= viewer->GetSliceMin())
        {
            printDebug("Decrementing from Slice " + QString::number(viewer->GetSlice()) + " to Slice " + QString::number(viewer->GetSlice()-1));
            viewer->SetSlice(viewer->GetSlice()-1);
        }
    }

    viewer->GetRenderWindow()->InvokeEvent(vtkCommand::ModifiedEvent, NULL);
    updateCoords(iren);
}

void milxQtImage::updateTrackedView(vtkObject *obj)
{
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);
    vtkRenderer* irenderer = iren->FindPokedRenderer(0,0);
    printDebug("Update Tracked View");

    // Get a shortcut to the pixel data.
    vtkImageData* pImageData = viewer->GetInput();

    if (dataPicker->Pick(iren->GetEventPosition()[0],
                         iren->GetEventPosition()[1],
                         0, //always zero
                         irenderer))
    {
        double ptMapped[3];
        dataPicker->GetMapperPosition(ptMapped);

        coordinate current( ptMapped ), cell(0.0); ///Get the pixel in real space
        const int zAxis = viewer->GetSliceOrientation();
        int actual[3];

//        printDebug("Prev. Current: " + QString::number(current[0]) + ", " + QString::number(current[1]) + ", " + QString::number(current[2]));

        ///Compute the actual coordinate for pixel access
        const bool isInside = pImageData->ComputeStructuredCoordinates(current.data_block(), actual, cell.data_block());

        if(isInside)
        {
            if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_XY) //z-axis is z-axis
            {
                viewer->SetSlice(actual[2]);
                printDebug("Axial");
            }
            else if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_XZ) //y-axis is z-axis
            {
                viewer->SetSlice(actual[1]);
                printDebug("Coronal");
            }
            else if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_YZ) //x-axis is z-axis
            {
                viewer->SetSlice(actual[0]);
                printDebug("Sagittal");
            }
        #if(VTK_MAJOR_VERSION > 5)
            viewer->GetCursor()->SetCenter(current[0], current[1], current[2]);
        #else
            viewer->GetCursor()->SetFocalPoint(current[0], current[1], current[2]);
        #endif
            viewer->GetCursor()->Update();
        }
    }

    viewer->GetRenderWindow()->InvokeEvent(vtkCommand::ModifiedEvent, NULL);
    viewer->Render();
    emit milxQtRenderWindow::modified(GetImageActor());
}

void milxQtImage::contour()
{
    ///Disable human Glyph
    disableOrient();

    if(!milxQtRenderWindow::contourWidget)
    {
        milxQtRenderWindow::contourWidget = vtkSmartPointer<vtkContourWidget>::New();
        milxQtRenderWindow::contourWidget->SetInteractor(QVTKWidget::GetInteractor());
        milxQtRenderWindow::contourWidget->FollowCursorOn();

        printInfo("Contour Image Mode enabled.\nLeft Click to place points, Right click to place end point.");
        printInfo("Delete key to delete point, Shift+Delete key to reset.");

        vtkSmartPointer<vtkOrientedGlyphContourRepresentation> rep = vtkOrientedGlyphContourRepresentation::SafeDownCast( milxQtRenderWindow::contourWidget->GetRepresentation() );
            rep->GetLinesProperty()->SetColor(1, 0.2, 0);
            rep->GetLinesProperty()->SetLineWidth(5.0);
//            rep->SetWorldTolerance(0.000005);

        vtkSmartPointer<vtkImageActorPointPlacer> pointPlacer = vtkSmartPointer<vtkImageActorPointPlacer>::New();
            pointPlacer->SetImageActor(viewer->GetImageActor());
        rep->SetPointPlacer(pointPlacer);

        ///Needs costs for the snapping contour
        // Gradient magnitude for edges
        vtkSmartPointer<vtkImageGradientMagnitude> grad = vtkSmartPointer<vtkImageGradientMagnitude>::New();
            grad->SetDimensionality(3);
            grad->HandleBoundariesOn();
        #if VTK_MAJOR_VERSION <= 5
            grad->SetInputConnection(imageData->GetProducerPort());
        #else
            grad->SetInputData(imageData);
        #endif
            grad->Update();

        // Invert the gradient magnitude so that low costs are
        // associated with strong edges and scale from 0 to 1
        double *range = grad->GetOutput()->GetScalarRange();
        vtkSmartPointer<vtkImageShiftScale> gradInvert = vtkSmartPointer<vtkImageShiftScale>::New();
            gradInvert->SetShift( -1.0*range[ 1 ] );
            gradInvert->SetScale( 1.0 /( range[ 0 ] - range[ 1 ] ) );
            gradInvert->SetOutputScalarTypeToFloat();
            gradInvert->SetInputConnection( grad->GetOutputPort() );
            gradInvert->Update();

        //Slice
        vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
            reslice->SetOutputExtent(viewer->GetImageActor()->GetDisplayExtent()); //get the slice info
            reslice->SetInputConnection(gradInvert->GetOutputPort());
            reslice->SetOutputSpacing(imageData->GetSpacing()); //needed
            reslice->SetOutputOrigin(imageData->GetOrigin()); //needed
//            reslice->SetOutputDimensionality(2);
            reslice->Update();

        vtkSmartPointer<vtkDijkstraImageContourLineInterpolator> interpolator = vtkSmartPointer<vtkDijkstraImageContourLineInterpolator>::New();
            interpolator->SetCostImage(reslice->GetOutput());
        vtkSmartPointer<vtkDijkstraImageGeodesicPath> path = interpolator->GetDijkstraImageGeodesicPath();
            path->StopWhenEndReachedOn();
            // prevent contour segments from overlapping
            path->RepelPathFromVerticesOn();
            // weights are scaled from 0 to 1 as are associated cost
            // components
            path->SetCurvatureWeight( 0.15 );
            path->SetEdgeLengthWeight( 0.8 );
            path->SetImageWeight( 1.0 );
        rep->SetLineInterpolator(interpolator);
    }

    if(milxQtRenderWindow::contourAct->isChecked())
    {
        milxQtRenderWindow::contourWidget->EnabledOn();
        printInfo("Enabled Contouring");
    }
    else
    {
        milxQtRenderWindow::contourWidget->EnabledOff();
        printInfo("Disabled Contouring");
    }

    refresh();
}

void milxQtImage::updateData(const bool orient)
{
    if(!usingVTKImage)
    {
        printDebug("Updating Data");

        floatImageType::DirectionType direction;
        floatImageType::PointType origin;
        floatImageType::SpacingType spacing;
        vtkSmartPointer<vtkImageData> newImageData = vtkSmartPointer<vtkImageData>::New();
        /*if(eightbit)
        {
            /// ITK to VTK image (unsigned char)
            if(orient)
                imageChar = milx::Image<charImageType>::ApplyOrientationToITKImage<charImageType, float>(imageChar, imageChar, true, flipped);
            imageData->DeepCopy(milx::Image<charImageType>::ConvertITKImageToVTKImage(imageChar));
            direction = imageChar->GetDirection();
            origin = imageChar->GetOrigin();
            spacing = imageChar->GetSpacing();
        }
        else if(rgb)
        {
            /// ITK to VTK image (RGB)
            if(orient)
                imageRGB = milx::Image<rgbImageType>::ApplyOrientationToITKImage<rgbImageType, float>(imageRGB, imageRGB, true, flipped);
            imageData->DeepCopy(milx::Image<rgbImageType>::ConvertITKImageToVTKImage(imageRGB));
            direction = imageRGB->GetDirection();
            origin = imageRGB->GetOrigin();
            spacing = imageRGB->GetSpacing();
        }
        else //if float and/or vector (which also generates float magnitude image)
        {
            /// ITK to VTK image (Float)
            if(orient)
                imageFloat = milx::Image<floatImageType>::ApplyOrientationToITKImage<floatImageType, float>(imageFloat, imageFloat, true, flipped);
            imageData->DeepCopy(milx::Image<floatImageType>::ConvertITKImageToVTKImage(imageFloat));
            direction = imageFloat->GetDirection();
            origin = imageFloat->GetOrigin();
            spacing = imageFloat->GetSpacing();
        }*/
        if(eightbit)
        {
            /// ITK to VTK image (unsigned char)
            newImageData->DeepCopy( milx::Image<charImageType>::ConvertITKImageToVTKImage(imageChar) );
            direction = imageChar->GetDirection();
            origin = imageChar->GetOrigin();
            spacing = imageChar->GetSpacing();
            if(orient)
                imageData = milx::Image<charImageType>::ApplyOrientationToVTKImage(newImageData, imageChar, transformMatrix, true, flipped);
            else
                imageData = newImageData;
            //Labelled image flag is set as true to avoid artefacts in resampling within the ApplyOrientationToVTKImage member
            printDebug("Updated Internal Char Image Data");
        }
        else if(integer)
        {
          /// ITK to VTK image (RGB)
          newImageData->DeepCopy(milx::Image<intImageType>::ConvertITKImageToVTKImage(imageInt));
          direction = imageInt->GetDirection();
          origin = imageInt->GetOrigin();
          spacing = imageInt->GetSpacing();
          if(orient)
            imageData = milx::Image<intImageType>::ApplyOrientationToVTKImage(newImageData, imageInt, transformMatrix, true, flipped);
          else
            imageData = newImageData;
          //Labelled image flag is set as true to avoid artefacts in resampling within the ApplyOrientationToVTKImage member
          printDebug("Updated Internal Integer Image Data");
        }
        else if(rgb)
        {
            /// ITK to VTK image (RGB)
            newImageData->DeepCopy( milx::Image<rgbImageType>::ConvertITKImageToVTKImage(imageRGB) );
            direction = imageRGB->GetDirection();
            origin = imageRGB->GetOrigin();
            spacing = imageRGB->GetSpacing();
            if(orient)
                imageData = milx::Image<rgbImageType>::ApplyOrientationToVTKImage(newImageData, imageRGB, transformMatrix, true, flipped);
            else
                imageData = newImageData;
            //Labelled image flag is set as true to avoid artefacts in resampling within the ApplyOrientationToVTKImage member
            printDebug("Updated Internal RGB Image Data");
        }
        else //if float and/or vector (which also generates float magnitude image)
        {
            /// ITK to VTK image (Float)
            newImageData->DeepCopy( milx::Image<floatImageType>::ConvertITKImageToVTKImage(imageFloat) );
            direction = imageFloat->GetDirection();
            origin = imageFloat->GetOrigin();
            spacing = imageFloat->GetSpacing();
            if(orient)
                imageData = milx::Image<floatImageType>::ApplyOrientationToVTKImage(newImageData, imageFloat, transformMatrix, false, flipped);
            else
                imageData = newImageData;
            //Labelled image flag is set as true to avoid artefacts in resampling within the ApplyOrientationToVTKImage member
            printDebug("Updated Internal Float Image Data");
        }
    }

    imageData->Modified();
#if VTK_MAJOR_VERSION <= 5
    imageData->Update();
#endif

    printDebug("Updated Image Data as " + QString(imageData->GetScalarTypeAsString()));
}

void milxQtImage::setupEvents()
{
    //Do not move, needs to be connected after setting up viewer!
    milxQtRenderWindow::Connector->Connect(milxQtRenderWindow::GetInteractor(),
                                  vtkCommand::KeyPressEvent,
                                  this,
                                  SLOT( updateSlice(vtkObject *) ),
                                  NULL, 1.0); //High Priority
}

void milxQtImage::autoLevel(float percentile)
{
    printInfo("Auto Updating Window Level");
    printDebug("Current Window:" + QString::number(GetIntensityWindow()));
    printDebug("Current Level:" + QString::number(GetIntensityLevel()));

    int bins = 256;
    float belowValue = -1000, aboveValue = 4000, lowerPercentile = 1.0-percentile, upperPercentile = percentile;
    if(maxValue == minValue)
        histogram(bins, belowValue, aboveValue, false); //above and below unused here, uses image min/max automatically
    belowValue = minValue;
    aboveValue = maxValue;

    emit working(-1);
    //Compute the percentile contributions to trim levels for better contrast
    size_t k = 0, kMin = 0, kMax = 0, kMid = 0;
    double binSpacing = (aboveValue-belowValue)/bins, accummulator = 0.0, proportion = 0.0;
    //Upper limit of histgram
    k = 0;
    accummulator = 0.0;
    for(double j = belowValue; j < aboveValue; j += binSpacing, k ++)
      {
        proportion = accummulator/hist->GetVoxelCount();
        if(proportion >= lowerPercentile)
          break;
        accummulator += hist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(k); //freq k
      }
    kMin = k-1;
//    printDebug("k low:" + QString::number(kMin));
//    printDebug("accummulator low:" + QString::number(accummulator));
    double lowLevel = minValue+(kMin)*binSpacing;
    //Lower limit of histgram
    k = bins-1;
    accummulator = 0.0;
    for(double j = aboveValue; j > belowValue; j -= binSpacing, k --)
    {
        proportion = accummulator/hist->GetVoxelCount();
        if(proportion >= lowerPercentile)
            break;
        accummulator += hist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(k); //freq k
    }
    kMax = k+1;
//    printDebug("k high:" + QString::number(kMax));
//    printDebug("accummulator high:" + QString::number(accummulator));
    double maxLevel = minValue+(kMax)*binSpacing;
    double windowLevel = maxLevel-lowLevel;
    double level = lowLevel+windowLevel/2; //center the new window from low limit
    emit done(-1);

    printDebug("Histogram Low Level:" + QString::number(lowLevel));
    printDebug("Histogram High Level:" + QString::number(maxLevel));
    printInfo("Window:" + QString::number(windowLevel));
    printInfo("Level:" + QString::number(level));

    viewer->SetColorWindow(windowLevel);
    viewer->SetColorLevel(level);
    viewer->Render();
}

void milxQtImage::setLevel(int level)
{
//    printInfo("Updating Window Level to "+QString::number(level)+"%");

    double range[2];
    imageData->GetScalarRange(range);
    double belowValue = range[0];
    double aboveValue = range[1];
    if(imageData->GetNumberOfScalarComponents() > 2) //Use RGB range
    {
        belowValue = 0;
        aboveValue = 255;
    }

    float lvl = (aboveValue-belowValue)*level/100 + belowValue;
    printDebug("Contrast Level is "+QString::number(lvl));

    viewer->SetColorWindow(aboveValue-lvl);
    if(imageData->GetNumberOfScalarComponents() > 2) //Use RGB range
        viewer->SetColorWindow(255);
    viewer->SetColorLevel(lvl);

    refresh();
}

double milxQtImage::GetIntensityWindow()
{
    return viewer->GetWindowLevel()->GetWindow();
}

double milxQtImage::GetIntensityLevel()
{
    return viewer->GetWindowLevel()->GetLevel();
}

void milxQtImage::SetIntensityWindow(double window)
{
    viewer->GetWindowLevel()->SetWindow(window);
}

void milxQtImage::SetIntensityLevel(double level)
{
    viewer->GetWindowLevel()->SetLevel(level);
}

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
void milxQtImage::overlay(QString filename)
{
    printInfo("Overlaying Labelled Image");

    if(filename.isEmpty())
        filename = getOpenFilename("Select Labelled Image");

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> labelledImage = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, labelledImage);

    if(!success)
        return;

//    binaryImage->contour(); //!< generate contour from binary image

    printDebug("Overlaying on Image");
    emit working(-1);
    if(eightbit)
    {
        if(labelledImage->is8BitImage())
            imageRGB = milx::Image<charImageType>::Overlay<charImageType>(imageChar, labelledImage->GetCharImage());
        else if(labelledImage->is32BitImage())
            imageRGB = milx::Image<charImageType>::Overlay<intImageType>(imageChar, labelledImage->GetIntImage());
        else
            imageRGB = milx::Image<charImageType>::Overlay<floatImageType>(imageChar, labelledImage->GetFloatImage());
    }
    else if(integer)
    {
        if(labelledImage->is8BitImage())
            imageRGB = milx::Image<intImageType>::Overlay<charImageType>(imageInt, labelledImage->GetCharImage());
        else if(labelledImage->is32BitImage())
            imageRGB = milx::Image<intImageType>::Overlay<intImageType>(imageInt, labelledImage->GetIntImage());
        else
            imageRGB = milx::Image<intImageType>::Overlay<floatImageType>(imageInt, labelledImage->GetFloatImage());
    }
//    else if(rgb)
//        imageRGB = milx::Image::Overlay<rgbImageType>(imageRGB, labelledImage->GetRGBImage());
    else
    {
        if(labelledImage->is8BitImage())
            imageRGB = milx::Image<floatImageType>::Overlay<charImageType>(imageFloat, labelledImage->GetCharImage());
        else if(labelledImage->is32BitImage())
            imageRGB = milx::Image<floatImageType>::Overlay<intImageType>(imageFloat, labelledImage->GetIntImage());
        else
            imageRGB = milx::Image<floatImageType>::Overlay<floatImageType>(imageFloat, labelledImage->GetFloatImage());
    }
    printWarning("Images Cast to RGB Image for colour!");
    printWarning("Origin may be reset!");
    emit done(-1);

    eightbit = false;
    integer = false;
    rgb = true;
    generateImage();
}

void milxQtImage::overlayContour(QString filename)
{
    printInfo("Overlaying Labelled Image as Contour");

    if(filename.isEmpty())
    {
        filename = getOpenFilename("Select Labelled Image");
    }

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> labelledImage = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, labelledImage);

    if(!success)
        return;

//    binaryImage->contour(); //!< generate contour from binary image

    printDebug("Overlaying Contour on Image");
    emit working(-1);
    if(eightbit)
    {
        if(labelledImage->is8BitImage())
            imageRGB = milx::Image<charImageType>::OverlayContour<charImageType>(imageChar, labelledImage->GetCharImage());
        else if(labelledImage->is32BitImage())
            imageRGB = milx::Image<charImageType>::OverlayContour<intImageType>(imageChar, labelledImage->GetIntImage());
        else
            imageRGB = milx::Image<charImageType>::OverlayContour<floatImageType>(imageChar, labelledImage->GetFloatImage());
    }
    else if(integer)
    {
        if(labelledImage->is8BitImage())
            imageRGB = milx::Image<intImageType>::OverlayContour<charImageType>(imageInt, labelledImage->GetCharImage());
        else if(labelledImage->is32BitImage())
            imageRGB = milx::Image<intImageType>::OverlayContour<intImageType>(imageInt, labelledImage->GetIntImage());
        else
            imageRGB = milx::Image<intImageType>::OverlayContour<floatImageType>(imageInt, labelledImage->GetFloatImage());
    }
//    else if(rgb)
//        imageRGB = milx::Image::OverlayContour<rgbImageType>(imageRGB, labelledImage->GetRGBImage());
    else
    {
        if(labelledImage->is8BitImage())
            imageRGB = milx::Image<floatImageType>::OverlayContour<charImageType>(imageFloat, labelledImage->GetCharImage());
        else if(labelledImage->is32BitImage())
            imageRGB = milx::Image<floatImageType>::OverlayContour<intImageType>(imageFloat, labelledImage->GetIntImage());
        else
            imageRGB = milx::Image<floatImageType>::OverlayContour<floatImageType>(imageFloat, labelledImage->GetFloatImage());
    }
    printWarning("Images Cast to RGB Image for colour!");
    printWarning("Origin may be reset!");
    emit done(-1);

    eightbit = false;
    integer = false;
    rgb = true;
    generateImage();
}

void milxQtImage::computeContour()
{
    if(usingVTKImage)
    {
        printError("Contour from VTK image not support yet.");
        return;
    }

    printInfo("Computing Contour of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::BinaryContour(imageChar, minValue, maxValue);
    else if(integer)
        imageInt = milx::Image<intImageType>::BinaryContour(imageInt, minValue, maxValue);
    else if(rgb)
        imageRGB = milx::Image<rgbImageType>::BinaryContour(imageRGB, minValue, maxValue);
    else
        imageFloat = milx::Image<floatImageType>::BinaryContour(imageFloat, minValue, maxValue);
    emit done(-1);

    generateImage();
}
#endif

void milxQtImage::blend(QString filename, float opacity)
{
    if(usingVTKImage)
    {
        printError("Information from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    printInfo("Displaying Blending of the Images");
//    emit working(-1);

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
    {
        emit done(-1);
        return;
    }

    imageToMatch->setName(filename);
//    imageToMatch->disableDefaultView();
    imageToMatch->generateImage();
    imageToMatch->setView(getView());

    blend(imageToMatch, opacity);
}

void milxQtImage::blend(milxQtImage *imageToMatch, float opacity)
{
    bool ok1 = true;
    if(opacity < 0.0)
    {
        opacity = QInputDialog::getDouble(this, tr("Please Provide the opacity of the blending"),
                         tr("Opacity:"), 0.5, 0, 2147483647, 2, &ok1);
    }

    if(!ok1)
        return;

    if(!GetLookupTable() || !imageToMatch->GetLookupTable())
        printWarning("Lookuptable for one of the images was not set.");

    emit working(-1);
    vtkSmartPointer<vtkImageData> ucharData1 = GetWindowLevel()->GetOutput();

    printInfo("Blending with Image: " + imageToMatch->strippedName());
    vtkSmartPointer<vtkImageData> ucharData2 = imageToMatch->GetWindowLevel()->GetOutput();
    printInfo("Number of components in Blending Image: " + QString::number(ucharData2->GetNumberOfScalarComponents()));

    // Combine the images (blend takes multiple connections on the 0th input port)
    vtkSmartPointer<vtkImageBlend> blend = vtkSmartPointer<vtkImageBlend>::New();
    #if VTK_MAJOR_VERSION <= 5
        blend->AddInput(ucharData1);
        blend->AddInput(ucharData2);
    #else
        blend->AddInputData(ucharData1);
        blend->AddInputData(ucharData2);
    #endif
        blend->SetOpacity(0,opacity);
        blend->SetOpacity(1,opacity);
        linkProgressEventOf(blend);
//        blend->SetBlendModeToCompound();
        blend->SetBlendModeToNormal();
        blend->Update();

    printInfo("Number of components in Blended Image: " + QString::number(blend->GetOutput()->GetNumberOfScalarComponents()));
    printInfo("Type of pixel in Blended Image: " + QString(blend->GetOutput()->GetScalarTypeAsString()));

    setData(blend->GetOutput());
    emit done(-1);

    //Force RGB colouring
    GetWindowLevel()->SetLookupTable(NULL);
    lookupTable = NULL;
    GetWindowLevel()->SetWindow(255);
    GetWindowLevel()->SetLevel(127.5);
    viewer->GetRenderer()->ResetCamera();

    generateImage();
}

void milxQtImage::volumeRendering()
{
    printInfo("Displaying Volume Rendering of Image");

    vtkSmartPointer<vtkImageData> img;
    if(flipped)
    {
        vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
        #if VTK_MAJOR_VERSION <= 5
            imageReorient->SetInput(imageData);
        #else
            imageReorient->SetInputData(imageData);
        #endif
            imageReorient->SetFilteredAxis(1);
            imageReorient->FlipAboutOriginOn();
            linkProgressEventOf(imageReorient);
            imageReorient->Update(); //ITK image would have been flipped
            img = imageReorient->GetOutput();
    }
    else
        img = imageData;

    ///Pass on the data so that plot class can be used without including plot class.
    ///Volume rendering is generated in main class or caller.
    emit imageToVolume(img, eightbit);
}

void milxQtImage::imageInformation()
{
    printInfo("Computing Information of Image");
    emit working(-1);
    //~ histogram(256, 0, 255, false);

    floatImageType::DirectionType direction;
    floatImageType::PointType origin;
    floatImageType::SpacingType spacing;
    floatImageType::RegionType::SizeType imageSize;
    if(usingVTKImage)
    {
        printWarning("Information from VTK image, as data not loaded as ITK image.");
        int extents[6];

        imageData->GetOrigin(origin.GetDataPointer());
        imageData->GetSpacing(spacing.GetDataPointer());
        imageData->GetExtent(extents);

        imageSize[0] = extents[1]-extents[0];
        imageSize[1] = extents[3]-extents[2];
        imageSize[2] = extents[5]-extents[4];
    }
    else
    {
        if(eightbit)
        {
            origin = imageChar->GetOrigin();
            spacing = imageChar->GetSpacing();
            imageSize = imageChar->GetLargestPossibleRegion().GetSize();
            direction = imageChar->GetDirection();
            printInfo("Image loaded as 8-bit image.");
        }
        else if(integer)
        {
          origin = imageInt->GetOrigin();
          spacing = imageInt->GetSpacing();
          imageSize = imageInt->GetLargestPossibleRegion().GetSize();
          direction = imageInt->GetDirection();
          printInfo("Image loaded as Integer image.");
        }
        else if(rgb)
        {
            origin = imageRGB->GetOrigin();
            spacing = imageRGB->GetSpacing();
            imageSize = imageRGB->GetLargestPossibleRegion().GetSize();
            direction = imageRGB->GetDirection();
            printInfo("Image loaded as RGB image.");
        }
        else
        {
            origin = imageFloat->GetOrigin();
            spacing = imageFloat->GetSpacing();
            imageSize = imageFloat->GetLargestPossibleRegion().GetSize();
            direction = imageFloat->GetDirection();
            printInfo("Image loaded as floating-point image.");
        }
    }

    printInfo("Origin: (" + QString::number(origin[0]) + ", " + QString::number(origin[1]) + ", " + QString::number(origin[2]) + ")");
    printInfo("Spacing: (" + QString::number(spacing[0]) + ", " + QString::number(spacing[1]) + ", " + QString::number(spacing[2]) + ")");
    printInfo("Size: (" + QString::number(imageSize[0]) + ", " + QString::number(imageSize[1]) + ", " + QString::number(imageSize[2]) + ")");
    printInfo("Real Size: (" + QString::number(spacing[0]*imageSize[0]) + ", " + QString::number(spacing[1]*imageSize[1]) + ", " + QString::number(spacing[2]*imageSize[2]) + ")");
    if(!usingVTKImage)
    {
        printInfo("Direction/Orientation: ");
        printInfo("|" + QString::number(direction(0,0)) + ", " + QString::number(direction(0,1)) + ", " + QString::number(direction(0,2)) + "|");
        printInfo("|" + QString::number(direction(1,0)) + ", " + QString::number(direction(1,1)) + ", " + QString::number(direction(1,2)) + "|");
        printInfo("|" + QString::number(direction(2,0)) + ", " + QString::number(direction(2,1)) + ", " + QString::number(direction(2,2)) + "|");
        
        QString orientFlagStr;
        if(eightbit)
          orientFlagStr = milx::Image<charImageType>::ImageOrientation(imageChar).c_str();
        else if(integer)
          orientFlagStr = milx::Image<intImageType>::ImageOrientation(imageInt).c_str();
        else if(rgb)
          orientFlagStr = milx::Image<rgbImageType>::ImageOrientation(imageRGB).c_str();
        else
          orientFlagStr = milx::Image<floatImageType>::ImageOrientation(imageFloat).c_str();
        printInfo("Orientation Flag: " + orientFlagStr);
    }

    emit done(-1);
}

void milxQtImage::rescale()
{
    if(usingVTKImage)
    {
        printError("Rescale from VTK image not support yet.");
        return;
    }

    bool ok1, ok2;
    float newMinValue = QInputDialog::getDouble(this, tr("Please Provide the minimum value of new intensities"),
                     tr("Minimum:"), minValue, -2147483647, 2147483647, 1, &ok1);
    float newMaxValue = QInputDialog::getDouble(this, tr("Please Provide the maximum value of new intensities"),
                     tr("Maximum:"), maxValue, -2147483647, 2147483647, 1, &ok2);

    if(!ok1 || !ok2)
        return;

    printInfo("Rescaling Intensities of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::RescaleIntensity(imageChar, newMinValue, newMaxValue);
    else if(integer)
        imageInt = milx::Image<intImageType>::RescaleIntensity(imageInt, newMinValue, newMaxValue);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::RescaledIntensity(imageRGB, newMinValue, newMaxValue);
    else
        imageFloat = milx::Image<floatImageType>::RescaleIntensity(imageFloat, newMinValue, newMaxValue);
    emit done(-1);

    generateImage();
}

void milxQtImage::relabel()
{
    if(usingVTKImage)
    {
        printError("Relabelling from VTK image not support yet.");
        return;
    }

    printInfo("Relabelling Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::RelabelImage(imageChar);
    else
    {
        printError("Relabelling only supports labelled images.");
        return;
    }
    emit done(-1);

    generateImage();
}

void milxQtImage::histogramEqualisation()
{
    if(usingVTKImage)
    {
        printError("Histogram equalisation from VTK image not support yet.");
        return;
    }

    bool ok1, ok2, ok3;
    float alphaValue = QInputDialog::getDouble(this, tr("Please Provide the alpha value (towards 1 produces unsharp mask)"),
                     tr("Alpha:"), 0.3, 0, 1.0, 5, &ok1);
    float betaValue = QInputDialog::getDouble(this, tr("Please Provide the beta value (towards 0 produces unsharp mask)"),
                     tr("Beta:"), 0.3, 0, 1.0, 5, &ok2);
    int radiusValue = QInputDialog::getInteger(this, tr("Please Provide the radius value (smaller for fine detail)"),
                     tr("Radius:"), 5, 1, 2147483647, 1, &ok3);

    if(!ok1 || !ok2 || !ok3)
        return;

    printInfo("Histogram Equalisation of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::HistogramEqualisation(imageChar, alphaValue, betaValue);
    else if(integer)
        imageInt = milx::Image<intImageType>::HistogramEqualisation(imageInt, alphaValue, betaValue);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::HistogramEqualisation(imageRGB);
    else
        imageFloat = milx::Image<floatImageType>::HistogramEqualisation(imageFloat, alphaValue, betaValue);
    emit done(-1);

    generateImage();
}

void milxQtImage::gradientMagnitude()
{
    if(usingVTKImage)
    {
        printError("Gradient Magnitude from VTK image not support yet.");
        return;
    }

    printInfo("Computing Gradient Magnitude of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::GradientMagnitude(imageChar);
    else if(integer)
        imageInt = milx::Image<intImageType>::GradientMagnitude(imageInt);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::GradientMagnitude(imageRGB);
    else
        imageFloat = milx::Image<floatImageType>::GradientMagnitude(imageFloat);
    emit done(-1);

    generateImage();
}

void milxQtImage::sobelEdges()
{
    if(usingVTKImage)
    {
        printError("Sobel Edges from VTK image not support yet.");
        return;
    }

    printInfo("Computing Sobel Edges of Image");
    emit working(-1);
    if(eightbit)
    {
        imageFloat = milx::Image<floatImageType>::SobelEdges( milx::Image<charImageType>::CastImage<floatImageType>(imageChar) );
        eightbit = false;
    }
    else if(integer)
    {
        imageFloat = milx::Image<floatImageType>::SobelEdges( milx::Image<intImageType>::CastImage<floatImageType>(imageInt) );
        integer = false;
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::SobelEdges(imageRGB);
    else
        imageFloat = milx::Image<floatImageType>::SobelEdges(imageFloat);
    emit done(-1);

    generateImage();
}

void milxQtImage::cannyEdges()
{
    if(usingVTKImage)
    {
        printError("Canny Edges from VTK image not support yet.");
        return;
    }

    bool ok1, ok2, ok3;
    float variance = QInputDialog::getDouble(this, tr("Please Provide the variance"),
                     tr("Variance:"), 2.0, 0.0, 1000.0, 5, &ok1);
    float upper = QInputDialog::getDouble(this, tr("Please Provide the upper threshold"),
                                          tr("Upper Threshold:"), maxValue, minValue, maxValue, 5, &ok2);
    float lower = QInputDialog::getDouble(this, tr("Please Provide the lower threshold"),
                                          tr("Lower Threshold:"), minValue, minValue, maxValue, 5, &ok3);

    if(!ok1 || !ok2 || !ok3)
        return;

    printInfo("Computing Canny Edges of Image");
    emit working(-1);
    if(eightbit)
    {
        //must be float type image
        imageFloat = milx::Image<floatImageType>::CannyEdges( milx::Image<charImageType>::CastImage<floatImageType>(imageChar), variance, lower, upper );
        eightbit = false;
    }
    else if(integer)
    {
        //must be float type image
        imageFloat = milx::Image<floatImageType>::CannyEdges(milx::Image<intImageType>::CastImage<floatImageType>(imageInt), variance, lower, upper);
        integer = false;
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::CannyEdges(imageRGB, variance, lower, upper);
    else
        imageFloat = milx::Image<floatImageType>::CannyEdges(imageFloat, variance, lower, upper);
    emit done(-1);

    generateImage();
}

void milxQtImage::laplacian()
{
    if(usingVTKImage)
    {
        printError("Laplacian from VTK image not support yet.");
        return;
    }

    printInfo("Computing Laplacian of Image");
    emit working(-1);
    if(eightbit)
    {
        //must be float type image
        imageFloat = milx::Image<floatImageType>::Laplacian( milx::Image<charImageType>::CastImage<floatImageType>(imageChar) );
        eightbit = false;
    }
    else if(integer)
    {
        //must be float type image
        imageFloat = milx::Image<floatImageType>::Laplacian(milx::Image<intImageType>::CastImage<floatImageType>(imageInt));
        integer = false;
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::Laplacian(imageRGB);
    else
        imageFloat = milx::Image<floatImageType>::Laplacian(imageFloat);
    emit done(-1);

    generateImage();
}

void milxQtImage::normalize()
{
    if(usingVTKImage)
    {
        printError("Normalizing from VTK image not support yet.");
        return;
    }

    printInfo("Normalizing Image");
    emit working(-1);
    if(eightbit)
    {
        imageFloat = milx::Image<charImageType>::Normalization(imageChar);
        eightbit = false;
    }
    else if(integer)
    {
        imageFloat = milx::Image<intImageType>::Normalization(imageInt);
        integer = false;
    }
//    else if(rgb)
//    {
//        imageFloat = milx::Image<rgbImageType>::Normalization(imageRGB);
//        rgb = false;
//    }
    else
        imageFloat = milx::Image<floatImageType>::Normalization(imageFloat);
    emit done(-1);

    generateImage();
}

void milxQtImage::invertIntensity()
{
    if(usingVTKImage)
    {
        printError("Invert Intensity from VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false); //needed to determine max value

    printInfo("Inverting Intensity of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::InvertIntensity(imageChar, maxValue);
    else if(integer)
        imageInt = milx::Image<intImageType>::InvertIntensity(imageInt, maxValue);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::InvertIntensity(imageRGB, maxValue);
    else
        imageFloat = milx::Image<floatImageType>::InvertIntensity(imageFloat, maxValue);
    emit done(-1);

    generateImage();
}

void milxQtImage::matchInfo(milxQtImage *imageToMatch)
{
    printInfo("Matching Info of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::MatchInformation(imageChar, imageToMatch->GetCharImage());
    else if(integer)
        imageInt = milx::Image<intImageType>::MatchInformation(imageInt, imageToMatch->GetIntImage());
    else if(rgb)
        imageRGB = milx::Image<rgbImageType>::MatchInformation(imageRGB, imageToMatch->GetRGBImage());
    else
        imageFloat = milx::Image<floatImageType>::MatchInformation(imageFloat, imageToMatch->GetFloatImage());
    emit done(-1);

    generateImage();
}

void milxQtImage::matchInfo(QString filename)
{
    if(usingVTKImage)
    {
        printError("Matching Info from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

    matchInfo(imageToMatch);
}

void milxQtImage::matchHistogram(milxQtImage *imageToMatch)
{
    printInfo("Matching Histogram of Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::MatchHistogram(imageChar, imageToMatch->GetCharImage());
    else if(integer)
        imageInt = milx::Image<intImageType>::MatchHistogram(imageInt, imageToMatch->GetIntImage());
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::MatchInformation(imageRGB, imageToMatch->GetRGBImage());
    else
        imageFloat = milx::Image<floatImageType>::MatchHistogram(imageFloat, imageToMatch->GetFloatImage());
    emit done(-1);

    generateImage();
}

void milxQtImage::matchHistogram(QString filename)
{
    if(usingVTKImage)
    {
        printError("Matching Histogram from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

    matchHistogram(imageToMatch);
}

void milxQtImage::resample(QString filename)
{
    if(usingVTKImage)
    {
        printError("Resampling from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

    printInfo("Resampling Image");
    emit working(-1);
    if(eightbit)
    {
        if(imageToMatch->is8BitImage())
            imageChar = milx::Image<charImageType>::ResampleImage<charImageType>(imageChar, imageToMatch->GetCharImage());
        else if(imageToMatch->is32BitImage())
        {
            imageInt = milx::Image<charImageType>::ResampleImage<intImageType>(imageChar, imageToMatch->GetIntImage());
            eightbit = false;
            integer = true;
        }
        else
        {
            imageFloat = milx::Image<charImageType>::ResampleImage<floatImageType>(imageChar, imageToMatch->GetFloatImage());
            eightbit = false;
        }
    }
    else if(integer)
    {
        if(imageToMatch->is8BitImage())
        {
            imageChar = milx::Image<intImageType>::ResampleImage<charImageType>(imageInt, imageToMatch->GetCharImage());
            eightbit = true;
            integer = false;
        }
        else if(imageToMatch->is32BitImage())
            imageInt = milx::Image<intImageType>::ResampleImage<intImageType>(imageInt, imageToMatch->GetIntImage());
        else
        {
            imageFloat = milx::Image<intImageType>::ResampleImage<floatImageType>(imageInt, imageToMatch->GetFloatImage());
            integer = false;
        }
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::ResampleImage(imageRGB, imageToMatch->GetRGBImage());
    else
    {
        if(imageToMatch->is8BitImage())
        {
            imageChar = milx::Image<floatImageType>::ResampleImage<charImageType>(imageFloat, imageToMatch->GetCharImage());
            eightbit = true;
        }
        else
            imageFloat = milx::Image<floatImageType>::ResampleImage<floatImageType>(imageFloat, imageToMatch->GetFloatImage());
    }
    emit done(-1);

    generateImage();
}

void milxQtImage::mask(QString filename)
{
    if(usingVTKImage)
    {
        printError("Masking from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename("Select Mask");

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

    printInfo("Masking Image");
    emit working(-1);
    if(eightbit)
    {
        if(imageToMatch->is8BitImage())
            imageChar = milx::Image<charImageType>::MaskImage<charImageType>(imageChar, imageToMatch->GetCharImage());
        else
            imageChar = milx::Image<charImageType>::MaskImage<floatImageType>(imageChar, imageToMatch->GetFloatImage());
    }
    else if(integer)
    {
        if(imageToMatch->is8BitImage())
            imageInt = milx::Image<intImageType>::MaskImage<charImageType>(imageInt, imageToMatch->GetCharImage());
        else
            imageInt = milx::Image<intImageType>::MaskImage<intImageType>(imageInt, imageToMatch->GetIntImage());
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::MaskImage<rgbImageType>(imageRGB, imageToMatch->GetRGBImage());
#if ITK_VERSION_MAJOR > 3
    else if(vectorised)
    {
        if(imageToMatch->is8BitImage())
            imageVector = milx::Image<vectorImageType>::MaskImage<charImageType>(imageVector, imageToMatch->GetCharImage());
        else
            imageVector = milx::Image<vectorImageType>::MaskImage<floatImageType>(imageVector, imageToMatch->GetFloatImage());

        magnitude();
    }
#endif
    else
    {
        if(imageToMatch->is8BitImage())
            imageFloat = milx::Image<floatImageType>::MaskImage<charImageType>(imageFloat, imageToMatch->GetCharImage());
        else
            imageFloat = milx::Image<floatImageType>::MaskImage<floatImageType>(imageFloat, imageToMatch->GetFloatImage());
    }
    emit done(-1);

    generateImage();
}

void milxQtImage::subsample(size_t xSampleFactor, size_t ySampleFactor, size_t zSampleFactor)
{
    if(usingVTKImage)
    {
        printError("Resampling from VTK image not support yet.");
        return;
    }

    bool ok1 = false, ok2 = false, ok3 = false;
    if(xSampleFactor == 0 && ySampleFactor == 0 && zSampleFactor == 0)
    {
        xSampleFactor = QInputDialog::getInt(this, tr("Please Provide the downsample factors"),
                                               tr("X Factor:"), 2, 0, 1000, 1, &ok1);
        ySampleFactor = QInputDialog::getInt(this, tr("Please Provide the downsample factors"),
                                               tr("Y Factor:"), 2, 0, 1000, 1, &ok2);
        zSampleFactor = QInputDialog::getInt(this, tr("Please Provide the downsample factors"),
                                               tr("Z Factor:"), 2, 0, 1000, 1, &ok3);
    }

    if(!ok1 || !ok2 || !ok3)
        return;

    printInfo("Subsampling Image");
    emit working(-1);
    floatImageType::SizeType factors;
    factors[0] = xSampleFactor;
    factors[1] = ySampleFactor;
    factors[2] = zSampleFactor;
    if(eightbit)
        imageChar = milx::Image<charImageType>::SubsampleImage(imageChar, factors);
    else if(integer)
        imageInt = milx::Image<intImageType>::SubsampleImage(imageInt, factors);
    else if(rgb)
        imageRGB = milx::Image<rgbImageType>::SubsampleImage(imageRGB, factors);
    else if(vectorised)
        imageVector = milx::Image<vectorImageType>::SubsampleImage(imageVector, factors);
    else
        imageFloat = milx::Image<floatImageType>::SubsampleImage(imageFloat, factors);
    emit done(-1);

    generateImage();
}

void milxQtImage::crop(QString filename)
{
    if(usingVTKImage)
    {
        printError("Masking from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename("Select Mask");

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
    printInfo("Masking and Cropping Image");
    emit working(-1);
    if(eightbit)
    {
        if(imageToMatch->is8BitImage())
            imageChar = milx::Image<charImageType>::MaskAndCropImage<charImageType>(imageChar, imageToMatch->GetCharImage());
        else
            imageChar = milx::Image<charImageType>::MaskAndCropImage<floatImageType>(imageChar, imageToMatch->GetFloatImage());
    }
    else if(integer)
    {
      if(imageToMatch->is8BitImage())
        imageInt = milx::Image<intImageType>::MaskAndCropImage<charImageType>(imageInt, imageToMatch->GetCharImage());
      else
        imageInt = milx::Image<intImageType>::MaskAndCropImage<floatImageType>(imageInt, imageToMatch->GetFloatImage());
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::MaskAndCropImage<rgbImageType>(imageRGB, imageToMatch->GetRGBImage());
    else if(vectorised)
    {
        if(imageToMatch->is8BitImage())
            imageVector = milx::Image<vectorImageType>::MaskAndCropImage<charImageType>(imageVector, imageToMatch->GetCharImage());
        else
            imageVector = milx::Image<vectorImageType>::MaskAndCropImage<floatImageType>(imageVector, imageToMatch->GetFloatImage());

        magnitude();
    }
    else
    {
        if(imageToMatch->is8BitImage())
            imageFloat = milx::Image<floatImageType>::MaskAndCropImage<charImageType>(imageFloat, imageToMatch->GetCharImage());
        else
            imageFloat = milx::Image<floatImageType>::MaskAndCropImage<floatImageType>(imageFloat, imageToMatch->GetFloatImage());
    }
    emit done(-1);
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)

    generateImage();
}

void milxQtImage::resampleLabel(QString filename)
{
    if(usingVTKImage)
    {
        printError("Resampling from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename("Select Reference Image");

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

    printInfo("Resampling as Labelled Image");
    emit working(-1);
    if(eightbit)
    {
        if(imageToMatch->is8BitImage())
            imageChar = milx::Image<charImageType>::ResampleLabel<charImageType>(imageChar, imageToMatch->GetCharImage());
        else
        {
            imageFloat = milx::Image<charImageType>::ResampleLabel<floatImageType>(imageChar, imageToMatch->GetFloatImage());
            eightbit = false;
        }
    }
    else if(integer)
    {
      if(imageToMatch->is8BitImage())
      {
           imageChar = milx::Image<intImageType>::ResampleLabel<charImageType>(imageInt, imageToMatch->GetCharImage());
           eightbit = true;
           integer = false;
      }
      else
      {
          imageFloat = milx::Image<intImageType>::ResampleLabel<floatImageType>(imageInt, imageToMatch->GetFloatImage());
          integer = false;
      }
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::ResampleLabel(imageRGB, imageToMatch->GetRGBImage());
    else
    {
        if(imageToMatch->is8BitImage())
        {
            imageChar = milx::Image<floatImageType>::ResampleLabel<charImageType>(imageFloat, imageToMatch->GetCharImage());
            eightbit = true;
        }
        else
            imageFloat = milx::Image<floatImageType>::ResampleLabel<floatImageType>(imageFloat, imageToMatch->GetFloatImage());
    }
    emit done(-1);

    generateImage();
}

void milxQtImage::transform(QString filename, QString refImgFilename, bool inverse)
{
    if(usingVTKImage)
    {
        printError("Checker Board from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
    {
        filename = getOpenFilename("Select Transform File", "ITK Transform Files (*.txt *.tfm);;VTK Transform Files (*.trsf)");

        QMessageBox msgBox;
        msgBox.setText("Use Reference Image");
        msgBox.setInformativeText("Do you want to resample to a reference image?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        int ret1 = msgBox.exec();

        if(ret1 == QMessageBox::Yes)
            refImgFilename = getOpenFilename("Select Reference Image");

        QMessageBox msgBox2;
        msgBox2.setText("Invert Transform");
        msgBox2.setInformativeText("Do you want to invert the transform?");
        msgBox2.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox2.setDefaultButton(QMessageBox::Yes);
        int ret2 = msgBox2.exec();

        if(ret2 == QMessageBox::Yes)
            inverse = true;
    }

    if(filename.isEmpty())
        return;

    typedef double transformType;
    typedef itk::Transform<transformType> TransformType;
    TransformType::Pointer transf = milx::File::OpenTransform<transformType>(filename.toStdString()); //!< Use ITK transform function

    cout << "Transform to be used: ";
    transf.Print(cout);

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = false;
    if(!refImgFilename.isEmpty())
    {
        success = reader->openImage(refImgFilename, imageToMatch);

        if(!success)
            return;
    }

    printInfo("Transforming Image with " + QString(transf->GetTransformTypeAsString().c_str()));
    emit working(-1);
    if(!refImgFilename.isEmpty())
    {
        if(eightbit)
        {
            if(imageToMatch->is8BitImage())
                imageChar = milx::Image<charImageType>::TransformImage<charImageType, TransformType, transformType>(imageChar, imageToMatch->GetCharImage(), transf, inverse, 0); //NN Interp
            else
            {
                imageFloat = milx::Image<charImageType>::TransformImage<floatImageType, TransformType, transformType>(imageChar, imageToMatch->GetFloatImage(), transf, inverse, 0); //NN interp
                eightbit = false;
            }
        }
        else if(integer)
        {
            if(imageToMatch->is8BitImage())
            {
                imageChar = milx::Image<intImageType>::TransformImage<charImageType, TransformType, transformType>(imageInt, imageToMatch->GetCharImage(), transf, inverse, 0); //NN Interp
                eightbit = true;
                integer = false;
            }
            else
            {
                imageFloat = milx::Image<intImageType>::TransformImage<floatImageType, TransformType, transformType>(imageInt, imageToMatch->GetFloatImage(), transf, inverse, 0); //NN interp
                integer = false;
            }
        }
    //    else if(rgb)
    //        imageRGB = milx::Image<rgbImageType>::TransformImage(imageRGB, imageToMatch->GetRGBImage(), transf);
        else
        {
            if(imageToMatch->is8BitImage())
            {
                imageChar = milx::Image<floatImageType>::TransformImage<charImageType, TransformType, transformType>(imageFloat, imageToMatch->GetCharImage(), transf, inverse);
                eightbit = true;
            }
            else
                imageFloat = milx::Image<floatImageType>::TransformImage<floatImageType, TransformType, transformType>(imageFloat, imageToMatch->GetFloatImage(), transf, inverse);
        }
    }
    else
    {
        printDebug("Transforming without reference image");
        if(eightbit)
        {
            if(imageToMatch->is8BitImage())
                imageChar = milx::Image<charImageType>::TransformImage<charImageType, TransformType, transformType>(imageChar, transf, inverse, 0); //NN Interp
            else
            {
                imageFloat = milx::Image<charImageType>::TransformImage<floatImageType, TransformType, transformType>(imageChar, transf, inverse, 0); //NN Interp
                eightbit = false;
            }
        }
        else if(integer)
        {
            if(imageToMatch->is8BitImage())
            {
                imageChar = milx::Image<intImageType>::TransformImage<charImageType, TransformType, transformType>(imageInt, transf, inverse);
                eightbit = true;
                integer = false;
            }
            else
            {
                imageFloat = milx::Image<intImageType>::TransformImage<floatImageType, TransformType, transformType>(imageInt, transf, inverse); 
                integer = false;
            }
        }
    //    else if(rgb)
    //        imageRGB = milx::Image<rgbImageType>::TransformImage(imageRGB, imageToMatch->GetRGBImage(), transf);
        else
        {
            if(imageToMatch->is8BitImage())
            {
                imageChar = milx::Image<floatImageType>::TransformImage<charImageType, TransformType, transformType>(imageFloat, transf, inverse);
                eightbit = true;
            }
            else
                imageFloat = milx::Image<floatImageType>::TransformImage<floatImageType, TransformType, transformType>(imageFloat, transf, inverse);
        }
    }
    emit done(-1);
    printDebug("Done Transforming");

    generateImage();
}

void milxQtImage::checkerBoard(milxQtImage *img, int numberOfSquares)
{
    if(usingVTKImage)
    {
        printError("Checker Board from VTK image not support yet.");
        return;
    }

    bool ok = false;
    if(numberOfSquares == 0)
    {
        numberOfSquares = QInputDialog::getInt(this, tr("Please Provide the number of squares"),
                                               tr("Squares:"), 10, 0, 100, 1, &ok);
    }

    if(!ok)
        return;

    printInfo("Checker Boarding Image");
    emit working(-1);
    if(actualNumberOfDimensions > 2)
    {
        if(eightbit)
            imageChar = milx::Image<charImageType>::CheckerBoard(imageChar, img->GetCharImage(), numberOfSquares);
        else if(eightbit)
            imageInt = milx::Image<intImageType>::CheckerBoard(imageInt, img->GetIntImage(), numberOfSquares);
    //    else if(rgb)
    //        imageRGB = milx::Image<rgbImageType>::CheckerBoard(imageRGB, img->GetRGBImage());
        else
            imageFloat = milx::Image<floatImageType>::CheckerBoard(imageFloat, img->GetFloatImage(), numberOfSquares);
    }
    else if(actualNumberOfDimensions == 2) //need to handle 2D images explicitly, floating point exception otherwise
    {
        int extent[6];
        viewer->GetImageActor()->GetDisplayExtent(extent);

        if(flipped)
        {
            int actualExtent[6];
            imageData->GetExtent(actualExtent);

            if(extent[3]-extent[2] == 0) //flip y extent
            {
                extent[2] = actualExtent[3]-extent[2];
                extent[3] = actualExtent[3]-extent[3];
            }
        }

        if(eightbit)
        {
            typedef itk::Image<charPixelType, 2> charImage2DType;
            charImage2DType::Pointer imageChar2D = milx::Image<charImageType>::ExtractSlice<charImage2DType>(imageChar, extent);
            charImage2DType::Pointer imageChar2DToChecker = milx::Image<charImageType>::ExtractSlice<charImage2DType>(img->GetCharImage(), extent);
            imageChar2D = milx::Image<charImage2DType>::CheckerBoard(imageChar2D, imageChar2DToChecker, numberOfSquares);
            imageChar = milx::Image<charImage2DType>::CastImage<charImageType>(imageChar2D);
        }
        else if(integer)
        {
            typedef itk::Image<intPixelType, 2> intImage2DType;
            intImage2DType::Pointer imageInt2D = milx::Image<charImageType>::ExtractSlice<intImage2DType>(imageChar, extent);
            intImage2DType::Pointer imageInt2DToChecker = milx::Image<intImageType>::ExtractSlice<intImage2DType>(img->GetIntImage(), extent);
            imageInt2D = milx::Image<intImage2DType>::CheckerBoard(imageInt2D, imageInt2DToChecker, numberOfSquares);
            imageInt = milx::Image<intImage2DType>::CastImage<intImageType>(imageInt2D);
        }
    //    else if(rgb)
    //        imageRGB = milx::Image<rgbImageType>::CheckerBoard(imageRGB, img->GetRGBImage());
        else
        {
            typedef itk::Image<floatPixelType, 2> floatImage2DType;
            floatImage2DType::Pointer imageFloat2D = milx::Image<floatImageType>::ExtractSlice<floatImage2DType>(imageFloat, extent);
            floatImage2DType::Pointer imageFloat2DToChecker = milx::Image<floatImageType>::ExtractSlice<floatImage2DType>(img->GetFloatImage(), extent);
            imageFloat2D = milx::Image<floatImage2DType>::CheckerBoard(imageFloat2D, imageFloat2DToChecker, numberOfSquares);
            imageFloat = milx::Image<floatImage2DType>::CastImage<floatImageType>(imageFloat2D);
        }
    }
    else
    {
        printError("Actual Dimensions of data is" + QString::number(actualNumberOfDimensions) + ", which is not supported. Ignoring.");
    }
    emit done(-1);

    generateImage();
}

void milxQtImage::checkerBoard(QString filename, int numberOfSquares)
{
    if(usingVTKImage)
    {
        printError("Checker Board from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToMatch = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToMatch);

    if(!success)
        return;

    checkerBoard(imageToMatch, numberOfSquares);
}

void milxQtImage::distanceMap(bool signedDistance, bool inside)
{
    if(usingVTKImage)
    {
        printError("Distance Map from VTK image not support yet.");
        return;
    }

    /*if(signedDistance && inside)
    {
        ///Ask if use signed distances
        QMessageBox msgBox;
        msgBox.setText("Distance Map Type");
        msgBox.setInformativeText("Do you wish to compute signed distances?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        int ret = msgBox.exec();
        if(ret == QMessageBox::Yes)
        {
          msgBox.setText("Distance Map Inside or Outside");
          msgBox.setInformativeText("Do you wish to compute inside (the object) distances?");
          msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
          msgBox.setDefaultButton(QMessageBox::No);
        }
        else
          signedDistance = false;
        int retInside = msgBox.exec();
    }*/

    emit working(-1);
    if(eightbit)
    {
        imageFloat = milx::Image<charImageType>::DistanceMap<floatImageType>(imageChar, true, signedDistance, inside);
        eightbit = false;
    }
    else if(integer)
    {
        imageFloat = milx::Image<intImageType>::DistanceMap<floatImageType>(imageInt, true, signedDistance, inside);
        integer = false;
    }
//    else if(rgb)
//        imageRGB = milx::Image::DistanceMap<rgbImageType>(imageRGB);
    else
    {
        imageFloat = milx::Image<floatImageType>::DistanceMap<floatImageType>(imageFloat, false, signedDistance, inside);
    }
    emit done(-1);

    generateImage();
}

void milxQtImage::thresholdAbove(float value, float level)
{
    if(usingVTKImage)
    {
        printError("Threshold above for VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false);

    bool ok1 = false, ok2 = false;
    if(value == 0 && level == 0)
    {
        value = QInputDialog::getDouble(this, tr("Please Provide Outside Value"),
                                        tr("Outside Value:"), 0, minValue, maxValue, 5, &ok1);
        level = QInputDialog::getDouble(this, tr("Please Provide the threshold upper level"),
                                        tr("Level:"), maxValue, minValue, maxValue, 5, &ok2);

        if(!ok1 || !ok2)
            return;
    }

    printInfo("Thresholding Image from above");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::ThresholdAboveImage(imageChar, value, level);
    else if(integer)
        imageInt = milx::Image<intImageType>::ThresholdAboveImage(imageInt, value, level);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::ThresholdAboveImage(imageRGB);
    else
        imageFloat = milx::Image<floatImageType>::ThresholdAboveImage(imageFloat, value, level);
    emit done(-1);

    generateImage();
}

void milxQtImage::thresholdBelow(float value, float level)
{
    if(usingVTKImage)
    {
        printError("Threshold below for VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false);

    bool ok1 = false, ok2 = false;
    if(value == 0 && level == 0)
    {
        value = QInputDialog::getDouble(this, tr("Please Provide Outside Value"),
                                        tr("Outside Value:"), 0, minValue, maxValue, 5, &ok1);
        level = QInputDialog::getDouble(this, tr("Please Provide the threshold lower level"),
                                        tr("Level:"), minValue, minValue, maxValue, 5, &ok2);

        if(!ok1 || !ok2)
            return;
    }

    printInfo("Thresholding Image from below");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::ThresholdBelowImage(imageChar, value, level);
    else if(integer)
        imageInt = milx::Image<intImageType>::ThresholdBelowImage(imageInt, value, level);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::ThresholdBelowImage(imageRGB);
    else
        imageFloat = milx::Image<floatImageType>::ThresholdBelowImage(imageFloat, value, level);
    emit done(-1);

    generateImage();
}

void milxQtImage::threshold(float value, float blevel, float alevel)
{
    if(usingVTKImage)
    {
        printError("Threshold for VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false);

    bool ok1 = false, ok2 = false, ok3 = false;
    if(value == 0 && blevel == 0 && alevel == 0)
    {
        value = QInputDialog::getDouble(this, tr("Please Provide Outside Value"),
                                        tr("Outside Value:"), 0, minValue, maxValue, 5, &ok1);
        blevel = QInputDialog::getDouble(this, tr("Please Provide the threshold lower level"),
                                         tr("Lower Level:"), minValue, minValue, maxValue, 5, &ok2);
        alevel = QInputDialog::getDouble(this, tr("Please Provide the threshold upper level"),
                                         tr("Upper Level:"), maxValue, minValue, maxValue, 5, &ok3);

        if(!ok1 || !ok2 || !ok3)
            return;
    }

    printInfo("Thresholding Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::ThresholdImage(imageChar, value, blevel, alevel);
    else if(integer)
        imageInt = milx::Image<intImageType>::ThresholdImage(imageInt, value, blevel, alevel);
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::ThresholdImage(imageRGB, blevel, alevel);
    else
        imageFloat = milx::Image<floatImageType>::ThresholdImage(imageFloat, value, blevel, alevel);
    emit done(-1);

    generateImage();
}

void milxQtImage::otsu(int bins)
{
    if(usingVTKImage)
    {
        printError("Otsu Threshold for VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false);

    bool ok = false;
    if(bins == 0)
    {
      bins = QInputDialog::getInt(this, tr("Please Provide the histogram bins to use"),
                                        tr("Bins:"), 128, 0, 8192, 1, &ok);
    }

    if(!ok)
        return;

    printInfo("Otsu Thresholding Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::OtsuThresholdImage<charImageType>(imageChar, bins);
    else if(integer)
        imageInt = milx::Image<intImageType>::OtsuThresholdImage<intImageType>(imageInt, bins);
    //~ else if(rgb)
        //~ imageChar = milx::Image<rgbImageType>::OtsuThresholdImage<charImageType>(imageRGB, bins);
    else
        imageChar = milx::Image<floatImageType>::OtsuThresholdImage<charImageType>(imageFloat, bins);
    emit done(-1);

    eightbit = true;
    rgb = false;
    generateImage();
}

void milxQtImage::otsuMultiple(int bins, int labels)
{
    if(usingVTKImage)
    {
        printError("Otsu Multiple Threshold for VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false);

    bool ok = false;
    if(bins == 0)
    {
      bins = QInputDialog::getInt(this, tr("Please Provide the histogram bins to use"),
                                        tr("Bins:"), 128, 0, 8192, 1, &ok);
      labels = QInputDialog::getInt(this, tr("Please Provide the number of labels to produce"),
                                        tr("Labels:"), 1, 0, 8192, 1, &ok);
    }

    if(!ok)
        return;

    printInfo("Otsu Multiple Thresholding Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::OtsuMultipleThresholdImage<charImageType>(imageChar, bins, labels);
    else if(integer)
        imageInt = milx::Image<intImageType>::OtsuMultipleThresholdImage<intImageType>(imageInt, bins, labels);
    //~ else if(rgb)
        //~ imageChar = milx::Image<rgbImageType>::OtsuMultipleThresholdImage<charImageType>(imageRGB, bins, labels);
    else
        imageChar = milx::Image<floatImageType>::OtsuMultipleThresholdImage<charImageType>(imageFloat, bins, labels);
    emit done(-1);

    eightbit = true;
    rgb = false;
    generateImage();
}

void milxQtImage::binaryThreshold(float value, float blevel, float alevel)
{
    if(usingVTKImage)
    {
        printError("Binary Threshold for VTK image not support yet.");
        return;
    }

    histogram(256, 0, 255, false);

    bool ok1 = false, ok2 = false, ok3 = false;
    if(value == 0 && blevel == 0 && alevel == 0)
    {
        value = QInputDialog::getDouble(this, tr("Please Provide Inside Value"),
                                        tr("Inside Value:"), 1.0, 0, 255, 1, &ok1);
        blevel = QInputDialog::getDouble(this, tr("Please Provide the threshold lower level"),
                                         tr("Lower Level:"), minValue, minValue, maxValue, 5, &ok2);
        alevel = QInputDialog::getDouble(this, tr("Please Provide the threshold upper level"),
                                         tr("Upper Level:"), maxValue, minValue, maxValue, 5, &ok3);

        if(!ok1 || !ok2 || !ok3)
            return;
    }

    printInfo("Binary Thresholding Image");
    emit working(-1);
    if(eightbit)
        imageChar = milx::Image<charImageType>::BinaryThresholdImage<charImageType>(imageChar, 0, value, blevel, alevel);
    else if(integer)
        imageChar = milx::Image<intImageType>::BinaryThresholdImage<charImageType>(imageInt, 0, value, blevel, alevel);
//    else if(rgb)
//        imageChar = milx::Image<rgbImageType>::TBinaryThresholdImage<charImageType>(imageRGB, 0, value, blevel, alevel);
    else
        imageChar = milx::Image<floatImageType>::BinaryThresholdImage<charImageType>(imageFloat, 0, value, blevel, alevel);
    emit done(-1);

    eightbit = true;
    rgb = false;
    generateImage();
}

void milxQtImage::flip(bool xAxis, bool yAxis, bool zAxis, bool aboutOrigin)
{
    if(usingVTKImage)
    {
        printError("Flipping from VTK image not support yet.");
        return;
    }

    if(!xAxis && !yAxis && !zAxis)
    {
        QMessageBox msgBoxX, msgBoxY, msgBoxZ, msgBoxOrigin;
        msgBoxX.setText("The image will be flipped");
        msgBoxX.setInformativeText("Do you want to flip the x-axis?");
        msgBoxX.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxX.setDefaultButton(QMessageBox::No);
        msgBoxY.setText("The image will be flipped");
        msgBoxY.setInformativeText("Do you want to flip the y-axis?");
        msgBoxY.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxY.setDefaultButton(QMessageBox::Yes);
        msgBoxZ.setText("The image will be flipped");
        msgBoxZ.setInformativeText("Do you want to flip the z-axis?");
        msgBoxZ.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxZ.setDefaultButton(QMessageBox::No);
        msgBoxOrigin.setText("The image will be flipped");
        msgBoxOrigin.setInformativeText("Do you want to flip about the origin?");
        msgBoxOrigin.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxOrigin.setDefaultButton(QMessageBox::Yes);
        int retX = msgBoxX.exec();
        int retY = msgBoxY.exec();
        int retZ = msgBoxZ.exec();
        int retOrigin = msgBoxOrigin.exec();

        if(retX == QMessageBox::Yes)
            xAxis = true;
        if(retY == QMessageBox::Yes)
            yAxis = true;
        if(retZ == QMessageBox::Yes)
            zAxis = true;
        if(retOrigin == QMessageBox::Yes)
            aboutOrigin = true;
        else
            aboutOrigin = false;
    }

    printInfo("Flipping Image Data");
    emit working(-1);
    //~ floatImageType::DirectionType direction;
    if(eightbit)
    {
        imageChar = milx::Image<charImageType>::FlipImage(imageChar, xAxis, yAxis, zAxis, aboutOrigin);
        //~ direction = imageChar->GetDirection();
        //~ cerr << "Flipped Direction: " << imageChar->GetDirection() << endl;
        //~ flipped = !flipped;
        //~ imageChar->GetDirection()(1,1) *= -1;
    }
    else if(integer)
    {
        imageInt = milx::Image<intImageType>::FlipImage(imageInt, xAxis, yAxis, zAxis, aboutOrigin);
        //~ direction = imageInt->GetDirection();
        //~ cout << "Flipped Direction: " << imageInt->GetDirection() << endl;
        //~ flipped = !flipped;
        //~ imageInt->GetDirection()(1,1) *= -1;
    }
    else if(rgb)
    {
        imageRGB = milx::Image<rgbImageType>::FlipImage(imageRGB, xAxis, yAxis, zAxis, aboutOrigin);
        //~ direction = imageRGB->GetDirection();
        //~ flipped = !flipped;
        //~ imageRGB->GetDirection()(1,1) *= -1;
    }
    else
    {
        imageFloat = milx::Image<floatImageType>::FlipImage(imageFloat, xAxis, yAxis, zAxis, aboutOrigin);
        //~ direction = imageFloat->GetDirection();
        //~ cerr << "Flipped Direction: " << imageFloat->GetDirection() << endl;
        //~ flipped = !flipped;
        //~ imageFloat->GetDirection()(1,1) *= -1;
    }

    //~ imageFloat::DirectionType flipMatrix;
    //~ flipMatrix->SetIdentity();

    emit done(-1);

    generateImage();
    viewer->GetRenderer()->ResetCamera();
}

void milxQtImage::surface(const float value)
{
    updateData(orientAct->isChecked());

    emit imageToSurface(imageData, value);
}

void milxQtImage::polyData()
{
    updateData(orientAct->isChecked());

    emit imageToPolyData(imageData);
}

void milxQtImage::magnitude()
{
    if(!vectorised)
        return;

    updateData(orientAct->isChecked());

    //For display we want to show the vector magnitudes only
    printInfo("Loaded image as vector image.");
    printDebug("Computing magnitude of field image.");
    imageFloat = milx::Image<vectorImageType>::VectorMagnitude<floatImageType>(imageVector);

    generateImage();
}

void milxQtImage::component(int index)
{
    if(!vectorised)
        return;

    updateData(orientAct->isChecked());

    bool ok1 = false;
    if(index < 0)
    {
        index = QInputDialog::getInt(this, tr("Please Provide the component index to display"),
                                          tr("Component:"), 0, 0, imageVector->GetNumberOfComponentsPerPixel(), 1, &ok1);
    }

    if(!ok1)
        return;

    imageFloat = milx::Image<vectorImageType>::ExtractComponent<floatImageType>(imageVector, index);

    generateImage();
}

void milxQtImage::pseudoImage()
{
    if(!vectorised)
        return;

    updateData(orientAct->isChecked());

    emit imageToPseudoImage(imageVector);
}

void milxQtImage::vectorField(int subsampleFactor, float scaling)
{
    if(!vectorised)
        return;

    updateData(orientAct->isChecked());

    const size_t components = imageVector->GetNumberOfComponentsPerPixel();
    int ret = QMessageBox::No;
    if(components > 3)
    {
        QMessageBox msgBox;
        msgBox.setText("Tensor Field?");
        msgBox.setInformativeText("Found more than 3 components. Do you wish to compute a tensor field?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        ret = msgBox.exec();
    }

//    if(components == 6 || components == 9)
    if(ret == QMessageBox::Yes)
        emit imageToTensorField(imageVector, imageFloat, subsampleFactor, scaling);
    else
        emit imageToVectorField(imageVector, imageFloat, subsampleFactor, scaling);
}

void milxQtImage::streamLines()
{
    if(!vectorised)
        return;

    updateData(orientAct->isChecked());

    const size_t components = imageVector->GetNumberOfComponentsPerPixel();
    if(components != 3)
    {
        printError("Can only integrate vector fields. Number of Components must be 3.");
        return;
    }

    int extent[6];
    viewer->GetImageActor()->GetDisplayExtent(extent);

    if(flipped)
    {
        int actualExtent[6];
        imageData->GetExtent(actualExtent);

        if(extent[3]-extent[2] == 0) //flip y extent
        {
            extent[2] = actualExtent[3]-extent[2];
            extent[3] = actualExtent[3]-extent[3];
        }
    }

    ///Extract the viewed slice
    floatImageType::Pointer sliceFloat = milx::Image<floatImageType>::ExtractSlice<floatImageType>(imageFloat, extent);

    emit imageToStreamLines(imageVector, sliceFloat);
}

void milxQtImage::anisotropicDiffusion()
{
    if(usingVTKImage)
    {
        printError("Anisotropic Diffusion from VTK image not support yet.");
        return;
    }

    //auto set timestep initially
    floatImageType::SizeType imgSize;
    if(eightbit)
        imgSize = imageChar->GetLargestPossibleRegion().GetSize();
    else
        imgSize = imageFloat->GetLargestPossibleRegion().GetSize();

    size_t maxDimension = 0;
    for(size_t k = 0; k < imgSize.GetSizeDimension(); k ++)
    {
      if(imgSize[k] > maxDimension)
        maxDimension = imgSize[k];
    }
    printDebug("Timestep will be initially the reciprocal of " + QString::number(maxDimension));
    float timestep = 2.0/maxDimension;

    bool ok1 = false, ok2 = false;
    int iterations = QInputDialog::getInt(this, tr("Please Provide the number of Iterations"),
                                          tr("Iterations:"), 5, 0, 1000, 1, &ok1);
    timestep = QInputDialog::getDouble(this, tr("Please Provide the timestep"),
                     tr("Timestep:"), timestep, 0.0, 100.0, 5, &ok2);

    if(ok1 && ok2)
    {
        printInfo("Computing Anisotropic Diffusion of Image");
        emit working(-1);
        if(eightbit)
        {
            imageFloat = milx::Image<charImageType>::AnisotropicDiffusion<floatImageType>(imageChar, iterations, timestep);
            eightbit = false;
        }
        else if(integer)
        {
            imageFloat = milx::Image<intImageType>::AnisotropicDiffusion<floatImageType>(imageInt, iterations, timestep);
            integer = false;
        }
//        else if(rgb)
//            imageRGB = milx::Image<rgbImageType>::AnisotropicDiffusion(imageRGB, iterations, timestep);
        else
            imageFloat = milx::Image<floatImageType>::AnisotropicDiffusion<floatImageType>(imageFloat, iterations, timestep);
        emit done(-1);

        generateImage();
    }
}

void milxQtImage::gaussianSmooth()
{
    if(usingVTKImage)
    {
        printError("Gaussian Smoothin from VTK image not support yet.");
        return;
    }

    bool ok1 = false;
    float variance = QInputDialog::getDouble(this, tr("Please Provide the variance of the Gaussian to use"),
                     tr("Variance:"), 0.5, 0.0, 2147483647, 5, &ok1);

    if(ok1)
    {
        printInfo("Computing Gaussian Smoothing of Image");
        emit working(-1);
        if(eightbit)
            imageChar = milx::Image<charImageType>::GaussianSmooth(imageChar, variance);
        else if(integer)
            imageInt = milx::Image<intImageType>::GaussianSmooth(imageInt, variance);
//        else if(rgb)
//            imageRGB = milx::Image<rgbImageType>::GaussianSmooth(imageRGB, variance);
        else
            imageFloat = milx::Image<floatImageType>::GaussianSmooth(imageFloat, variance);
        emit done(-1);

        generateImage();
    }
}

void milxQtImage::bilateral()
{
  if(usingVTKImage)
    {
      printError("Bilateral Smoothin from VTK image not support yet.");
      return;
    }

  bool ok1 = false, ok2 = false;
  float sigmaRange = QInputDialog::getDouble(this, tr("Please Provide the range sigma to use"),
                                           tr("Range Sigma:"), 0.5, 0.0, 2147483647, 5, &ok1);
  float sigmaSpatial = QInputDialog::getDouble(this, tr("Please Provide the domain/spatial sigma to use"),
                                           tr("Domain Sigma:"), 5, 0.0, 2147483647, 5, &ok2);

  if(ok1)
  {
      printInfo("Computing Bilateral Smoothing of Image");
      emit working(-1);
      if(eightbit)
          imageChar = milx::Image<charImageType>::Bilateral(imageChar, sigmaRange, sigmaSpatial);
      else if(integer)
          imageInt = milx::Image<intImageType>::Bilateral(imageInt, sigmaRange, sigmaSpatial);
      //        else if(rgb)
      //            imageRGB = milx::Image<rgbImageType>::Bilateral(imageRGB, sigmaRange, sigmaSpatial);
      else
          imageFloat = milx::Image<floatImageType>::Bilateral(imageFloat, sigmaRange, sigmaSpatial);
      emit done(-1);

      generateImage();
  }
}

void milxQtImage::median()
{
    if(usingVTKImage)
    {
        printError("Median from VTK image not support yet.");
        return;
    }

    bool ok = false;
    int radius = QInputDialog::getInt(this, tr("Please Provide the radius to use"),
                                      tr("Radius:"), 5, 0, 1000, 1, &ok);

    if(ok)
    {
        printInfo("Computing Median of Image");
        emit working(-1);
        if(eightbit)
            imageChar = milx::Image<charImageType>::Median(imageChar, radius);
        else if(integer)
            imageInt = milx::Image<intImageType>::Median(imageInt, radius);
//        else if(rgb)
//            imageRGB = milx::Image<rgbImageType>::Median(imageRGB, radius);
        else
            imageFloat = milx::Image<floatImageType>::Median(imageFloat, radius);
        emit done(-1);

        generateImage();
    }
}

void milxQtImage::zeros(const unsigned long xSize, const unsigned long ySize, const unsigned long zSize, milxQtImage *refImage)
{
    if(usingVTKImage)
    {
        printError("Adding of VTK image not support yet.");
        return;
    }

    printInfo("Blank Image");
    emit working(-1);
    if(eightbit)
    {
        charImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        imageChar = milx::Image<charImageType>::BlankImage(0.0, blankSize);

        if(refImage)
        {
            imageChar->SetOrigin(refImage->GetCharImage()->GetOrigin());
            imageChar->SetSpacing(refImage->GetCharImage()->GetSpacing());
            imageChar->SetDirection(refImage->GetCharImage()->GetDirection());
        }
    }
    else if(integer)
    {
        intImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        imageInt = milx::Image<intImageType>::BlankImage(0.0, blankSize);

        if(refImage)
        {
            imageInt->SetOrigin(refImage->GetIntImage()->GetOrigin());
            imageInt->SetSpacing(refImage->GetIntImage()->GetSpacing());
            imageInt->SetDirection(refImage->GetIntImage()->GetDirection());
        }
    }
    else if(rgb)
    {
        rgbImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        imageRGB = milx::Image<rgbImageType>::BlankImage(0.0, blankSize);

        if(refImage)
        {
            imageRGB->SetOrigin(refImage->GetRGBImage()->GetOrigin());
            imageRGB->SetSpacing(refImage->GetRGBImage()->GetSpacing());
            imageRGB->SetDirection(refImage->GetRGBImage()->GetDirection());
        }
    }
    else if(!eightbit && !rgb)
    {
        floatImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        imageFloat = milx::Image<floatImageType>::BlankImage(0.0, blankSize);

        if(refImage)
        {
            imageFloat->SetOrigin(refImage->GetFloatImage()->GetOrigin());
            imageFloat->SetSpacing(refImage->GetFloatImage()->GetSpacing());
            imageFloat->SetDirection(refImage->GetFloatImage()->GetDirection());
        }
    }
    else
        printError("Blank image of different types not supported.");
    emit done(-1);

    generateImage();
}

void milxQtImage::resize(double outputSpacing)
{
    floatImageType::SpacingType newSpacing;
    floatImageType::SizeType newSize;
    floatImageType::PointType newOrigin;
    floatImageType::DirectionType newDirection;
    if(eightbit)
    {
        newSpacing = imageChar->GetSpacing();
        newSize = imageChar->GetLargestPossibleRegion().GetSize();
        newOrigin = imageChar->GetOrigin();
        newDirection = imageChar->GetDirection();
    }
    else if(integer)
    {
        newSpacing = imageInt->GetSpacing();
        newSize = imageInt->GetLargestPossibleRegion().GetSize();
        newOrigin = imageInt->GetOrigin();
        newDirection = imageInt->GetDirection();
    }
    else if(!eightbit && !rgb && !vectorised)
    {
        newSpacing = imageFloat->GetSpacing();
        newSize = imageFloat->GetLargestPossibleRegion().GetSize();
        newOrigin = imageFloat->GetOrigin();
        newDirection = imageFloat->GetDirection();
    }
    else
    {
        printError("Resize image of different types not supported.");
        return;
    }

    bool ok1 = false;
    if(outputSpacing == 0.0) //none provided
    {
        outputSpacing = QInputDialog::getDouble(this, tr("Please Provide the new isotropic spacing"),
                                         tr("Isotropic Spacing:"), 0.5, 0.001, 10, 3, &ok1);

        if(!ok1)
            return;
    }

    printInfo("Resizing Image based on new spacing");
    emit working(-1);
    floatImageType::SpacingType scaleFactor;
      scaleFactor[0] = newSpacing[0]/outputSpacing;
      scaleFactor[1] = newSpacing[1]/outputSpacing;
      scaleFactor[2] = newSpacing[2]/outputSpacing;
      newSpacing.Fill(outputSpacing);
      newSize[0] = static_cast<unsigned>(newSize[0]*scaleFactor[0]+0.5);  // size along X
      newSize[1] = static_cast<unsigned>(newSize[1]*scaleFactor[1]+0.5);  // size along Y
      newSize[2] = static_cast<unsigned>(newSize[2]*scaleFactor[2]+0.5);  // size along Z
    std::cout << "Rescaling image by [" << scaleFactor[0] << ", " << scaleFactor[1] << ", " << scaleFactor[2] << "]" << std::endl;
    std::cout << "Resampling image to [" << newSize[0] << ", " << newSize[1] << ", " << newSize[2] << "]" << std::endl;

    if(eightbit)
        imageChar = milx::Image<charImageType>::ResizeImage(imageChar, newSize, newSpacing, newOrigin, newDirection);
    else if(integer)
        imageInt = milx::Image<intImageType>::ResizeImage(imageInt, newSize, newSpacing, newOrigin, newDirection);
    else if(!eightbit && !rgb && !vectorised)
        imageFloat = milx::Image<floatImageType>::ResizeImage(imageFloat, newSize, newSpacing, newOrigin, newDirection);
    emit done(-1);

    generateImage();
}

void milxQtImage::resize(const unsigned long xSize, const unsigned long ySize, const unsigned long zSize, milxQtImage *refImage)
{
    if(usingVTKImage)
    {
        printError("Adding of VTK image not support yet.");
        return;
    }

    printInfo("Resizing Image");
    emit working(-1);
    if(eightbit)
    {
        charImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        if(refImage)
            imageChar = milx::Image<charImageType>::ResizeImage(imageChar, blankSize, refImage->GetCharImage()->GetSpacing(), refImage->GetCharImage()->GetOrigin(), refImage->GetCharImage()->GetDirection());
        else
            imageChar = milx::Image<charImageType>::ResizeImage(imageChar, blankSize, imageChar->GetSpacing(), imageChar->GetOrigin(), imageChar->GetDirection());
    }
    else if(integer)
    {
        intImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        if(refImage)
            imageInt = milx::Image<intImageType>::ResizeImage(imageInt, blankSize, refImage->GetIntImage()->GetSpacing(), refImage->GetIntImage()->GetOrigin(), refImage->GetIntImage()->GetDirection());
        else
            imageInt = milx::Image<intImageType>::ResizeImage(imageInt, blankSize, imageInt->GetSpacing(), imageInt->GetOrigin(), imageInt->GetDirection());
    }
//    else if(rgb)
//    {
//        rgbImageType::SizeType blankSize;
//          blankSize[0]  = xSize;  // size along X
//          blankSize[1]  = ySize;  // size along Y
//          blankSize[2]  = zSize;  // size along Z
//
//        if(refImage)
//            imageRGB = milx::Image<rgbImageType>::ResizeImage(imageRGB, blankSize, refImage->GetRGBImage()->GetSpacing());
//        else
//            imageRGB = milx::Image<rgbImageType>::ResizeImage(imageRGB, blankSize, imageRGB->GetSpacing());
//    }
    else if(!eightbit && !rgb && !vectorised)
    {
        floatImageType::SizeType blankSize;
          blankSize[0]  = xSize;  // size along X
          blankSize[1]  = ySize;  // size along Y
          blankSize[2]  = zSize;  // size along Z

        if(refImage)
            imageFloat = milx::Image<floatImageType>::ResizeImage(imageFloat, blankSize, refImage->GetFloatImage()->GetSpacing(), refImage->GetCharImage()->GetOrigin(), refImage->GetCharImage()->GetDirection());
        else
            imageFloat = milx::Image<floatImageType>::ResizeImage(imageFloat, blankSize, imageFloat->GetSpacing(), imageFloat->GetOrigin(), imageFloat->GetDirection());
    }
    else
        printError("Blank image of different types not supported.");
    emit done(-1);

    generateImage();
}

void milxQtImage::add(milxQtImage *img)
{
    if(usingVTKImage)
    {
        printError("Adding of VTK image not support yet.");
        return;
    }

    printInfo("Adding Image");
    emit working(-1);
    if(eightbit && img->is8BitImage())
        imageChar = milx::Image<charImageType>::AddImages(imageChar, img->GetCharImage());
    else if(integer && img->is32BitImage())
        imageInt = milx::Image<intImageType>::AddImages(imageInt, img->GetIntImage());
    else if(rgb && img->isRGBImage())
        imageRGB = milx::Image<rgbImageType>::AddImages(imageRGB, img->GetRGBImage());
    else if(vectorised && img->isVectorImage())
    {
        imageVector = milx::Image<vectorImageType>::AddImages(imageVector, img->GetVectorImage());
        magnitude();
    }
    else if(!eightbit && !rgb && img->isFloatingPointImage())
        imageFloat = milx::Image<floatImageType>::AddImages(imageFloat, img->GetFloatImage());
    else
        printError("Adding images of different types not supported.");
    emit done(-1);

    generateImage();
}

void milxQtImage::add(QString filename)
{
    if(usingVTKImage)
    {
        printError("Adding from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToSubtract = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToSubtract);

    if(!success)
        return;

    add(imageToSubtract);
}

void milxQtImage::subtract(milxQtImage *img)
{
    if(usingVTKImage)
    {
        printError("Subtracting of VTK image not support yet.");
        return;
    }

    printInfo("Subtracting Image");
    emit working(-1);
    if(eightbit && img->is8BitImage())
        imageChar = milx::Image<charImageType>::SubtractImages(imageChar, img->GetCharImage());
    else if(integer && img->is32BitImage())
        imageInt = milx::Image<intImageType>::SubtractImages(imageInt, img->GetIntImage());
    else if(rgb && img->isRGBImage())
        imageRGB = milx::Image<rgbImageType>::SubtractImages(imageRGB, img->GetRGBImage());
    else if(vectorised && img->isVectorImage())
    {
        imageVector = milx::Image<vectorImageType>::SubtractImages(imageVector, img->GetVectorImage());
        magnitude();
    }
    else if(!eightbit && !rgb && img->isFloatingPointImage())
        imageFloat = milx::Image<floatImageType>::SubtractImages(imageFloat, img->GetFloatImage());
    else
        printError("Subtracting images of different types not supported.");
    emit done(-1);

    generateImage();
}

void milxQtImage::subtract(QString filename)
{
    if(usingVTKImage)
    {
        printError("Subtracting from VTK image not support yet.");
        return;
    }

    if(filename.isEmpty())
        filename = getOpenFilename();

    if(filename.isEmpty())
        return;

    QPointer<milxQtImage> imageToSubtract = new milxQtImage;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openImage(filename, imageToSubtract);

    if(!success)
        return;

    subtract(imageToSubtract);
}

void milxQtImage::multiply(milxQtImage *img)
{
  if (usingVTKImage)
  {
      printError("Multiplying of VTK image not support yet.");
      return;
  }

  printInfo("Multiplying Image");
  emit working(-1);
  if (eightbit && img->is8BitImage())
      imageChar = milx::Image<charImageType>::MultiplyImages(imageChar, img->GetCharImage());
  else if (integer && img->is32BitImage())
      imageInt = milx::Image<intImageType>::MultiplyImages(imageInt, img->GetIntImage());
  else if (!eightbit && !rgb && img->isFloatingPointImage())
      imageFloat = milx::Image<floatImageType>::MultiplyImages(imageFloat, img->GetFloatImage());
  else
      printError("Multiplying images of Vector, RGB or different types not supported.");
  emit done(-1);

  generateImage();
}

void milxQtImage::multiply(QString filename)
{
  if (usingVTKImage)
  {
      printError("Multiplying from VTK image not support yet.");
      return;
  }

  if (filename.isEmpty())
      filename = getOpenFilename();

  if (filename.isEmpty())
      return;

  QPointer<milxQtImage> imageToMultiply = new milxQtImage;
  QPointer<milxQtFile> reader = new milxQtFile;
  bool success = reader->openImage(filename, imageToMultiply);

  if (!success)
      return;

  multiply(imageToMultiply);
}

void milxQtImage::scale(float scaling)
{
    if(usingVTKImage)
    {
        printError("Scaling of VTK image not support yet.");
        return;
    }

    printInfo("Scaling Image");
    emit working(-1);
    if(eightbit)
    {
        imageFloat = milx::Image<charImageType>::ScaleImage<floatImageType>(imageChar, scaling);
        eightbit = false;
    }
    else if(integer)
    {
        imageFloat = milx::Image<intImageType>::ScaleImage<floatImageType>(imageInt, scaling);
        integer = false;
    }
//    else if(rgb)
//        imageRGB = milx::Image<rgbImageType>::ScaleImage(imageRGB);
    else if(vectorised)
    {
        imageVector = milx::Image<vectorImageType>::ScaleVectorImage(imageVector, scaling, imageVector->GetNumberOfComponentsPerPixel());
        magnitude();
    }
    else
        imageFloat = milx::Image<floatImageType>::ScaleImage<floatImageType>(imageFloat, scaling);
    emit done(-1);

    generateImage();
}

void milxQtImage::convolve(milxQtImage *img)
{
    if(usingVTKImage)
    {
        printError("Convolving of VTK images not support yet.");
        return;
    }
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
    printInfo("Convolving Images");
    emit working(-1);
    if(eightbit && img->is8BitImage())
        imageChar = milx::Image<charImageType>::ConvolveImages(imageChar, img->GetCharImage());
    else if(integer && img->is32BitImage())
        imageInt = milx::Image<intImageType>::ConvolveImages(imageInt, img->GetIntImage());
//    else if(rgb && img->isRGBImage())
//        imageRGB = milx::Image<rgbImageType>::ConvolveImages(imageRGB, img->GetRGBImage());
//    else if(vectorised && img->isVectorImage())
//    {
//        imageVector = milx::Image<vectorImageType>::ConvolveImages(imageVector, img->GetVectorImage());
//        magnitude();
//    }
    else if(!eightbit && !rgb && img->isFloatingPointImage())
        imageFloat = milx::Image<floatImageType>::ConvolveImages(imageFloat, img->GetFloatImage());
    else
        printError("Convolving images of different types not supported.");
    emit done(-1);

    generateImage();
#endif
}

void milxQtImage::highpass()
{
    imageData->DeepCopy( butterWorthHighPass(imageData) );

    typedef itk::VTKImageToImageFilter<floatImageType> ConvertImageType;

    ConvertImageType::Pointer convertFilter = ConvertImageType::New();
    convertFilter->SetInput(imageData);
    convertFilter->AddObserver(itk::ProgressEvent(), observeProgress);
    try
    {
        convertFilter->Update();
    }
    catch (itk::ExceptionObject & ex )
    {
        printError("Failed Converting VTK Image to ITK Image");
        printError(ex.GetDescription());
    }

    imageFloat = milx::Image<floatImageType>::DuplicateImage(convertFilter->GetOutput());
    flipped = false;

    generateImage();
}

void milxQtImage::interpolateDisplay(const bool quietly)
{
    if(loaded && viewerSetup)
    {
        if(viewer->GetImageActor()->GetInterpolate())
        {
            printInfo("Disabling Interpolation");
            viewer->GetImageActor()->InterpolateOff();
            interpolateAct->setChecked(false);
        }
        else
        {
            printInfo("Enabling Interpolation");
            viewer->GetImageActor()->InterpolateOn();
            interpolateAct->setChecked(true);
        }

        emit milxQtRenderWindow::modified(GetImageActor());
        if(!quietly)
            emit modified(this);
    }
}

void milxQtImage::applyOrientDisplay(const bool quietly)
{
    if(loaded && viewerSetup)
    {
        generateImage();
        viewer->GetRenderer()->ResetCamera();

        emit milxQtRenderWindow::modified(GetImageActor()); //caution: crash without it
        if(!quietly)
            emit modified(this);
    }
}

void milxQtImage::setDefaultOrientation(int orientMode)
{
    if(orientMode == RADIOLOGICAL)
        viewer->NeurologicalViewOff();
    if(orientMode == NEUROLOGICAL)
        viewer->NeurologicalViewOn();

    milxQtRenderWindow::setDefaultOrientation(orientMode);
}

void milxQtImage::enableScale(QString title, const bool quiet, double minRange, double maxRange, int noOfLabels)
{
    if(!milxQtRenderWindow::scale)
        milxQtRenderWindow::scale = vtkSmartPointer<vtkScalarBarActor>::New();
    if(!milxQtRenderWindow::scalarBar)
        milxQtRenderWindow::scalarBar = vtkSmartPointer<vtkScalarBarWidget>::New();

    ///Ask if use default scalar LUT
    QMessageBox msgBox;
    msgBox.setText("An auto adjusted bar is about to be created");
    msgBox.setInformativeText("Would you like to customise the bar?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);

    int ret = QMessageBox::No;
    if(!quiet)
      ret = msgBox.exec();

    const float barWidth = 0.1, barHeight = 0.7;

    double range[2];
    imageData->GetScalarRange(range);

    lookupTable->SetRange(range[0], range[1]);

    vtkSmartPointer<vtkLogLookupTable> logLookupTable;
    if(milxQtRenderWindow::logScale)
    {
        printInfo("Detected log scale.");
        logLookupTable = vtkSmartPointer<vtkLogLookupTable>::New();
        logLookupTable->DeepCopy(vtkLookupTable::SafeDownCast(lookupTable));
    }

    if(ret == QMessageBox::Yes && !quiet)
    {
        bool ok1 = false, ok2 = false;

        noOfLabels = QInputDialog::getInt(this, tr("How many labels to show"),
                                          tr("Labels:"), noOfLabels, 0, 99, 1, &ok1);
        title = QInputDialog::getText(this, tr("Title of Bar"),
                                          tr("Title:"), QLineEdit::Normal,
                                          title, &ok2);

        if(!ok1 || !ok2)
            return;

        if(milxQtRenderWindow::logScale)
            milxQtRenderWindow::scale->SetLookupTable(logLookupTable);
        else
            milxQtRenderWindow::scale->SetLookupTable(lookupTable);
        milxQtRenderWindow::scale->SetNumberOfLabels(noOfLabels);

        vtkImageMapToWindowLevelColors *filterColorsOverlay = viewer->GetWindowLevel();
          filterColorsOverlay->SetLookupTable(lookupTable);
          filterColorsOverlay->PassAlphaToOutputOn();
          filterColorsOverlay->Update();

        milxQtRenderWindow::customScalarBar = true;
    }
    else if(quiet && minRange != maxRange)
    {
        printInfo("Using custom scalar range for scalars.");
        if(milxQtRenderWindow::logScale)
            milxQtRenderWindow::scale->SetLookupTable(logLookupTable);
        else
            milxQtRenderWindow::scale->SetLookupTable(lookupTable);
        milxQtRenderWindow::scale->SetNumberOfLabels(noOfLabels);

        vtkImageMapToWindowLevelColors *filterColorsOverlay = viewer->GetWindowLevel();
          filterColorsOverlay->SetLookupTable(lookupTable);
          filterColorsOverlay->PassAlphaToOutputOn();
          filterColorsOverlay->Update();

        milxQtRenderWindow::customScalarBar = true;
    }
    else
    {
        printInfo("Using scalar range from image.");
        if(milxQtRenderWindow::logScale)
            milxQtRenderWindow::scale->SetLookupTable(logLookupTable);
        else
            milxQtRenderWindow::scale->SetLookupTable(lookupTable);
        milxQtRenderWindow::scale->SetNumberOfLabels(3);

        vtkImageMapToWindowLevelColors *filterColorsOverlay = viewer->GetWindowLevel();
          filterColorsOverlay->SetLookupTable(lookupTable);
          filterColorsOverlay->PassAlphaToOutputOn();
          filterColorsOverlay->Update();

        milxQtRenderWindow::customScalarBar = false;
    }

    milxQtRenderWindow::scale->SetTitle(title.toStdString().c_str());
    milxQtRenderWindow::scale->GetLabelTextProperty()->SetFontFamilyToArial();
    milxQtRenderWindow::scale->GetLabelTextProperty()->SetFontSize(8);
    milxQtRenderWindow::scale->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    milxQtRenderWindow::scale->GetPositionCoordinate()->SetValue(.2,.05);
    milxQtRenderWindow::scale->SetWidth( barWidth );
    milxQtRenderWindow::scale->SetHeight( barHeight );
    milxQtRenderWindow::scale->SetPosition( 0.99 - barWidth, 0.1 );
    milxQtRenderWindow::scale->SetLabelFormat("%-#6.3f");
    milxQtRenderWindow::scale->GetTitleTextProperty()->SetFontFamilyToArial();
    milxQtRenderWindow::scale->GetTitleTextProperty()->SetFontSize(8);
    milxQtRenderWindow::scale->GetLabelTextProperty()->SetJustificationToCentered();

    if(milxQtRenderWindow::backgroundAct->isChecked())
    {
        milxQtRenderWindow::scale->GetLabelTextProperty()->SetColor(0, 0, 0);
        milxQtRenderWindow::scale->GetTitleTextProperty()->SetColor(0, 0, 0);
    }

    //Add scale to scale widget
    milxQtRenderWindow::scalarBar->SetInteractor(QVTKWidget::GetInteractor());
    milxQtRenderWindow::scalarBar->SetScalarBarActor(milxQtRenderWindow::scale);
    milxQtRenderWindow::scalarBar->EnabledOn();

    milxQtRenderWindow::scaleBefore = true;
    milxQtRenderWindow::scaleAct->setChecked(true);
}

void milxQtImage::scaleDisplay(const bool forceDisplay)
{
    if(!viewerSetup)
    {
        printError("Image Data not generated. Ignoring scalar bar operation.");
        return;
    }

    if(!lookupTable)
    {
      QMessageBox msgBox;
      msgBox.setText("Currently not using a colour map");
      msgBox.setInformativeText("Please choose a colour map first");
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.exec();
      scaleAct->setChecked(false);
      return;
    }

    if(scaleAct->isChecked() || forceDisplay)
        enableScale(imageData->GetPointData()->GetScalars()->GetName(), forceDisplay);
    else
        disableScale();

    Render();
}

#if VTK_MAJOR_VERSION > 5
void milxQtImage::resliceMode(const bool quietly)
{
    if(resliceAct->isChecked())
        enableResliceMode();
    else
        disableResliceMode();

    if(!quietly)
        emit modified(this);
}
#endif

void milxQtImage::showCrosshair(const bool quietly)
{
    if(cursorAct->isChecked())
        enableCrosshair();
    else
        disableCrosshair();

    if(!quietly)
        emit modified(this);
}

void milxQtImage::setView(int viewMode)
{
    if(!volume)
    {
        printDebug("Volume is 2D. Not changing view to " + QString::number(viewMode));
        viewToXYPlane();
        return;
    }

    milxQtRenderWindow::setView(viewMode);
}

void milxQtImage::viewToXYPlane()
{
    if(viewerSetup)
    {
        viewer->SetSliceOrientationToXY();
        currentView = AXIAL;
    }
}

void milxQtImage::viewToZXPlane()
{
    if(viewerSetup)
    {
        viewer->SetSliceOrientationToXZ();
        currentView = CORONAL;
    }
}

void milxQtImage::viewToZYPlane()
{
    if(viewerSetup)
    {
        viewer->SetSliceOrientationToYZ();
        currentView = SAGITTAL;
    }
}

void milxQtImage::updateLookupTable()
{
    double range[2];
    imageData->GetScalarRange(range);

    lookupTable->SetRange(range[0], range[1]);

    vtkImageMapToWindowLevelColors *filterColorsOverlay = viewer->GetWindowLevel();
        filterColorsOverlay->SetLookupTable(lookupTable);
        filterColorsOverlay->PassAlphaToOutputOn();
        filterColorsOverlay->Update();

    scaleDisplay();
}

void milxQtImage::histogram(int bins, float belowValue, float aboveValue, bool plotHistogram)
{
    printInfo("Computing Histogram of Image");
    updateData(orientAct->isChecked());

    double range[2];
    imageData->GetScalarRange(range);

    int ret = QMessageBox::No;
    if(plotHistogram)
    {
        ///ask user number of bins
        bool ok1 = false, ok2 = false, ok3 = false;
        bins = QInputDialog::getInt(this, tr("Please Provide the number of bins"),
                                          tr("Number of Bins:"), 256, 2, 16384, 1, &ok1);
        belowValue = QInputDialog::getDouble(this, tr("Please Provide the lower level"),
                                         tr("Lower Level:"), range[0], -2147483647, 2147483647, 1, &ok2);
        aboveValue = QInputDialog::getDouble(this, tr("Please Provide the upper level"),
                                         tr("Upper Level:"), range[1], belowValue, 2147483647, 1, &ok3);

        if(!ok1 || !ok2 || !ok3)
            return;

        QMessageBox msgBox;
        msgBox.setText("Ignore Zero Values in Statistics?");
        msgBox.setInformativeText("Do you wish to ignore zero values in Max/Min/Sum/Mean statistics?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        ret = msgBox.exec();
    }
    else
    {
        belowValue = range[0];
        aboveValue = range[1];
    }

    double binSpacing = (aboveValue-belowValue)/bins;
    printInfo("Bin Spacing: " + QString::number(binSpacing));
    hist = vtkSmartPointer<vtkImageAccumulate>::New(); //!< Histogram of the image
    #if VTK_MAJOR_VERSION <= 5
        hist->SetInput(imageData);
    #else
        hist->SetInputData(imageData);
    #endif
        //Only set first values as single component image assumed
        hist->SetComponentExtent(0, bins-1, 0, 0, 0, 0); //bins
        hist->SetComponentOrigin(belowValue, 0, 0); //offset
        hist->SetComponentSpacing(binSpacing, 0, 0); //spacing
        if(ret == QMessageBox::Yes)
        {
            printWarning("Ignore zero values in summary statistics");
            hist->IgnoreZeroOn();
        }
        else
            hist->IgnoreZeroOff();
        linkProgressEventOf(hist);
        hist->Update();
        hist->Print(cout);

    const coordinate meanTuple( hist->GetMean() );
    const coordinate stddevTuple(hist->GetStandardDeviation());
    const coordinate minTuple( hist->GetMin() );
    const coordinate maxTuple( hist->GetMax() );
    if(imageData->GetNumberOfScalarComponents() == 1)
    {
        meanValue = meanTuple[0];
        stddevValue = stddevTuple[0];
        minValue = minTuple[0];
        maxValue = maxTuple[0];
    }
    else
    {
        meanValue = meanTuple.mean();
        stddevValue = stddevTuple.mean();
        minValue = minTuple.mean();
        maxValue = maxTuple.mean();
    }
    printInfo("Mean of the data: " + QString::number(meanValue));
    printInfo("Std Deviation of the data: " + QString::number(stddevValue));
    printInfo("Min/Max of the data: " + QString::number(minValue) + "/" + QString::number(maxValue));

    if(plotHistogram)
    {
    #if VTK_MAJOR_VERSION <= 5
        hist->GetOutput()->Update();
    #endif

        vtkSmartPointer<vtkFloatArray> binsArray = vtkSmartPointer<vtkFloatArray>::New();
            binsArray->SetNumberOfComponents(1);
            binsArray->SetNumberOfTuples(bins);
            binsArray->SetName("Greyscale Bins");
        vtkSmartPointer<vtkIntArray> freqArray = vtkSmartPointer<vtkIntArray>::New();
            freqArray->SetNumberOfComponents(1);
            freqArray->SetNumberOfTuples(bins);
            freqArray->SetName("Frequency");
            freqArray->FillComponent(0, 0);
        int k = 0;
        for(float j = belowValue; j < aboveValue; j += binSpacing, k ++)
        {
            binsArray->SetTuple1(k, j);
            freqArray->SetTuple1(k, hist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(k));
        }

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            table->AddColumn(binsArray);
            table->AddColumn(freqArray);

        emit tableToPlot(table, "Histogram");
    }
}

void milxQtImage::surfacePlot()
{
    int extent[6];
    viewer->GetImageActor()->GetDisplayExtent(extent);
    printDebug("Currently Display Extent: [" +
        QString::number(extent[0]) + ", " + QString::number(extent[1]) + ", " +
        QString::number(extent[2]) + ", " + QString::number(extent[3]) + ", " +
        QString::number(extent[4]) + ", " + QString::number(extent[5]) + "]");

    vtkSmartPointer<vtkImageData> img;
    if(flipped)
    {
        vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
        #if VTK_MAJOR_VERSION <= 5
            imageReorient->SetInput(imageData);
        #else
            imageReorient->SetInputData(imageData);
        #endif
            imageReorient->SetFilteredAxis(1);
            imageReorient->FlipAboutOriginOn();
            linkProgressEventOf(imageReorient);
            imageReorient->Update();
            img = imageReorient->GetOutput();

        int actualExtent[6];
        imageData->GetExtent(actualExtent);

        if(extent[3]-extent[2] == 0) //flip y extent
        {
            extent[2] = actualExtent[3]-extent[2];
            extent[3] = actualExtent[3]-extent[3];
        }
    }
    else
        img = imageData;

    //Pull out the current slice being viewed
    vtkSmartPointer<vtkImageReslice> slice = vtkSmartPointer<vtkImageReslice>::New();
        slice->SetOutputExtent(extent);
    #if VTK_MAJOR_VERSION <= 5
        slice->SetInput(img);
    #else
        slice->SetInputData(img);
    #endif
        slice->SetOutputSpacing(img->GetSpacing()); //needed
        slice->SetOutputOrigin(img->GetOrigin()); //needed
        linkProgressEventOf(slice);
        slice->Update();

    //check slice orientation
    const int zAxis = viewer->GetSliceOrientation();
    int axisToDisplace = 2; //z
    if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_XZ) //y
        axisToDisplace = 1;
    else if(zAxis == vtkImageViewer2::SLICE_ORIENTATION_YZ) //x
        axisToDisplace = 0;

    emit imageToPlot(slice->GetOutput(), axisToDisplace);
}

void milxQtImage::refresh()
{
    viewer->UpdateCursor();
    viewer->GetInteractorStyle()->InvokeEvent(vtkCommand::ResetWindowLevelEvent); //Reset window level as if pressing 'r'
    viewer->Render();
    milxQtRenderWindow::refresh();
}

void milxQtImage::reset()
{
    updateData(orientAct->isChecked());
    viewer->GetRenderer()->ResetCamera(); //Reset window view as if pressing 'Shift+r'
    refresh();
}

void milxQtImage::createMenu(QMenu *menu)
{
    if(!menu)
        return;

    menu->clear();
    menu->addMenu(basicContextMenu());

    if(!extActionsToAdd.empty())
        menu->addSeparator()->setText(tr("Extensions"));
    foreach(QAction *currAct, extActionsToAdd)
    {
        menu->addAction(currAct);
    }

    menu->addSeparator();
    menu->addAction(milxQtRenderWindow::scaleAct);
    menu->addMenu(milxQtRenderWindow::contourMenu);
    menu->addMenu(milxQtRenderWindow::windowPropertiesMenu);

    menu->addSeparator();
    foreach(QAction *currAct, milxQtWindow::actionsToAppend)
    {
      menu->addAction(currAct);
    }
    foreach(QMenu *currMenu, milxQtWindow::menusToAppend)
    {
      menu->addMenu(currMenu);
    }
    menu->addAction(milxQtRenderWindow::refreshAct);
    menu->addAction(milxQtRenderWindow::resetAct);

    //disabling
    milxQtRenderWindow::contourPolyDataAct->setDisabled(!contourAct->isChecked());
    milxQtRenderWindow::contourNodePolyDataAct->setDisabled(!contourAct->isChecked());
    milxQtRenderWindow::contourInitAct->setDisabled(!contourAct->isChecked());
}

void milxQtImage::customOperation()
{
    milxQtImage *currentImgWindow = qobject_cast<milxQtImage *>( windowActionGroup->checkedAction()->parent() );

    if(currentImgWindow == 0) //Not image so...
    {
        milxQtRenderWindow *currentWindow = qobject_cast<milxQtRenderWindow *>( windowActionGroup->checkedAction()->parent() );

        if(currentWindow == 0)
        {
            printError("Error in determining parent of action");
            return;
        }

        vtkActorCollection *actors = currentWindow->GetActors();
        printDebug("Checked Window is: " + currentWindow->strippedNamePrefix());

        size_t n = actors->GetNumberOfItems();

        actors->InitTraversal();
        for(size_t j = 0; j < n; j ++)
        {
            printDebug("Adding Generic Actor into scene.");
            vtkActor *currentActor = actors->GetNextItem();

            //Check if already in view
            foreach(ModelActorItem item, modelActors)
            {
                if(currentActor == item.parentActor)
                    return; //Dont need to do anything.
            }

            ModelActorItem item;
            item.parentActor = currentActor;
            item.modelActor = vtkSmartPointer<vtkActor>::New();
            item.modelActor->ShallowCopy(currentActor);
            item.modelActor->GetProperty()->LightingOn();
            //Fixed discoloration caused by transforming the actor.
            //http://vtk.1045678.n5.nabble.com/diffuse-lighting-after-scaling-with-negative-values-td1240907.html
            milxQtRenderWindow::lightingAct->setChecked(false); //Disable two-sided lighting to display scalars properly on meshes
            milxQtRenderWindow::lighting();
            item.modelActor->Modified();

            Connector->Connect(currentActor,
                               vtkCommand::ModifiedEvent,
                               this,
                               SLOT( updateModelActor(vtkObject*, unsigned long, void*, void*, vtkCommand*) ),
                               NULL, 1.0); //High Priority

            if(transformMatrix)
            {
                printDebug("Transforming Actor to Image Space");
                item.modelActor->SetUserMatrix( transformMatrix );
                //Fix the actor properties, as they get reset to default ones. \todo fix this
            }
            modelActors.append(item);
            AddActor(item.modelActor); //!< transfer model actor to current model

            qApp->processEvents(); ///Keep UI responsive
        }
    }

    Render();
}

void milxQtImage::updateModelActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
{
    vtkSmartPointer<vtkActor> actor = vtkActor::SafeDownCast(obj);

    //Update the actor which is somehere in list
    foreach(ModelActorItem item, modelActors)
    {
        if(actor == item.parentActor)
        {
            item.modelActor->ShallowCopy(actor);
            item.modelActor->Modified();

            break; //Only one of each is allowed
        }
    }

    Render();
}

void milxQtImage::updateDisplay(QPointer<milxQtImage> img)
{
    printDebug("Updating dependent display");
    setDisplayData(img);
    generateImage(true);

    interpolateAct->setChecked(img->isInterpolated());
    if(img->isInterpolated() ^ isInterpolated()) //XOR
        interpolateDisplay(true);

    if(img->isOriented() ^ isOriented()) //XOR
    {
        orientAct->setChecked(img->isOriented());
        applyOrientDisplay(true);
    }

    if(img->isCrosshair() ^ isCrosshair())
    {
        cursorAct->setChecked(img->isCrosshair());
        showCrosshair(true);
    }

    Render();
}

void milxQtImage::createActions()
{
    //Filters
    operateMenu = new QMenu(this);
    operateMenu->setTitle(QApplication::translate("Image", "Operations", 0, QApplication::UnicodeUTF8));
    rescaleAct = new QAction(this);
    rescaleAct->setText(QApplication::translate("Image", "Rescale Intensities", 0, QApplication::UnicodeUTF8));
    rescaleAct->setShortcut(tr("Alt+r"));
    equaliseAct = new QAction(this);
    equaliseAct->setText(QApplication::translate("Image", "Histogram Equalisation", 0, QApplication::UnicodeUTF8));
    equaliseAct->setShortcut(tr("Alt+h"));
    computeContourAct = new QAction(this);
    computeContourAct->setText(QApplication::translate("Image", "Compute Contour", 0, QApplication::UnicodeUTF8));
    computeContourAct->setShortcut(tr("Alt+m"));
    smoothAct = new QAction(this);
    smoothAct->setText(QApplication::translate("Image", "Smooth via Anisotropic Diffusion", 0, QApplication::UnicodeUTF8));
    smoothAct->setShortcut(tr("Alt+s"));
    gaussianAct = new QAction(this);
    gaussianAct->setText(QApplication::translate("Image", "Smooth via Gaussian Convolution", 0, QApplication::UnicodeUTF8));
    gaussianAct->setShortcut(tr("Alt+c"));
    bilateralAct = new QAction(this);
    bilateralAct->setText(QApplication::translate("Image", "Smooth via Bilateral Filter", 0, QApplication::UnicodeUTF8));
    bilateralAct->setShortcut(tr("Alt+b"));
    medianAct = new QAction(this);
    medianAct->setText(QApplication::translate("Image", "Smooth via Median", 0, QApplication::UnicodeUTF8));
    medianAct->setShortcut(tr("Shift+Alt+s"));
    gradMagAct = new QAction(this);
    gradMagAct->setText(QApplication::translate("Image", "Gradient Magnitude", 0, QApplication::UnicodeUTF8));
    gradMagAct->setShortcut(tr("Alt+g"));
    sobelAct = new QAction(this);
    sobelAct->setText(QApplication::translate("Image", "Sobel Edge Detection", 0, QApplication::UnicodeUTF8));
    sobelAct->setShortcut(tr("Alt+e"));
    cannyAct = new QAction(this);
    cannyAct->setText(QApplication::translate("Image", "Canny Edge Detection", 0, QApplication::UnicodeUTF8));
    cannyAct->setShortcut(tr("Shift+Alt+e"));
    laplacianAct = new QAction(this);
    laplacianAct->setText(QApplication::translate("Image", "Apply Laplacian", 0, QApplication::UnicodeUTF8));
    laplacianAct->setShortcut(tr("Alt+l"));
    highPassAct = new QAction(this);
    highPassAct->setText(QApplication::translate("Image", "Butterworth High-Pass Filter", 0, QApplication::UnicodeUTF8));
    highPassAct->setShortcut(tr("Shift+Alt+b"));
    normAct = new QAction(this);
    normAct->setText(QApplication::translate("Image", "Normalize", 0, QApplication::UnicodeUTF8));
    normAct->setShortcut(tr("Alt+n"));
    invertAct = new QAction(this);
    invertAct->setText(QApplication::translate("Image", "Invert Intensities", 0, QApplication::UnicodeUTF8));
    invertAct->setShortcut(tr("Alt+v"));
    relabelAct = new QAction(this);
    relabelAct->setText(QApplication::translate("Image", "Relabel", 0, QApplication::UnicodeUTF8));
    relabelAct->setShortcut(tr("Shift+Alt+l"));
    //Transform
    transformMenu = new QMenu(this);
    transformMenu->setTitle(QApplication::translate("Image", "Transforms", 0, QApplication::UnicodeUTF8));
    matchAct = new QAction(this);
    matchAct->setText(QApplication::translate("Image", "Match Information to ...", 0, QApplication::UnicodeUTF8));
    matchAct->setShortcut(tr("Shift+Alt+i"));
    matchHistAct = new QAction(this);
    matchHistAct->setText(QApplication::translate("Image", "Match Histogram to ...", 0, QApplication::UnicodeUTF8));
    matchHistAct->setShortcut(tr("Shift+Alt+h"));
    resampleSpacingAct = new QAction(this);
    resampleSpacingAct->setText(QApplication::translate("Image", "Resample to spacing ...", 0, QApplication::UnicodeUTF8));
    resampleSpacingAct->setShortcut(tr("Ctrl+Alt+r"));
    resampleAct = new QAction(this);
    resampleAct->setText(QApplication::translate("Image", "Resample Image to ...", 0, QApplication::UnicodeUTF8));
    resampleAct->setShortcut(tr("Shift+Alt+r"));
    resampleLabelAct = new QAction(this);
    resampleLabelAct->setText(QApplication::translate("Image", "Resample as Labelled Image to ...", 0, QApplication::UnicodeUTF8));
    resampleLabelAct->setShortcut(tr("Shift+Alt+l"));
    subsampleAct = new QAction(this);
    subsampleAct->setText(QApplication::translate("Image", "Subsample Image", 0, QApplication::UnicodeUTF8));
    subsampleAct->setShortcut(tr("Shift+Alt+s"));
    transformAct = new QAction(this);
    transformAct->setText(QApplication::translate("Image", "Transform via File ...", 0, QApplication::UnicodeUTF8));
    transformAct->setShortcut(tr("Shift+Alt+t"));
    maskAct = new QAction(this);
    maskAct->setText(QApplication::translate("Image", "Mask Image with ...", 0, QApplication::UnicodeUTF8));
    maskAct->setShortcut(tr("Shift+Alt+m"));
    cropAct = new QAction(this);
    cropAct->setText(QApplication::translate("Image", "Mask and Crop Image with ...", 0, QApplication::UnicodeUTF8));
    cropAct->setShortcut(tr("Shift+Alt+a"));
    checkerAct = new QAction(this);
    checkerAct->setText(QApplication::translate("Image", "Compare as Checkboard to ...", 0, QApplication::UnicodeUTF8));
    checkerAct->setShortcut(tr("Shift+Alt+c"));
    distMapAct = new QAction(this);
    distMapAct->setText(QApplication::translate("Image", "Distance Map", 0, QApplication::UnicodeUTF8));
    distMapAct->setShortcut(tr("Alt+d"));
    flipAct = new QAction(this);
    flipAct->setText(QApplication::translate("Image", "Flip", 0, QApplication::UnicodeUTF8));
    flipAct->setShortcut(tr("Alt+f"));
    //Threshold
    thresholdMenu = new QMenu(this);
    thresholdMenu->setTitle(QApplication::translate("Image", "Thresholds", 0, QApplication::UnicodeUTF8));
    otsuAct = new QAction(this);
    otsuAct->setText(QApplication::translate("Image", "Otsu Threshold", 0, QApplication::UnicodeUTF8));
    otsuAct->setShortcut(tr("Alt+u"));
    otsuMultipleAct = new QAction(this);
    otsuMultipleAct->setText(QApplication::translate("Image", "Otsu Multiple Threshold", 0, QApplication::UnicodeUTF8));
    otsuMultipleAct->setShortcut(tr("Shift+Alt+u"));
    binaryAct = new QAction(this);
    binaryAct->setText(QApplication::translate("Image", "Binary Threshold", 0, QApplication::UnicodeUTF8));
    binaryAct->setShortcut(tr("Shift+Alt+b"));
    bandAct = new QAction(this);
    bandAct->setText(QApplication::translate("Image", "Threshold Outside Band", 0, QApplication::UnicodeUTF8));
    bandAct->setShortcut(tr("Alt+t"));
    aboveAct = new QAction(this);
    aboveAct->setText(QApplication::translate("Image", "Threshold Above", 0, QApplication::UnicodeUTF8));
    aboveAct->setShortcut(tr("Shift+Alt+a"));
    belowAct = new QAction(this);
    belowAct->setText(QApplication::translate("Image", "Threshold Below", 0, QApplication::UnicodeUTF8));
    belowAct->setShortcut(tr("Shift+Alt+b"));
    //Vector imaging
    vectorMenu = new QMenu(this);
    vectorMenu->setTitle(QApplication::translate("Image", "Complex/Vector/4D Imaging", 0, QApplication::UnicodeUTF8));
    vectorMenu->setDisabled(false);
    vectorMagnitudeAct = new QAction(this);
    vectorMagnitudeAct->setText(QApplication::translate("Image", "Display Magnitude", 0, QApplication::UnicodeUTF8));
    vectorMagnitudeAct->setShortcut(tr("Alt+m"));
    vectorComponentAct = new QAction(this);
    vectorComponentAct->setText(QApplication::translate("Image", "Display Component ...", 0, QApplication::UnicodeUTF8));
    vectorComponentAct->setShortcut(tr("Alt+c"));
    pseudoImageAct = new QAction(this);
    pseudoImageAct->setText(QApplication::translate("Image", "Display Pseudo-Image", 0, QApplication::UnicodeUTF8));
    pseudoImageAct->setShortcut(tr("Alt+p"));
    vectorFieldAct = new QAction(this);
    vectorFieldAct->setText(QApplication::translate("Image", "Display Vector/Tensor Field", 0, QApplication::UnicodeUTF8));
    vectorFieldAct->setShortcut(tr("Alt+f"));
    streamLinesAct = new QAction(this);
    streamLinesAct->setText(QApplication::translate("Image", "Display Streamlines from Slice", 0, QApplication::UnicodeUTF8));
    streamLinesAct->setShortcut(tr("Shift+Alt+s"));

    //Display
    levelAct = new QAction(this);
    levelAct->setText(QApplication::translate("Image", "Auto-Level Display", 0, QApplication::UnicodeUTF8));
    levelAct->setShortcut(tr("Alt+o"));
    overlayAct = new QAction(this);
    overlayAct->setText(QApplication::translate("Image", "Overlay Labelled Image from ...", 0, QApplication::UnicodeUTF8));
    overlayAct->setShortcut(tr("Alt+o"));
    overlayContourAct = new QAction(this);
    overlayContourAct->setText(QApplication::translate("Image", "Overlay Labelled Image as Contour from ...", 0, QApplication::UnicodeUTF8));
    overlayContourAct->setShortcut(tr("Shift+Alt+o"));
    blendAct = new QAction(this);
    blendAct->setText(QApplication::translate("Image", "Blend Image with ...", 0, QApplication::UnicodeUTF8));
    blendAct->setShortcut(tr("Shift+Alt+b"));
    volRenderAct = new QAction(this);
    volRenderAct->setText(QApplication::translate("Image", "Display as Volume Rendering", 0, QApplication::UnicodeUTF8));
    volRenderAct->setShortcut(tr("Alt+v"));
    histogramAct = new QAction(this);
    histogramAct->setText(QApplication::translate("Image", "Display Histogram", 0, QApplication::UnicodeUTF8));
    histogramAct->setShortcut(tr("Alt+h"));
    surfacePlotAct = new QAction(this);
    surfacePlotAct->setText(QApplication::translate("Image", "Display Slice Surface Plot", 0, QApplication::UnicodeUTF8));
    surfacePlotAct->setShortcut(tr("Alt+s"));
    surfaceAct = new QAction(this);
    surfaceAct->setText(QApplication::translate("Image", "Display Iso-surface", 0, QApplication::UnicodeUTF8));
    surfaceAct->setShortcut(tr("Shift+Alt+s"));
    polyDataAct = new QAction(this);
    polyDataAct->setText(QApplication::translate("Image", "Generate Polygonal Data", 0, QApplication::UnicodeUTF8));
    polyDataAct->setShortcut(tr("Shift+Alt+p"));
    polyDataAct->setDisabled(true); //!< \todo Disabled because feature is broken, fix
    infoAct = new QAction(this);
    infoAct->setText(QApplication::translate("Image", "Image Information", 0, QApplication::UnicodeUTF8));
    infoAct->setShortcut(tr("Alt+i"));
    interpolateAct = new QAction(this);
    interpolateAct->setText(QApplication::translate("Image", "Interpolation", 0, QApplication::UnicodeUTF8));
    interpolateAct->setShortcut(tr("Shift+Alt+i"));
    interpolateAct->setCheckable(true);
    interpolateAct->setChecked(true);
    orientAct = new QAction(this);
    orientAct->setText(QApplication::translate("Image", "Apply Orientation", 0, QApplication::UnicodeUTF8));
    orientAct->setShortcut(tr("Shift+Alt+o"));
    orientAct->setCheckable(true);
    orientAct->setChecked(true);
    resliceAct = new QAction(this);
    resliceAct->setText(QApplication::translate("Image", "3D Slice View Mode", 0, QApplication::UnicodeUTF8));
    resliceAct->setShortcut(tr("Shift+Ctrl+s"));
    resliceAct->setCheckable(true);
    resliceAct->setChecked(false);
#if VTK_MAJOR_VERSION <= 5
    resliceAct->setDisabled(true);
#endif
    cursorAct = new QAction(this);
    cursorAct->setText(QApplication::translate("Image", "Show Cursor", 0, QApplication::UnicodeUTF8));
    cursorAct->setShortcut(tr("Shift+Alt+c"));
    cursorAct->setCheckable(true);
    cursorAct->setChecked(false);

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
    cropAct->setDisabled(false);
    computeContourAct->setDisabled(false);
    overlayAct->setDisabled(false);
    overlayContourAct->setDisabled(false);
#else
    cropAct->setDisabled(true);
    computeContourAct->setDisabled(true);
    overlayAct->setDisabled(true);
    overlayContourAct->setDisabled(true);
#endif

    milxQtRenderWindow::actionDefault->setChecked(true);
}

void milxQtImage::createConnections()
{
    connect(this, SIGNAL(mouseEvent(QMouseEvent*)), this, SLOT(userEvent(QMouseEvent*)));

    //Filter
    connect(rescaleAct, SIGNAL(triggered()), this, SLOT(rescale()));
    connect(equaliseAct, SIGNAL(triggered()), this, SLOT(histogramEqualisation()));
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
    connect(computeContourAct, SIGNAL(triggered()), this, SLOT(computeContour()));
#endif
    connect(smoothAct, SIGNAL(triggered()), this, SLOT(anisotropicDiffusion()));
    connect(gaussianAct, SIGNAL(triggered()), this, SLOT(gaussianSmooth()));
    connect(bilateralAct, SIGNAL(triggered()), this, SLOT(bilateral()));
    connect(medianAct, SIGNAL(triggered()), this, SLOT(median()));
    connect(gradMagAct, SIGNAL(triggered()), this, SLOT(gradientMagnitude()));
    connect(sobelAct, SIGNAL(triggered()), this, SLOT(sobelEdges()));
    connect(cannyAct, SIGNAL(triggered()), this, SLOT(cannyEdges()));
    connect(laplacianAct, SIGNAL(triggered()), this, SLOT(laplacian()));
    connect(highPassAct, SIGNAL(triggered()), this, SLOT(highpass()));
    connect(normAct, SIGNAL(triggered()), this, SLOT(normalize()));
    connect(invertAct, SIGNAL(triggered()), this, SLOT(invertIntensity()));
    connect(relabelAct, SIGNAL(triggered()), this, SLOT(relabel()));

    connect(matchAct, SIGNAL(triggered()), this, SLOT(matchInfo()));
    connect(matchHistAct, SIGNAL(triggered()), this, SLOT(matchHistogram()));
    connect(resampleSpacingAct, SIGNAL(triggered()), this, SLOT(resize()));
    connect(resampleAct, SIGNAL(triggered()), this, SLOT(resample()));
    connect(resampleLabelAct, SIGNAL(triggered()), this, SLOT(resampleLabel()));
    connect(subsampleAct, SIGNAL(triggered()), this, SLOT(subsample()));
    connect(transformAct, SIGNAL(triggered()), this, SLOT(transform()));
    connect(maskAct, SIGNAL(triggered()), this, SLOT(mask()));
    connect(checkerAct, SIGNAL(triggered()), this, SLOT(checkerBoard()));
    connect(distMapAct, SIGNAL(triggered()), this, SLOT(distanceMap()));
    connect(flipAct, SIGNAL(triggered()), this, SLOT(flip()));
    connect(surfaceAct, SIGNAL(triggered()), this, SLOT(surface()));
    connect(polyDataAct, SIGNAL(triggered()), this, SLOT(polyData()));

    connect(aboveAct, SIGNAL(triggered()), this, SLOT(thresholdAbove()));
    connect(belowAct, SIGNAL(triggered()), this, SLOT(thresholdBelow()));
    connect(bandAct, SIGNAL(triggered()), this, SLOT(threshold()));
    connect(otsuAct, SIGNAL(triggered()), this, SLOT(otsu()));
    connect(otsuMultipleAct, SIGNAL(triggered()), this, SLOT(otsuMultiple()));
    connect(binaryAct, SIGNAL(triggered()), this, SLOT(binaryThreshold()));

    connect(vectorMagnitudeAct, SIGNAL(triggered()), this, SLOT(magnitude()));
    connect(vectorComponentAct, SIGNAL(triggered()), this, SLOT(component()));
    connect(pseudoImageAct, SIGNAL(triggered()), this, SLOT(pseudoImage()));
    connect(vectorFieldAct, SIGNAL(triggered()), this, SLOT(vectorField()));
    connect(streamLinesAct, SIGNAL(triggered()), this, SLOT(streamLines()));

    //Display
    connect(levelAct, SIGNAL(triggered()), this, SLOT(autoLevel()));
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
    connect(cropAct, SIGNAL(triggered()), this, SLOT(crop()));
    connect(overlayAct, SIGNAL(triggered()), this, SLOT(overlay()));
    connect(overlayContourAct, SIGNAL(triggered()), this, SLOT(overlayContour()));
#endif
    connect(blendAct, SIGNAL(triggered()), this, SLOT(blend()));
    connect(volRenderAct, SIGNAL(triggered()), this, SLOT(volumeRendering()));
    connect(histogramAct, SIGNAL(triggered()), this, SLOT(histogram()));
    connect(surfacePlotAct, SIGNAL(triggered()), this, SLOT(surfacePlot()));
    connect(infoAct, SIGNAL(triggered()), this, SLOT(imageInformation()));
    connect(interpolateAct, SIGNAL(triggered()), this, SLOT(interpolateDisplay()));
    connect(orientAct, SIGNAL(triggered()), this, SLOT(applyOrientDisplay()));
#if VTK_MAJOR_VERSION > 5
    connect(resliceAct, SIGNAL(triggered()), this, SLOT(resliceMode()));
#endif
    connect(cursorAct, SIGNAL(triggered()), this, SLOT(showCrosshair()));
    connect(milxQtRenderWindow::refreshAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(milxQtRenderWindow::resetAct, SIGNAL(triggered()), this, SLOT(reset()));

    milxQtWindow::createConnections(); //consume right click events etc.
}

void milxQtImage::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    createMenu(contextMenu);

    contextMenu->exec(currentEvent->globalPos());
}

QMenu* milxQtImage::basicContextMenu()
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection
    contextMenu->setTitle(QApplication::translate("MainWindow", "Imaging", 0, QApplication::UnicodeUTF8));

    foreach(QAction *currAct, milxQtWindow::actionsToAdd)
    {
        contextMenu->addAction(currAct);
    }
    foreach(QMenu *currMenu, milxQtWindow::menusToAdd)
    {
        contextMenu->addMenu(currMenu);
    }
    contextMenu->addMenu(operationsMenu());
    contextMenu->addMenu(thresholdsMenu());
    contextMenu->addMenu(transformsMenu());
    contextMenu->addMenu(vectorsMenu());

    contextMenu->addSeparator()->setText(tr("Display"));
    contextMenu->addAction(levelAct);
    contextMenu->addAction(overlayAct);
    contextMenu->addAction(overlayContourAct);
    contextMenu->addAction(blendAct);
    contextMenu->addAction(checkerAct);
    contextMenu->addAction(volRenderAct);
    contextMenu->addAction(histogramAct);
    contextMenu->addAction(surfacePlotAct);
    contextMenu->addAction(surfaceAct);
//    contextMenu->addAction(polyDataAct);
    contextMenu->addAction(infoAct);
    contextMenu->addAction(interpolateAct);
    contextMenu->addAction(orientAct);
    contextMenu->addAction(resliceAct);
    contextMenu->addAction(cursorAct);
    contextMenu->addAction(milxQtRenderWindow::humanAct);
    ///Change View of Volume
    contextMenu->addMenu(milxQtRenderWindow::viewMenu);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewXY);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZX);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZY);
    milxQtRenderWindow::viewGroup->setDisabled(!volume);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewFileAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewFileAct);
    contextMenu->addMenu(milxQtRenderWindow::colourMapMenu);

    milxQtRenderWindow::enableActionBasedOnView();

    return contextMenu;
}

QMenu* milxQtImage::operationsMenu()
{
    operateMenu->addAction(rescaleAct);
    operateMenu->addAction(equaliseAct);
    operateMenu->addAction(computeContourAct);
    operateMenu->addAction(smoothAct);
    operateMenu->addAction(gaussianAct);
    operateMenu->addAction(bilateralAct);
    operateMenu->addAction(medianAct);
    operateMenu->addAction(gradMagAct);
    operateMenu->addAction(sobelAct);
    operateMenu->addAction(cannyAct);
    operateMenu->addAction(laplacianAct);
    operateMenu->addAction(highPassAct);
    operateMenu->addAction(normAct);
    operateMenu->addAction(distMapAct);
    operateMenu->addAction(matchHistAct);
    operateMenu->addAction(invertAct);
    operateMenu->addAction(relabelAct);

    return operateMenu;
}

QMenu* milxQtImage::thresholdsMenu()
{
    thresholdMenu->addAction(otsuAct);
    thresholdMenu->addAction(otsuMultipleAct);
    thresholdMenu->addAction(binaryAct);
    thresholdMenu->addAction(bandAct);
    thresholdMenu->addAction(aboveAct);
    thresholdMenu->addAction(belowAct);

    return thresholdMenu;
}

QMenu* milxQtImage::transformsMenu()
{
    transformMenu->addAction(matchAct);
    transformMenu->addAction(resampleSpacingAct);
    transformMenu->addAction(resampleAct);
    transformMenu->addAction(resampleLabelAct);
    transformMenu->addAction(subsampleAct);
    transformMenu->addAction(transformAct);
    transformMenu->addAction(maskAct);
    transformMenu->addAction(cropAct);
    transformMenu->addAction(flipAct);

    return transformMenu;
}

QMenu* milxQtImage::vectorsMenu()
{
    vectorMenu->addAction(vectorMagnitudeAct);
    vectorMenu->addAction(vectorComponentAct);
    vectorMenu->addAction(pseudoImageAct);
    vectorMenu->addAction(vectorFieldAct);
    vectorMenu->addAction(streamLinesAct);

    //Disable enable based on flags
    vectorMenu->setEnabled(vectorised);
    pseudoImageAct->setDisabled(true); //not working yet

    return vectorMenu;
}

void milxQtImage::dropEvent(QDropEvent *currentEvent)
{
    if(currentEvent->source() == this) //disallow self drops
        return;

    milxQtImage *sourceImage = qobject_cast<milxQtImage *>(currentEvent->source());
    if(currentEvent->spontaneous())
        currentEvent->ignore(); //dont accept drop

    QList<QUrl> urlsList = currentEvent->mimeData()->urls();
    QString tmp, typeString = currentEvent->mimeData()->text();

    for(int j = 0; j < urlsList.size(); j ++)
    {
        if(urlsList[j].isValid())
        {
#ifdef Q_WS_WIN
            tmp = urlsList[j].path().remove(0,1); //!< Remove leading forward slash
#else
            tmp = urlsList[j].path();
#endif
            printInfo("Dropped Path into Image: " + tmp);

            if(sourceImage)
            {
                currentEvent->acceptProposedAction();
                if(currentEvent->proposedAction() == Qt::MoveAction)
//                    blend(tmp);
                    blend(sourceImage);
            #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
                if(currentEvent->proposedAction() == Qt::LinkAction)
                    overlay(tmp);
            #endif
                else
//                    checkerBoard(tmp);
                    checkerBoard(sourceImage);
            }
            Render();
        }
    }
}

void milxQtImage::SetupWidgets(vtkRenderWindowInteractor *interactor)
{
    printInfo("Setting up image for Annotation");
    interactor->Initialize();

    double imageBounds[6], imgCenter[3];
    imageData->GetBounds(imageBounds);
    imageData->GetOrigin(imgCenter);
    double length = imageBounds[1]-imageBounds[0]/10;

    //line widget setup
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

    //sphere widget setup
    sphereRep = vtkSmartPointer<vtkSphereRepresentation>::New();
        sphereRep->SetRepresentationToWireframe();
        sphereRep->SetThetaResolution(16);
        sphereRep->SetPhiResolution(16);
        sphereRep->SetRadius(length);
        sphereRep->SetCenter(imgCenter);
        sphereRep->HandleTextOn();
        sphereRep->RadialLineOn();
        sphereRep->HandleVisibilityOn();
    sphereWidget = vtkSmartPointer<vtkSphereWidget2>::New();
        sphereWidget->SetInteractor(interactor);
//        sphereWidget->CreateDefaultRepresentation();
        sphereWidget->SetRepresentation(sphereRep);
        sphereWidget->ScalingEnabledOn();

    //Change colour
    sphereRep->GetSphereProperty()->SetColor(0.9, 0.9, 1.0);

    vtkSmartPointer<vtkPolyData> sphere = vtkSmartPointer<vtkPolyData>::New();
        sphereRep->GetPolyData(sphere);

    //box widget
    boxWidget = vtkSmartPointer<vtkBoxWidget2>::New();
        boxWidget->SetInteractor(interactor);
        boxWidget->RotationEnabledOn();
        //boxWidget->CreateDefaultRepresentation();

    vtkSmartPointer<vtkBoxRepresentation> boxRepresentation = vtkSmartPointer<vtkBoxRepresentation>::New();
        boxRepresentation->SetPlaceFactor( 1 );
        boxRepresentation->PlaceWidget(sphere->GetBounds());
        boxWidget->SetRepresentation(boxRepresentation);

    //plane widget setup
    double actorBounds[6], center[3];
//    viewer->GetImageActor()->GetBounds(actorBounds);
    sphere->GetBounds(actorBounds);
    center[0] = (actorBounds[0] + actorBounds[1])/2.0;
    center[1] = (actorBounds[2] + actorBounds[3])/2.0;
    center[2] = (actorBounds[4] + actorBounds[5])/2.0;
    planeWidget = vtkSmartPointer<vtkPlaneWidget>::New();
          planeWidget->SetInteractor(interactor);
          planeWidget->SetPlaceFactor( 1 );
//          planeWidget->SetRepresentationToOutline();
//          planeWidget->SetHandleSize(planeWidget->GetHandleSize()*4);
//          planeWidget->SetHandleSize((actorBounds[1]-actorBounds[0])/10.0);
//          planeWidget->NormalToXAxisOn();
//          planeWidget->NormalToYAxisOn();
          planeWidget->NormalToZAxisOn();
//          planeWidget->PlaceWidget(actorBounds[0], actorBounds[1], actorBounds[2], actorBounds[3], actorBounds[4], actorBounds[5]);
//          planeWidget->PlaceWidget(viewer->GetImageActor()->GetBounds());
//          planeWidget->PlaceWidget(sphere->GetBounds());
          //X
//          planeWidget->SetOrigin(center[0],actorBounds[2],actorBounds[4]);
//          planeWidget->SetPoint1(center[0],actorBounds[3],actorBounds[4]);
//          planeWidget->SetPoint2(center[0],actorBounds[2],actorBounds[5]);
          //Y
//          planeWidget->SetOrigin(actorBounds[0],center[1],actorBounds[4]);
//          planeWidget->SetPoint1(actorBounds[1],center[1],actorBounds[4]);
//          planeWidget->SetPoint2(actorBounds[0],center[1],actorBounds[5]);
          //Z
          planeWidget->SetOrigin(actorBounds[0],actorBounds[2],center[2]);
          planeWidget->SetPoint1(actorBounds[1],actorBounds[2],center[2]);
          planeWidget->SetPoint2(actorBounds[0],actorBounds[3],center[2]);
          planeWidget->SetResolution(2);
//          planeWidget->UpdatePlacement();
}

QString milxQtImage::getOpenFilename(const QString labelForDialog, QString exts)
{
    if(exts.isEmpty())
        exts = QString(openMedImageExts.c_str()) + ";;" + QString(openImageExts.c_str()) + ";;" + QString(openOtherImageExts.c_str());
    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString();
    QFileDialog *fileOpener = new QFileDialog(this);
    QString filename = fileOpener->getOpenFileName(this,
                                           tr(labelForDialog.toStdString().c_str()),
                                           path,
                                           tr(exts.toStdString().c_str()) );

    return filename;
}
