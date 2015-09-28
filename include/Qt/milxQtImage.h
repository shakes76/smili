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
#ifndef MILXQTIMAGE_H
#define MILXQTIMAGE_H

#include <QStatusBar>

//VTK Headers
#include <vtkSmartPointer.h>
#include <itkCommand.h>
#include <vtkImageData.h>
#include <vtkMatrix4x4.h>
#include <vtkImageAppend.h>
#include <vtkImagePermute.h>
#include <vtkImageMapToWindowLevelColors.h>
//ITK Imaging
#include <itkRGBPixel.h>
#include <itkImage.h>
#include <itkVectorImage.h>
//VXL
#include <vnl/vnl_matrix.h>

#include "milxQtRenderWindow.h"
#include "vtkImageViewer3.h" //smili vtk-ext
#include "milxImage.h"

//Typedefs
typedef unsigned char charPixelType;
typedef itk::Image<charPixelType, milx::imgDimension> charImageType;
typedef itk::RGBPixel<unsigned char> rgbPixelType;
typedef itk::Image<rgbPixelType, milx::imgDimension> rgbImageType;
typedef float floatPixelType;
typedef itk::Image<floatPixelType, milx::imgDimension> floatImageType;
typedef itk::VectorImage<floatPixelType, milx::imgDimension> vectorImageType;

//Struct for image actors
struct ModelActorItem
{
    vtkSmartPointer<vtkActor> modelActor; //image actor to display in model view
    vtkSmartPointer<vtkActor> parentActor; //The parent actor it's connected to
};

/**
    \class itkEventQtObserver
    \brief This class can be attached as an observer to trigger Qt to update the event loop

    This keeps the application responsive.
*/
class itkEventQtObserver : public itk::Command
{
public:
    itkNewMacro( itkEventQtObserver );

public:

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
        qApp->processEvents(); //Update Qt event loop to keep application responsive
    }
};

/*!
    \class milxQtImage
    \brief This class represents the MILX Qt Image Display object using VTK.
    \author Shekhar S. Chandra, 2013

    The class displays image values using OpenGL via the VTK library. itk::Image or vtkImageData can be passed to this class and
    a number of processes generated from the data.

    The rendering is encapsulated within a QVTK widget.

    Controls:
        - Left Mouse (hold): adjust the colour level.
        - Wheel: zoom the camera.
        - Middle Mouse (hold): translate image.
        - Right Mouse: context menu and zoom toggle.
        - Keypress f: fly to the picked point.
        - Keypress Shift+r: reset the camera view along the current view direction. Centers the actors and moves the camera so that all actors are visible.
        - Keypress r: reset the colour level to defaults.

    Usage Example:
    \code
    milxQtImage *image = new milxQtImage;

    image->setName(imgName);
        image->setData(newImg);
        image->generateImage();
        image->show();
    \endcode
*/
class MILXQT_EXPORT milxQtImage : public milxQtRenderWindow
{
    Q_OBJECT

public:
    /*!
        \fn milxQtImage::milxQtImage(QWidget *parent = 0, bool contextSystem = true)
        \brief The standard constructor
    */
    milxQtImage(QWidget *theParent = 0, bool contextSystem = true);
    /*!
        \fn milxQtImage::~milxQtImage()
        \brief The standard destructor
    */
    virtual ~milxQtImage();

    /*!
        \brief ITK interfacing member for exceptions.
    */
    virtual inline const char * GetNameOfClass() const
    {
      return "milxQtImage";
    }

    //Data Members
    /*!
        \fn milxQtImage::strippedNamePrefix()
        \brief Returns the stripped (path removed) name of the data with "Image" prefix.
    */
    inline QString strippedNamePrefix()
    {
        return prefix + QFileInfo(name).fileName();
    }

    /*!
        \fn milxQtImage::SetInput(vtkSmartPointer<vtkImageData> newImg)
        \brief VTK interface function: Assigns the image data internally. Same as setData() function.
    */
    inline void SetInput(vtkSmartPointer<vtkImageData> newImg)
    {
        setData(newImg);
    }
    /*!
        \fn milxQtImage::SetInput(charImageType::Pointer newImg, const bool flipY = true)
        \brief ITK interface function: Assigns the image data internally. Same as setData() function.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    inline void SetInput(charImageType::Pointer newImg, const bool flipY = true)
    {
        setData(newImg, flipY);
    }
    /*!
        \fn milxQtImage::SetInput(rgbImageType::Pointer newImg, const bool flipY = true)
        \brief ITK interface function: Assigns the image data internally. Same as setData() function.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    inline void SetInput(rgbImageType::Pointer newImg, const bool flipY = true)
    {
        setData(newImg, flipY);
    }
    /*!
        \fn milxQtImage::SetInput(floatImageType::Pointer newImg, const bool flipY = true)
        \brief ITK interface function: Assigns the image data internally. Same as setData() function.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    inline void SetInput(floatImageType::Pointer newImg, const bool flipY = true)
    {
        setData(newImg, flipY);
    }
    /*!
        \fn milxQtImage::SetInput(vectorImageType::Pointer newImg, const bool flipY = true)
        \brief ITK interface function: Assigns the image data internally. Same as setData() function.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    inline void SetInput(vectorImageType::Pointer newImg, const bool flipY = true)
    {
        setData(newImg, flipY);
    }
    /*!
		\fn milxQtImage::SetTransform(vtkSmartPointer<vtkTransform> transform)
		\brief Sets the transform for the image that will be generated. Must pass a vtkTransform objects, which is easy to use.

		Note that this does not change the internal ITK data inside.
    */
    void SetTransform(vtkSmartPointer<vtkTransform> transform);

    /*!
        \fn milxQtImage::setData(QPointer<milxQtImage> newImg, const bool forceDeepCopy = false)
        \brief Assigns the milxQtImage data to current image. You will need to call generate image after this.

        Image provided is generally deep copied except for vector images. Use forceDeepCopy to always deep copy regardless of image type.
    */
    void setData(QPointer<milxQtImage> newImg, const bool forceDeepCopy = false);
    /*!
        \fn milxQtImage::setData(vtkSmartPointer<vtkImageData> newImg)
        \brief Assigns the VTK array data to image. You will need to call generate image after this.
    */
    void setData(vtkSmartPointer<vtkImageData> newImg);
    /*!
        \fn milxQtImage::setData(charImageType::Pointer newImg, const bool flipY = true)
        \brief Assigns the ITK image data to image. You will need to call generate image after this.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    void setData(charImageType::Pointer newImg, const bool flipY = true);
    /*!
        \fn milxQtImage::setData(rgbImageType::Pointer newImg, const bool flipY = true)
        \brief Assigns the ITK image data to image. You will need to call generate image after this.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    void setData(rgbImageType::Pointer newImg, const bool flipY = true);
    /*!
        \fn milxQtImage::setData(floatImageType::Pointer newImg, const bool flipY = true)
        \brief Assigns the ITK image data to image. You will need to call generate image after this.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
    */
    void setData(floatImageType::Pointer newImg, const bool flipY = true);
    /*!
        \fn milxQtImage::setData(vectorImageType::Pointer newImg, const bool flipY = true, const bool deepCopy = false)
        \brief Assigns the ITK image data to image. You will need to call generate image after this.

        ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
        Note that the data itself is not flipped, only the display.
        Note also that unlike the other setData() members, this member does not copy the data since vector image tend to be very large unless deepCopy flag is set.
    */
    void setData(vectorImageType::Pointer newImg, const bool flipY = true, const bool deepCopy = false);
    /*!
        \fn milxQtImage::setData(vnl_matrix<double> &newData)
        \brief Assigns the VNL Matrix to image. You will need to call generate image after this.
    */
    void setData(vnl_matrix<double> &newData);
    /*!
        \fn milxQtImage::setData(const unsigned slice, vnl_matrix<double> &newData)
        \brief Assigns the VNL Matrix to a slice of the image. You will need to call generate image after this. Possibly BROKEN
    */
    void setData(const unsigned slice, vnl_matrix<double> &newData);
    /*!
      \fn milxQtImage::setDisplayData(QPointer<milxQtImage> newImg)
      \brief Shares the milxQtImage data to current image. You will need to call generate image after this.

      Useful when attempting to share image data just for display etc.
    */
    void setDisplayData(QPointer<milxQtImage> newImg);
    /*!
      \fn milxQtImage::setSharedData(QPointer<milxQtImage> newImg)
      \brief Shares the ITK image data to image (same as setDisplayData()). You will need to call generate image after this.
    */
    inline void setSharedData(QPointer<milxQtImage> newImg)
    {	setDisplayData(newImg);	}
    /*!
      \fn milxQtImage::setDisplayData(charImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image. You will need to call generate image after this.

      Useful when attempting to share image data just for display etc.

      ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
      Note that the data itself is not flipped, only the display.
    */
    void setDisplayData(charImageType::Pointer newImg, const bool flipY = true);
    /*!
      \fn milxQtImage::setSharedData(charImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image (same as setDisplayData()). You will need to call generate image after this.
    */
    inline void setSharedData(charImageType::Pointer newImg, const bool flipY = true)
    {	setDisplayData(newImg, flipY);	}
    /*!
      \fn milxQtImage::setDisplayData(rgbImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image. You will need to call generate image after this.

      Useful when attempting to share image data just for display etc.

      ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
      Note that the data itself is not flipped, only the display.
    */
    void setDisplayData(rgbImageType::Pointer newImg, const bool flipY = true);
    /*!
      \fn milxQtImage::setSharedData(rgbImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image (same as setDisplayData()). You will need to call generate image after this.
    */
    inline void setSharedData(rgbImageType::Pointer newImg, const bool flipY = true)
    {	setDisplayData(newImg, flipY);	}
    /*!
      \fn milxQtImage::setDisplayData(floatImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image. You will need to call generate image after this.

      Useful when attempting to share image data just for display etc.

      ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
      Note that the data itself is not flipped, only the display.
    */
    void setDisplayData(floatImageType::Pointer newImg, const bool flipY = true);
    /*!
      \fn milxQtImage::setSharedData(floatImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image (same as setDisplayData()). You will need to call generate image after this.
    */
    inline void setSharedData(floatImageType::Pointer newImg, const bool flipY = true)
    {	setDisplayData(newImg, flipY);	}
    /*!
      \fn milxQtImage::setDisplayData(vectorImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image. You will need to call generate image after this.

      Useful when attempting to share image data just for display etc.

      ITK has different orientation for images, so the flipY flag (not done by default) can be used to flip the image appropriately.
      Note that the data itself is not flipped, only the display.
    */
    void setDisplayData(vectorImageType::Pointer newImg, const bool flipY = true);
    /*!
      \fn milxQtImage::setSharedData(vectorImageType::Pointer newImg, const bool flipY = true)
      \brief Shares the ITK image data to image (same as setDisplayData()). You will need to call generate image after this.
    */
    inline void setSharedData(vectorImageType::Pointer newImg, const bool flipY = true)
    {	setDisplayData(newImg, flipY);	}

    /*!
        \fn milxQtImage::GetOutput()
        \brief Returns the image data object (ImageData) used internally VTK style.
    */
    inline vtkSmartPointer<vtkImageData> GetOutput()
    {
        return imageData;
    }
    /*!
        \fn milxQtImage::GetCharImage()
        \brief Returns the internal unsigned char image data.
    */
    inline charImageType::Pointer GetCharImage()
    {
        return imageChar;
    }
    /*!
        \fn milxQtImage::GetRGBImage()
        \brief Returns the internal RGB image data.
    */
    inline rgbImageType::Pointer GetRGBImage()
    {
        return imageRGB;
    }
    /*!
        \fn milxQtImage::GetFloatImage()
        \brief Returns the internal float image data.
    */
    inline floatImageType::Pointer GetFloatImage()
    {
        return imageFloat;
    }
    /*!
        \fn milxQtImage::GetVectorImage()
        \brief Returns the internal vector image data. Unlike the other Get*Image() members, the return value can be NULL.

        User needs to check if the return value is NULL or not.
    */
    inline vectorImageType::Pointer GetVectorImage()
    {
        return imageVector;
    }
#if VTK_MAJOR_VERSION <= 5
    /*!
          \fn milxQtImage::GetOutputPort()
          \brief Returns the image data object (ImageData) producer port used internally VTK style.
    */
    inline vtkAlgorithmOutput* GetOutputPort()
    {   return imageData->GetProducerPort();   }
#endif
    /*!
        \fn milxQtImage::GetImageActor()
        \brief Returns the internal image actor used for display.
    */
    inline vtkImageActor* GetImageActor()
    {
        return viewer->GetImageActor();
    }
    /*!
        \fn milxQtImage::GetDisplayExtent()
        \brief Returns the current display extent of the image data, i.e the current slice dimensions/extents.
    */
    inline int* GetDisplayExtent()
    {
        return GetImageActor()->GetDisplayExtent();
    }
    /*!
      \brief Get the data, return the vtkImageData object downcast to vtkDataSet, useful for getting scalar range etc.
    */
    inline virtual vtkDataSet* GetDataSet()
    {
        return imageData;
    }
    /*!
      \brief Get the interactor associated with the view rendering
    */
    inline virtual vtkRenderWindowInteractor* GetVTKInteractor()
    {
        return viewer->GetRenderWindow()->GetInteractor();
    }

    //Flags
    /*!
        \fn milxQtImage::set8BitImage()
        \brief Sets image as an 8-bit (unsigned char) image. Relevant only for when creating images, otherwise automatically set.
    */
    inline void set8BitImage()
    {
        eightbit = true;
        rgb = false;
        vectorised = false;
    }
    /*!
        \fn milxQtImage::is8BitImage()
        \brief Returns true if image is an 8-bit (unsigned char) image
    */
    inline bool is8BitImage()
    {
        return eightbit;
    }
    /*!
        \fn milxQtImage::setRGBImage()
        \brief Sets image as an RGB (3-vector unsigned char image) image. Relevant only for when creating images, otherwise automatically set.
    */
    inline void setRGBImage()
    {
        eightbit = false;
        rgb = true;
        vectorised = false;
    }
    /*!
        \fn milxQtImage::isRGBImage()
        \brief Returns true if image is an RGB (3-vector unsigned char image) image
    */
    inline bool isRGBImage()
    {
        return rgb;
    }
    /*!
        \fn milxQtImage::setFloatingPointImage()
        \brief Sets image is a floating point (float) image. Relevant only for when creating images, otherwise automatically set.
    */
    inline void setFloatingPointImage()
    {
        eightbit = false;
        rgb = false;
        vectorised = false;
    }
    /*!
        \fn milxQtImage::isFloatingPointImage()
        \brief Returns true if image is a floating point (float) image
    */
    inline bool isFloatingPointImage()
    {
        return (!eightbit && !rgb && !vectorised);
    }
    /*!
        \fn milxQtImage::setVectorImage()
        \brief Sets image as a vector image. Relevant only for when creating images, otherwise automatically set.
    */
    inline void setVectorImage()
    {
        eightbit = false;
        rgb = false;
        vectorised = true;
    }
    /*!
        \fn milxQtImage::isVectorImage()
        \brief Returns true if image is a vector image
    */
    inline bool isVectorImage()
    {
        return vectorised;
    }
    /*!
        \fn milxQtImage::isVTKImage()
        \brief Returns true if image is a VTK-type image
    */
    inline bool isVTKImage()
    {
        return usingVTKImage;
    }

    /*!
        \fn milxQtImage::isDisplayFlipped()
        \brief Returns true if image was flipped internally for display purposes only.

        Note: Data is never flipped unless the flip() member is called manually by the user.
    */
    inline bool isDisplayFlipped()
    {
        return flipped;
    }
    /*!
        \fn milxQtImage::isInterpolated()
        \brief Returns true if image display is interpolated
    */
    inline bool isInterpolated()
    {
        if(viewer->GetImageActor()->GetInterpolate())
            return true;
        else
            return false;
    }
    /*!
        \fn milxQtImage::isOriented()
        \brief Returns true if orientation is applied to the display of the image
    */
    inline bool isOriented()
    {
        return orientAct->isChecked();
    }

    /*!
        \fn milxQtImage::enableScale(QString title = "", const bool quiet = false, double minRange = 0.0, double maxRange = 0.0, int noOfLabels = 3)
        \brief Enable scale bar display with the title provided.

        Quiet Boolean is to prevent possible popups to ask user parameters.
    */
    virtual void enableScale(QString title = "", const bool quiet = false, double minRange = 0.0, double maxRange = 0.0, int noOfLabels = 3);
    /*!
        \fn milxQtImage::scaleDisplay(const bool forceDisplay = false)
        \brief Toggles the scale bar display.

        forceDisplay Boolean is to overide possible previous settings and display bar.
    */
    virtual void scaleDisplay(const bool forceDisplay = false);
    /*!
        \fn milxQtImage::enableCrosshair()
        \brief Enables the mouse pointer as a crosshair instead. Scene must be rendered before calling.
    */
    virtual inline void enableCrosshair()
    {
        if(rendered)
        {
            viewer->GetRenderWindow()->SetCurrentCursor(VTK_CURSOR_CROSSHAIR);
            crosshairAct->setChecked(true);
        }
    }
    /*!
        \fn milxQtImage::disableCrosshair()
        \brief Restores the mouse pointer to default.
    */
    virtual inline void disableCrosshair()
    {
        viewer->GetRenderWindow()->SetCurrentCursor(0);
        crosshairAct->setChecked(false);
    }

    /*!
        \fn milxQtImage::setActualNumberOfDimensions(const size_t dims)
        \brief Sets the actual number of dimensions of the image.

        All images loaded as 3D images or 3D vector images, this shows actual dimension thats maintained in file.
    */
    inline void setActualNumberOfDimensions(const size_t dims)
    {
        actualNumberOfDimensions = dims;
    }
    inline size_t getActualNumberOfDimensions()
    {
        return actualNumberOfDimensions;
    }

    /*!
        \fn milxQtImage::trackView(milxQtImage *windowToTrack, ViewType viewTo)
        \brief Enables tracking of the view (axial etc.) to imgToTrack provided.

        View tracking is linking the position of cursors in other milxQtImage objects to this one.
        This is useful for multi-view display where each window tracks different views (axial etc.) of the same data
    */
    void trackView(milxQtImage *windowToTrack, ViewType viewTo);

    //VTK Filters
    vtkSmartPointer<vtkImageData> butterWorthHighPass(vtkSmartPointer<vtkImageData> img);

public slots:
    void userEvent(QMouseEvent *event = NULL);
    void userEvent(QKeyEvent *event);
    void userEvent(QWheelEvent *event);
    /*!
        \fn milxQtImage::updateCoords(vtkObject *obj)
        \brief Picks the coordinates and pixel value from the current mouse position in the window.
    */
    virtual void updateCoords(vtkObject *obj);
    /*!
        \fn milxQtImage::updateSlice(vtkObject *obj)
        \brief Updates the display of the slice in the window.
    */
    void updateSlice(vtkObject *obj);
    /*!
        \fn milxQtImage::updateTrackedView(vtkObject *obj)
        \brief Updates the display of the slice in the window according to view tracking.
    */
    void updateTrackedView(vtkObject *obj);
    /*!
        \fn milxQtImage::contour()
        \brief Draw contour interactively on images
    */
    virtual void contour();

    /**
        \fn milxQtImage::autoLevel()
        \brief Auto window level the display. Uses Otsu threshold value.
    */
    void autoLevel();
/**
        \fn milxQtImage::setLevel(int level)
        \brief Set window level the display to proportion of maxValue.
    */
    void setLevel(int level);
    /*!
        \fn milxQtImage::GetWindowLevel()
        \brief Returns the internal image window levels data used for display.
    */
    inline vtkImageMapToWindowLevelColors* GetWindowLevel()
    {
        return viewer->GetWindowLevel();
    }
    /*!
        \fn milxQtImage::GetIntensityWindow()
        \brief Get the window for the image intensity tranfer function
    */
    double GetIntensityWindow();
    /*!
        \fn milxQtImage::SetIntensityWindow(double level)
        \brief Set the window for the image intensity tranfer function
    */
    void SetIntensityWindow(double window);
    /*!
        \fn milxQtImage::GetIntensityLevel()
        \brief Get the level for the image intensity tranfer function
    */
    double GetIntensityLevel();
    /*!
        \fn milxQtImage::SetIntensityLevel(double level)
        \brief Set the level for the image intensity tranfer function
    */
    void SetIntensityLevel(double level);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
    /**
        \fn milxQtImage::overlayContour(QString filename = "")
        \brief Overlays a labelled image on current image
    */
    void overlay(QString filename = "");
    /**
        \fn milxQtImage::overlayContour(QString filename = "")
        \brief Converts the labelled image into a contour and overlays on current image
    */
    void overlayContour(QString filename = "");
    /**
        \fn milxQtImage::computeContour()
        \brief Computes the contour of the image and displays it. Assumes image is binary currently.
    */
    void computeContour();
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
    /**
        \fn milxQtImage::blend(QString filename, float opacity = -1.0)
        \brief Display a blend of the current image with another image.
    */
    void blend(QString filename = "", float opacity = -1.0);
    /**
        \fn milxQtImage::blend(milxQtImage *imageToMatch, float opacity = -1.0)
        \brief Display a blend of the current image with another image. Image window can be directly passed into this member.
    */
    void blend(milxQtImage *imageToMatch, float opacity = -1.0);
    /**
        \fn milxQtImage::volumeRendering()
        \brief Display current image as volume rendering
    */
    void volumeRendering();
    /**
        \fn milxQtImage::imageInformation()
        \brief Displays the origin and other important information about the image (not the displayed image, but internal image used for computation)
    */
    void imageInformation();
    /**
        \fn milxQtImage::rescale()
        \brief Rescaled the intensities of an image to specified new max and min values.
    */
    void rescale();
    /**
        \fn milxQtImage::relabel()
        \brief Relabel the objects of a labelled image based on connectivity in consecutive order.
    */
    void relabel();
    /**
        \fn milxQtImage::histogramEqualisation()
        \brief Rescale the intensities of an image based on histogram equalisation.
    */
    void histogramEqualisation();
    /**
        \fn milxQtImage::gradientMagnitude()
        \brief Computes the gradient magnitude of the image and displays it.
    */
    void gradientMagnitude();
    /**
        \fn milxQtImage::sobelEdges()
        \brief Computes the Sobel edge detection of the image and displays it.
    */
    void sobelEdges();
    /**
        \fn milxQtImage::cannyEdges()
        \brief Computes the Canny edge detection of the image and displays it.
    */
    void cannyEdges();
    /**
        \fn milxQtImage::laplacian()
        \brief Computes the Laplacian of the image and displays it.
    */
    void laplacian();
    /**
        \fn milxQtImage::normalize()
        \brief Computes the normalization of the image and displays it.
    */
    void normalize();
    /**
        \fn milxQtImage::invertIntensity()
        \brief Computes the inverse intensities of the image and displays it.
    */
    void invertIntensity();
    /**
        \fn milxQtImage::matchInfo(milxQtImage *imageToMatch)
        \brief Matches the info of another image to the current image. See generateMatchedInformation() for more info.
    */
    void matchInfo(milxQtImage *imageToMatch);
    /**
        \fn milxQtImage::matchInfo(QString filename = "")
        \brief Matches the info of another image to the current image. See generateMatchedInformation() for more info.
    */
    void matchInfo(QString filename = "");
    /**
        \fn milxQtImage::matchHistogram(milxQtImage *imageToMatch)
        \brief Matches the histogram of another image to the current image. See generateMatchedHistogram() for more info.
    */
    void matchHistogram(milxQtImage *imageToMatch);
    /**
        \fn milxQtImage::matchHistogram(QString filename = "")
        \brief Matches the histogram of another image to the current image. See generateMatchedHistogram() for more info.
    */
    void matchHistogram(QString filename = "");
    /**
        \fn milxQtImage::resample(QString filename = "")
        \brief Resamples current image to another image.
    */
    void resample(QString filename = "");
    /**
        \fn milxQtImage::mask(QString filename = "")
        \brief Masks current image with another image.

        Supports Vector images also.
    */
    void mask(QString filename = "");
    /**
        \fn milxQtImage::subsample(size_t xSampleFactor = 0, size_t ySampleFactor = 0, size_t zSampleFactor = 0)
        \brief Downsamples current image by factors given.
    */
    void subsample(size_t xSampleFactor = 0, size_t ySampleFactor = 0, size_t zSampleFactor = 0);
    /**
        \fn milxQtImage::crop(QString filename = "")
        \brief Masks and crops the current image with another image.

        Supports Vector images also.
    */
    void crop(QString filename = "");
    /**
        \fn milxQtImage::resampleLabel(QString filename = "")
        \brief Resamples current image as a labelled image (toavoid interpolation artefacts) to another image.
    */
    void resampleLabel(QString filename = "");
    /**
        \fn milxQtImage::transform(QString filename, QString refImgFilename, bool inverse)
        \brief Transforms current image given ITK transform file.

        You need to provide not only the ITK transform file but the image to resample the current image to.
    */
    void transform(QString filename = "", QString refImgFilename = "", bool inverse = false);
    /**
        \fn milxQtImage::checkerBoard(QString filename = "", int numberOfSquares = 0)
        \brief Tiles another image with the current image using alternating squares for comparison. See generateCheckerBoard() for more info.
    */
    void checkerBoard(QString filename = "", int numberOfSquares = 0);
    /**
        \fn milxQtImage::checkerBoard(milxQtImage *img, int numberOfSquares = 0)
        \brief Tiles another image with the current image using alternating squares for comparison. See generateCheckerBoard() for more info.

        This version accepts image objects directly.
    */
    void checkerBoard(milxQtImage *img, int numberOfSquares = 0);
    /**
        \fn milxQtImage::distanceMap(bool signedDistance = true, bool inside = false)
        \brief Computes the distance map of the image and displays it.

        if signed distance then Signed Maurer distance maps is generated.
    */
    void distanceMap(bool signedDistance = true, bool inside = false);
    /**
        \fn milxQtImage::thresholdAbove(float value = 0, float level = 0)
        \brief Threshold the image for all values above a given value.
    */
    void thresholdAbove(float value = 0, float level = 0);
    /**
        \fn milxQtImage::thresholdBelow(float value = 0, float level = 0)
        \brief Threshold the image for all values below a given value.
    */
    void thresholdBelow(float value = 0, float level = 0);
    /**
        \fn milxQtImage::threshold(float value = 0, float blevel = 0, float alevel = 0)
        \brief Threshold the image for all values between a band to a given value.
    */
    void threshold(float value = 0, float blevel = 0, float alevel = 0);
    /**
        \fn milxQtImage::binaryThreshold(float value = 0, float blevel = 0, float alevel = 0)
        \brief Binary Threshold the image for all values between a band to a given value.
    */
    void binaryThreshold(float value = 0, float blevel = 0, float alevel = 0);
    /**
        \fn milxQtImage::otsu(int bins = 0)
        \brief Otsu Threshold the image for all values between a band to a given level.
    */
    void otsu(int bins = 0);
    /**
        \fn milxQtImage::otsuMultiple(int bins = 0, int labels = 0)
        \brief Otsu Multiple Threshold the image for all values between a band to a given level.
    */
    void otsuMultiple(int bins = 0, int labels = 0);
    /**
        \fn milxQtImage::flip(bool xAxis = false, bool yAxis = false, bool zAxis = false, bool aboutOrigin = true)
        \brief Flips the image and displays it. If all the arguments are false, then axes will be asked for using GUI.
    */
    void flip(bool xAxis = false, bool yAxis = false, bool zAxis = false, bool aboutOrigin = true);
    /**
        \fn milxQtImage::surface(const float value = numeric_limits<float>::max())
        \brief Converts the image to a surface using the Marching Cubes algorithm. Ideally suited for binary images.

        Note that this only emits a signal so that a class like milxQtMain can redirect it to the milxQtModel class for processing.
    */
    void surface(const float value = numeric_limits<float>::max());
    /**
        \fn milxQtImage::polyData()
        \brief Converts the image to polygonal data. Ideally suited for any type of image.

        Note that this only emits a signal so that a class like milxQtMain can redirect it to the milxQtModel class for processing.
    */
    void polyData();
    /**
        \fn milxQtImage::magnitude()
        \brief Shows the vector image as a scalar image given by the magnitude of the vectors.
    */
    void magnitude();
    /**
        \fn milxQtImage::component(int index = -1)
        \brief Shows the vector image as a scalar image given component index.
    */
    void component(int index = -1);
    /**
        \fn milxQtImage::pseudoImage()
        \brief Converts the vector image to a pseudo-image, where each vector component is 'coloured'.

        Note that this only emits a signal so that a class like milxQtMain can redirect it to the milxQtModel class for processing.
    */
    void pseudoImage();
    /**
        \fn milxQtImage::vectorField(int subsampleFactor = 0, float scaling = 0.0)
        \brief Converts the vector image to a vector field.

        Note that this only emits a signal so that a class like milxQtMain can redirect it to the milxQtModel class for processing.
    */
    void vectorField(int subsampleFactor = 0, float scaling = 0.0);
    /**
        \fn milxQtImage::streamLines()
        \brief Converts the vector image to a series of stream lines which all start at the current slice being viewed.

        Note that this only emits a signal so that a class like milxQtMain can redirect it to the milxQtModel class for processing.
    */
    void streamLines();
    /**
        \fn milxQtImage::anisotropicDiffusion()
        \brief Computes the anisotropic diffusion (smoothing) of the image and displays it.
    */
    void anisotropicDiffusion();
    /**
        \fn milxQtImage::gaussianSmooth()
        \brief Computes the Gaussian smoothing of the image and displays it.
    */
    void gaussianSmooth();
    /**
        \fn milxQtImage::median()
        \brief Computes the median image and displays it.
    */
    void median();
    /**
        \fn milxQtImage::zeros(const unsigned long xSize, const unsigned long ySize, const unsigned long zSize, milxQtImage *refImage = NULL)
        \brief Creates image of zeros of size given.
    */
    void zeros(const unsigned long xSize, const unsigned long ySize, const unsigned long zSize, milxQtImage *refImage = NULL);
    /**
        \fn milxQtImage::resize(double outputSpacing = 0.0)
        \brief Resizes image to spacing given. Currently only does isotropic. \todo Extend to support anisotropic sampling
    */
    void resize(double outputSpacing = 0.0);
    /**
        \fn milxQtImage::resize(const unsigned long xSize, const unsigned long ySize, const unsigned long zSize, milxQtImage *refImage = NULL)
        \brief Resizes image to size given by reference image.
    */
    void resize(const unsigned long xSize, const unsigned long ySize, const unsigned long zSize, milxQtImage *refImage = NULL);
    /**
        \fn milxQtImage::add(milxQtImage *img)
        \brief Adds img to current image.
    */
    void add(milxQtImage *img);
    /**
        \fn milxQtImage::add(QString filename = "")
        \brief Adds img to current image.
    */
    void add(QString filename = "");
    /**
        \fn milxQtImage::subtract(milxQtImage *img)
        \brief Subtracts img from current image.
    */
    void subtract(milxQtImage *img);
    /**
        \fn milxQtImage::subtract(QString filename = "")
        \brief Subtracts img from current image.
    */
    void subtract(QString filename = "");
    /**
        \fn milxQtImage::scale(float scaling)
        \brief scales the intensities of current images.
    */
    void scale(float scaling);
    /**
        \fn milxQtImage::convolve(milxQtImage *img)
        \brief Convolves img to current image.
    */
    void convolve(milxQtImage *img);

    //VTK Filters
    /**
        \fn milxQtImage::highpass()
        \brief Applies Butterworth high pass filter to image.

        The VTK image result is converted to an ITK image object at the end to allow further processing.
    */
    void highpass();

    //Display
    /*!
        \fn milxQtImage::interpolateDisplay(const bool quietly = false)
        \brief Toggles interpolation.

        Quietly makes sure the change is not reannounced to dependent windows. Useful for multi-view
    */
    void interpolateDisplay(const bool quietly = false);
    inline void disableInterpolateDisplay()
    {   viewer->GetImageActor()->InterpolateOn();   interpolateDisplay();   }
    inline void enableInterpolateDisplay()
    {   viewer->GetImageActor()->InterpolateOff();   interpolateDisplay();   }
    /*!
        \fn milxQtImage::applyOrientDisplay(const bool quietly = false)
        \brief Toggles applying the orientation matrix of the image.

        Not to be confused with showing the orientation marker with orientDisplay()
        Quietly makes sure the change is not reannounced to dependent windows. Useful for multi-view
    */
    void applyOrientDisplay(const bool quietly = false);
    inline void disableApplyOrientDisplay()
    {   orientAct->setChecked(false);   applyOrientDisplay();   }
    inline void enableApplyOrientDisplay()
    {   orientAct->setChecked(true);   applyOrientDisplay();   }
    /*!
        \fn milxQtImage::setDefaultOrientation(int orientMode)
        \brief Change orientation mode to one of the supported standards. Default: Radiological.

        0-Radiological: Feet first view
        1-Neurological: Head first view
    */
    void setDefaultOrientation(int orientMode);
    /*!
        \fn milxQtImage::setView(int viewMode)
        \brief Change view to view mode identified by number.
        0-axial, 1-coronal, 2-sagittal
    */
    void setView(int viewMode);
    /*!
        \fn milxQtImage::viewToXYPlane()
        \brief Change view to xy-plane.
    */
    virtual void viewToXYPlane();
    /*!
        \fn milxQtImage::viewToZXPlane()
        \brief Change view to zx-plane.
    */
    virtual void viewToZXPlane();
    /*!
        \fn milxQtImage::viewToZYPlane()
        \brief Change view to zy-plane.
    */
    virtual void viewToZYPlane();
    /*!
        \fn milxQtImage::updateLookupTable()
        \brief Sets the necessary LUTs to model view.
    */
    virtual void updateLookupTable();
    /**
        \fn milxQtImage::histogram(int bins = 256, float belowValue = 0, float aboveValue = 255, bool plotHistogram = true)
        \brief Computes the histogram of the image and prints info found.
    */
    void histogram(int bins = 256, float belowValue = 0, float aboveValue = 255, bool plotHistogram = true);
    /**
        \fn milxQtImage::surfacePlot()
        \brief Surface plot of slice
    */
    void surfacePlot();
    /*!
        \fn milxQtImage::setSlice(int slice)
        \brief Sets the slice of the volume to slice value given.

        Bounds is checked.
    */
    void setSlice(int slice);
    /*!
        \fn milxQtImage::getSlice()
        \brief Gets the current slice of the volume.
    */
    inline int getSlice()
    {
        return viewer->GetSlice();
    }
    /*!
        \fn milxQtImage::getIntensityWindow()
        \brief Get or set the window/level for the image intensity tranfert function
    */
    inline double getIntensityWindow()
    {
        return viewer->GetWindowLevel()->GetWindow();
    }
    /*!
        \fn milxQtImage::getIntensityLevel()
        \brief Get or set the window/level for the image intensity tranfert function
    */
    inline double getIntensityLevel()
    {
        return viewer->GetWindowLevel()->GetLevel();
    }
    /*!
        \fn milxQtImage::setIntensityWindow()
        \brief Get or set the window/level for the image intensity tranfert function
    */
    inline void setIntensityWindow(double window)
    {
        viewer->GetWindowLevel()->SetWindow(window);
    }
    /*!
        \fn milxQtImage::setIntensityLevel()
        \brief Get or set the window/level for the image intensity tranfert function
    */
    inline void setIntensityLevel(double level)
    {
        viewer->GetWindowLevel()->SetLevel(level);
    }
    /**
        \fn milxQtImage::refresh()
        \brief Refresh display
    */
    void refresh();
    /**
        \fn milxQtImage::reset()
        \brief Update data internally and refresh display
    */
    void reset();

    /*!
        \fn milxQtImage::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);

    /*!
        \fn milxQtImage::generateImage(const bool quietly = false)
        \brief Assigns the array data to the image and setups up the viewer.

        Quietly makes sure the change is not reannounced to dependent windows. Useful for multi-view
    */
    void generateImage(const bool quietly = false);
    /*!
        \fn milxQtImage::generateVoxelisedSurface(vtkSmartPointer<vtkPolyData> surfaceToVoxelise, double *bounds = NULL, double *spacing = NULL)
        \brief Converts or Voxelises a surface (vtkPolyData) to an image.
    */
    void generateVoxelisedSurface(vtkSmartPointer<vtkPolyData> surfaceToVoxelise, double *bounds = NULL, double *spacing = NULL);

    //Custom
    /*!
        \fn milxQtImage::customOperation()
        \brief Custom operation, data dependent for viewing in unified environment.

        For the image, this displays the image actor from the window selected.
    */
    virtual void customOperation();
    void updateModelActor(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command);
    void updateDisplay(QPointer<milxQtImage> img);

signals:
    /*!
        \fn milxQtImage::imageToSurface(vtkSmartPointer<vtkImageData>, const float )
        \brief Emit signal to compute image to surface process.
    */
    void imageToSurface(vtkSmartPointer<vtkImageData>, const float );
    /*!
        \fn milxQtImage::imageToPolyData(vtkSmartPointer<vtkImageData> )
        \brief Emit signal to compute image to poly data process.
    */
    void imageToPolyData(vtkSmartPointer<vtkImageData> );
    /*!
        \fn milxQtImage::imageToPseudoImage(vectorImageType::Pointer )
        \brief Emit signal to compute image to pseudo-image process.
    */
    void imageToPseudoImage(vectorImageType::Pointer );
    /*!
        \fn milxQtImage::imageToVectorField(vectorImageType::Pointer, floatImageType::Pointer, int, float )
        \brief Emit signal to compute image to vector field process.
    */
    void imageToVectorField(vectorImageType::Pointer, floatImageType::Pointer, int, float);
    /*!
        \fn milxQtImage::imageToTensorField(vectorImageType::Pointer, floatImageType::Pointer, int, float )
        \brief Emit signal to compute image to tensor field process.
    */
    void imageToTensorField(vectorImageType::Pointer, floatImageType::Pointer, int, float );
    /*!
        \fn milxQtImage::imageToStreamLines(vectorImageType::Pointer, floatImageType::Pointer )
        \brief Emit signal to compute image streamlines of vector field process. Every pixel in the float image is used to seed the lines.
    */
    void imageToStreamLines(vectorImageType::Pointer, floatImageType::Pointer );
    /*!
        \fn milxQtImage::imageToVolume(vtkSmartPointer<vtkImageData> , bool)
        \brief Emit signal to show the image as a volume.
    */
    void imageToVolume(vtkSmartPointer<vtkImageData> , bool);
    /*!
        \fn milxQtImage::imageToPlot(vtkSmartPointer<vtkImageData>, int)
        \brief Emit signal to show the image as a plot.
    */
    void imageToPlot(vtkSmartPointer<vtkImageData>, int);
    /*!
        \fn milxQtImage::modified(QPointer<milxQtImage> )
        \brief Emit signal to send updated image data. Used for updating dependent displays
    */
    void modified(QPointer<milxQtImage> );
    /*!
        \fn milxQtImage::coordinateChanged(int, int, int)
        \brief Emit signal with the coordinate user is changed pointing to.
    */
    void coordinateChanged(int, int, int);

protected:
    //Flags
    bool usingVTKImage; //!< using VTK image data?
    bool imported; //!< Imported before?
    bool appendedData; //!< Appended image data?
    bool eightbit; //!< Using eightbit data?
    bool rgb; //!< Using RGB data?
    bool vectorised; //!< Using Vector image data?
    bool viewerSetup; //!< has the viewer/window been setup (only done initial so is to not disturb users settings)
    bool volume; //!< is the image a volume?
    bool flipped; //!< Flip for display?
    bool track; //!< track the coordinates during user interaction

    //Image Related
    //ITK
    charImageType::Pointer imageChar; //!< Up to date 8-bit greyscale image data
    rgbImageType::Pointer imageRGB; //!< Up to date 32-bit image data (used only internally atm)
    floatImageType::Pointer imageFloat; //!< Up to date floating point image data
    vectorImageType::Pointer imageVector; //!< Up to date vector image data

    size_t actualNumberOfDimensions; //!< All images loaded as 3D images or 3D vector images, this shows actual dimension
    ViewType viewToTrack; //!< In tracking mode, what slice to show

    //VTK
    vtkSmartPointer<vtkImageViewer3> viewer; //!< VTK Viewer handler, Smart Pointer
    vtkSmartPointer<vtkImageData> imageData; //!< Points to the current VTK Image Data, Smart Pointer
    vtkSmartPointer<vtkImageAppend> imageDataAppended; //!< Appended Data
    vtkSmartPointer<vtkImagePermute> permute; //!< Permute axis class, Smart Pointer
    QList<ModelActorItem> modelActors; //!< Model actors being displayed in image view

    itkEventQtObserver::Pointer observeProgress; //!< Observer for the Qt event loop

    ///Other Variables
    double meanValue; //!< Average data value currently held
    double minValue; //!< min value in image
    double maxValue; //!< max value in image

    //Context Menu
    //menu defined in milxQtWindow
    //------------------
    QMenu* operateMenu; //!< Operate Menu
    QAction* rescaleAct; //!< Action for contouring image
    QAction* equaliseAct; //!< Action for contouring image
    QAction* computeContourAct; //!< Action for contouring image
    QAction* smoothAct; //!< Action for smoothing of image
    QAction* gaussianAct; //!< Action for Gaussian smoothing of image
    QAction* medianAct; //!< Action for median smoothing of image
    QAction* gradMagAct; //!< Action for gradient magnitude of image
    QAction* sobelAct; //!< Action for sobel edges of image
    QAction* cannyAct; //!< Action for canny edges of image
    QAction* laplacianAct; //!< Action for Laplacian of image
    QAction* highPassAct; //!< Action for high pass filtering of image
    QAction* normAct; //!< Action for normalization of image
    QAction* invertAct; //!< Action for invert intensity of image
    QAction* relabelAct; //!< Action for relabelling image
    //------------------
    QMenu* transformMenu; //!< Transform Menu
    QAction* matchAct; //!< Action for matching info of image to another image
    QAction* matchHistAct; //!< Action for matching histogram of image
    QAction* resampleSpacingAct; //!< Action for resampling image based on spacing
    QAction* resampleAct; //!< Action for resampling image
    QAction* resampleLabelAct; //!< Action for resampling image
    QAction* subsampleAct; //!< Action for downsampling image
    QAction* transformAct; //!< Action for transforming image
    QAction* maskAct; //!< Action for resampling image
    QAction* cropAct; //!< Action for auto cropping image
    QAction* checkerAct; //!< Action for checkerboard of image
    QAction* distMapAct; //!< Action for distance map of image
    QAction* flipAct; //!< Action for flipping image
    QAction* surfaceAct; //!< Action for image to surface
    QAction* polyDataAct; //!< Action for image to poly data
    //------------------
    QMenu* thresholdMenu; //!< Threshold Menu
    QAction* otsuAct; //!< Otsu Threshold
    QAction* otsuMultipleAct; //!< Otsu Threshold
    QAction* binaryAct; //!< Binary Threshold
    QAction* bandAct; //!< Threshold inbetween band
    QAction* aboveAct; //!< Threshold from above
    QAction* belowAct; //!< Threshold from below
    //------------------
    QMenu* vectorMenu; //!< Vector Menu
    QAction* vectorMagnitudeAct; //!< Action for displaying magnitude of vector images
    QAction* vectorComponentAct; //!< Action for displaying components of vector images
    QAction* pseudoImageAct; //!< Action for displaying vector/tensor fields as a pseudo image
    QAction* vectorFieldAct; //!< Action for displaying vector/tensor fields
    QAction* streamLinesAct; //!< Action for displaying stream lines
    //------------------
    QAction* levelAct; //!< Action for auto-levelling gamma for display
    QAction* overlayAct; //!< Action for overlaying labelled image
    QAction* overlayContourAct; //!< Action for overlaying labelled image as contour
    QAction* blendAct; //!< Action for blending images
    QAction* volRenderAct; //!< Action for volume rendering
    QAction* histogramAct; //!< Action for displaying histogram
    QAction* surfacePlotAct; //!< Action for displaying surface plot
    QAction* infoAct; //!< Action for displaying information about the image.
    QAction* interpolateAct; //!< Interpolate image?
    QAction* orientAct; //!< Orient image?

    /*!
    	\fn milxQtImage::updateData(const bool orient = true)
    	\brief Ensures the internal visualisation data is up to date.
    */
    void updateData(const bool orient = true);
    /*!
    	\fn milxQtImage::setupEvents()
        \brief Executes common events setup and connections code for image viewing.

        This includes keys pressed for window.
    */
    void setupEvents();
    //Internal Members
    /*!
        \fn milxQtImage::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
        \fn milxQtImage::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
    	\fn milxQtImage::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);
    /*!
    	\fn milxQtImage::basicContextMenu()
    	\brief Return the basic context menu with the milxQtImage class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* basicContextMenu();
    /*!
    	\fn milxQtImage::operationsMenu()
    	\brief Return the operations menu with the milxQtImage class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* operationsMenu();
    /*!
    	\fn milxQtImage::thresholdsMenu()
    	\brief Return the thresholds menu with the milxQtImage class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* thresholdsMenu();
    /*!
    	\fn milxQtImage::transformsMenu()
    	\brief Return the transforms menu with the milxQtImage class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* transformsMenu();
    /*!
    	\fn milxQtImage::vectorsMenu()
    	\brief Return the vector imaging menu with the milxQtImage class ordering. This is for the benefit of derived class context menus and/or extensions.
    */
    QMenu* vectorsMenu();
    /*!
        \fn milxQtImage::dropEvent(QDropEvent *event)
        \brief Part of the Drag and Drop feature members. Opens the dropped files.
    */
    void dropEvent(QDropEvent *event);

    /*!
        \fn milxQtImage::SetupWidgets(vtkRenderWindowInteractor *interactor)
        \brief Update the interactor so that the widgets are usable. Resizing widges for image objects.
    */
    virtual void SetupWidgets(vtkRenderWindowInteractor *interactor);
    /*!
        \fn milxQtImage::updateViewer(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
        \brief Update the viewer image data to reflect changes.
    */
    inline void updateViewer(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
    {
    #if VTK_MAJOR_VERSION <= 5
        viewer->SetInput(imageData);
    #else
        viewer->SetInputData(imageData);
    #endif
    }
    /*!
        \fn milxQtImage::getOpenFilename(const QString labelForDialog = "Select Image", QString exts = "")
        \brief Use the QFileDialog to get an image filename.
    */
    QString getOpenFilename(const QString labelForDialog = "Select Image", QString exts = "");

private:

};

#endif // MILXQTIMAGE_H
