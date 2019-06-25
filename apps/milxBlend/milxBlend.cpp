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
/**
    \brief This program blends image and labels into full 3D images for comparison.
    \author Shekhar S. Chandra, 2016

    Usage:
    
*/
//Qt
#include <QApplication>

#include "milxQtImage.h"
//#include "milxQtMain.h"
#include "milxQtFile.h"

#include <vtkNew.h>
#include <vtkMath.h>
#include <vtkImageBlend.h>
#include <vtkImageCast.h>

#include <tclap/CmdLine.h> //Command line parser library

typedef std::vector< std::string >::iterator stringiterator;

using namespace TCLAP;

//prototypes
QPointer<milxQtImage> imagesBlend(milxQtImage *firstImg, milxQtImage *secondImg, float opacity);
template<typename TImage>
void CorrectCESTMapping(itk::SmartPointer<TImage> img, float scaleFactor = 10);

int main(int argc, char* argv[])
{
    QApplication app(argc,argv);
    QMainWindow mainWindow;//app takes ownership so need to be on stack, not heap

    //---------------------------
    ///Program Info
    milx::PrintInfo("--------------------------------------------------------");
    milx::PrintInfo("SMILI Blend Tool for Medical Imaging");
    milx::PrintInfo("(c) Copyright Chandra et al., 2015.");
    milx::PrintInfo("University of Queensland, Australia.");
    milx::PrintInfo("Australian e-Health Research Centre, CSIRO, Australia.");
    milx::PrintInfo("SMILI Version: " + milx::NumberToString(milx::Version));
    milx::PrintInfo("milxQt Version: " + milx::NumberToString(milxQtVersion));
    milx::PrintInfo("Application Version: 1.0");
    milx::PrintInfo("--------------------------------------------------------\n");

    //---------------------------
    ///Process Arguments
    CmdLine cmd("An blend tool for images", ' ', milx::NumberToString(milx::Version));

    ///Optional
    ValueArg<std::string> outputArg("o", "output", "Output Screenshot name", false, "screenshot.png", "Output");
    ValueArg<std::string> maskArg("m", "mask", "Mask the blend image", false, "mask.nii.gz", "Mask");
    ValueArg<std::string> labelArg("l", "label", "Overlay binary image/label contour on the image", false, "labelling.nii.gz", "Label");
    //ValueArg<std::string> transformArg("t", "transform", "Transform (ITK Format) to apply to objects being rendered", false, "rigid.txt", "Transform");
    ValueArg<std::string> scalarBarArg("", "scalarbar", "Enable the scalar bar for display with given name.", false, "Scalars", "Scalar Bar");
    ValueArg<float> opacityArg("", "opacity", "Set the opacity of the images", false, 0.8, "Opacity");
    ValueArg<float> minArg("", "min", "Set the minimum value for scalars on image.", false, 0.0, "Minimum Scalar");
    ValueArg<float> maxArg("", "max", "Set the maximum value for scalars on image.", false, 1.0, "Maximum Scalar");
    ValueArg<float> redArg("", "redbackground", "Set the redness value (0-1) for the background.", false, 1.0, "Red Component");
    ValueArg<float> greenArg("", "greenbackground", "Set the greenness value (0-1) for the background.", false, 1.0, "Green Component");
    ValueArg<float> blueArg("", "bluebackground", "Set the blueness value (0-1) for the background.", false, 1.0, "Blue Component");
    ValueArg<int> heightArg("y", "height", "Set the height of the onscreen window.", false, 512, "Height");
    ValueArg<int> widthArg("x", "width", "Set the width of the onscreen window.", false, 512, "Width");
    ///Switches
    SwitchArg jetArg("", "jet", "Change colourmap of the scalars to the Jet map", false);
    SwitchArg vtkArg("", "vtk", "Change colourmap of the scalars to blue-red (rainbow VTK) map", false);
    SwitchArg hsvArg("", "hsv", "Change colourmap of the scalars to blue-red (rainbow HSV) map", false);
    SwitchArg rainbowArg("", "rainbow", "Change colourmap of the scalars to the rainbow map", false);
    SwitchArg spectralArg("", "spectral", "Change colourmap of the scalars to the spectral map", false);
    SwitchArg nihArg("", "NIH", "Change colourmap of the scalars to NIH", false);
    SwitchArg fireArg("", "NIH_FIRE", "Change colourmap of the scalars to NIH Fire", false);
    SwitchArg aalArg("", "AAL", "Change colourmap of the scalars to AAL", false);
    SwitchArg fsArg("", "FS", "Change colourmap of the scalars to FreeSurfer", false);
    SwitchArg hotArg("", "HOT", "Change colourmap of the scalars to HOT", false);
    SwitchArg coolArg("", "COOL", "Change colourmap of the scalars to COOL", false);
    SwitchArg onscreenArg("", "onscreen", "Enable on screen rendering, i.e. display the rendering as an interactive window.", false);
    SwitchArg whiteArg("", "white", "Make background white rather than default gradient colour.", false);
    //SwitchArg inverseArg("", "inverse", "Use the inverse transform when transform file is provided.", false);
    SwitchArg clampScalarsArg("", "clampscalars", "Clamp the scalars of the images to background colour (outside colourmap) from min/max values provided.", false);
    SwitchArg clampScalarsBelowArg("", "clampbelow", "Clamp the scalars of the images to background colour (outside colourmap) below min values provided.", false);
    SwitchArg autoColourArg("", "autocolour", "Colour each additional image automatically based on order.", false);
    SwitchArg axialArg("", "axial", "Show axial view onscreen based on position of first surface.", false);
    SwitchArg coronalArg("", "coronal", "Show coronal onscreen view based on position of first surface.", false);
    SwitchArg saggitalArg("", "saggital", "Show saggital onscreen view based on position of first surface.", false);

    //Scale modes
    SwitchArg gagCestArg("", "gagcest", "Apply scaling to image that blends on top. GagCEST scaling.", false);
    
    ///Mandatory
    UnlabeledMultiArg<std::string> multinames("images", "Images to be blended over", true, "Images");
    ValueArg<std::string> imageArg("i", "image", "Image to blend on top of input images", true, "image.nii.gz", "Image");

    ///Add argumnets
    cmd.add( multinames );
    cmd.add( outputArg );
    cmd.add( imageArg );
    cmd.add( maskArg );
#ifdef ITK_REVIEW //Review only members
    cmd.add( labelArg );
#endif
    //cmd.add( transformArg );
    cmd.add( minArg );
    cmd.add( maxArg );
    cmd.add( redArg );
    cmd.add( greenArg );
    cmd.add( blueArg );
    cmd.add( opacityArg );
    cmd.add( heightArg );
    cmd.add( widthArg );

    cmd.add( jetArg );
    cmd.add( vtkArg );
    cmd.add( hsvArg );
    cmd.add( rainbowArg );
    cmd.add( spectralArg );
    cmd.add( nihArg );
    cmd.add( fireArg );
    cmd.add( aalArg );
    cmd.add( fsArg );
    cmd.add( hotArg );
    cmd.add( coolArg );
    cmd.add( scalarBarArg );
    cmd.add( onscreenArg );
    cmd.add( whiteArg );
    //cmd.add( inverseArg );
    //cmd.add( clampScalarsArg );
    //cmd.add( clampScalarsBelowArg );
    //cmd.add( autoColourArg );
    cmd.add( axialArg );
    cmd.add( coronalArg );
    cmd.add( saggitalArg );

    cmd.add( gagCestArg );

    ///Parse the argv array.
    cmd.parse( argc, argv );

    ///Get the value parsed by each arg.
    //Filenames of surfaces
    std::vector<std::string> filenames = multinames.getValue();
    const std::string outName = outputArg.getValue();
    const std::string imageName = imageArg.getValue();
    const std::string maskName = maskArg.getValue();
    const std::string scalarMaskName = maskArg.getValue();
    const std::string labelName = labelArg.getValue();
    //const std::string transformName = transformArg.getValue();
    const std::string barName = scalarBarArg.getValue();
    const float opacity = opacityArg.getValue();
    const float minValue = minArg.getValue();
    const float maxValue = maxArg.getValue();
    const float redValue = redArg.getValue();
    const float greenValue = greenArg.getValue();
    const float blueValue = blueArg.getValue();
    const int windowHeight = heightArg.getValue();
    const int windowWidth = widthArg.getValue();

    //Check arguments
    //Most of the checking is done by TCLAP
    /*if(inverseArg.isSet())
    {
        if(!transformArg.isSet())
        {
            cerr << "Error in arguments! Inverse argument needs to be used with the transform argument." << endl;
            exit(EXIT_FAILURE);
        }
    }*/

    ///Read
    //Reader object
    QScopedPointer<milxQtFile> reader(new milxQtFile);
    bool errorReading = false;
    bool success = false;

    std::cout << "Total Images: " << filenames.size() << std::endl;
    if (filenames.empty())
    {
      milx::PrintError("No image file names provided. Exiting.");
      exit(EXIT_FAILURE);
    }
    for (stringiterator name = filenames.begin(); name != filenames.end(); name++)
      std::cout << *name << ", ";
    std::cout << std::endl;

    //Read info
    std::cerr << "Reading Header for type info etc." << std::endl;
    std::string pixelType, componentType;
    size_t dimensions = 3;
    if (!milx::File::ReadImageInformation(filenames[0], pixelType, componentType, dimensions))
    {
      milx::PrintError("Failed Reading First Image. Check the image type/file. Exiting.");
      exit(EXIT_FAILURE);
    }

    std::cerr << "Reading Images... ";
    std::vector< itk::SmartPointer<charImageType> > labelledCollection;
    std::vector< itk::SmartPointer<floatImageType> > collection;
    bool labelledImages = false;
    if (componentType == "unsigned_char" || componentType == "unsigned char")
    {
      milx::PrintInfo("Detected labelled images.");
      if (!milx::File::OpenImages<charImageType>(filenames, labelledCollection)) //Error NOT printed inside
      {
        milx::PrintError("Failed Reading Labelled Images. Exiting.");
        exit(EXIT_FAILURE);
      }
      labelledImages = true;
      std::cout << "Read " << labelledCollection.size() << " labels" << std::endl;
    }
    else
    {
      milx::PrintInfo("Detected floating point images.");
      if (!milx::File::OpenImages<floatImageType>(filenames, collection)) //Error NOT printed inside
      {
        milx::PrintError("Failed Reading Images. Exiting.");
        exit(EXIT_FAILURE);
      }
      std::cout << "Read " << collection.size() << " images" << std::endl;
    }

    //Read transform
    //vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New(); //transform matrix
    //    transform->Identity();
    //    transform->PostMultiply();
    //if(transformArg.isSet())
    //{
    //    vtkSmartPointer<vtkMatrix4x4> matrix = milx::File::OpenITKTransform(transformName);
    //    transform->Concatenate(matrix);
    //    if(inverseArg.isSet())
    //        transform->Inverse();
    //}
/*
    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    if(autoColourArg.isSet())
    {
        lookupTable->SetTableRange(0.0, n+1);
        lookupTable->Build();
    }*/

    if(labelledCollection.empty() && collection.empty())
    {
        milx::PrintError("No Images Read. Exiting.");
        exit(EXIT_FAILURE);
    }

    //Read Blend Image
    QScopedPointer<milxQtImage> img(new milxQtImage);  //hierarchy deletion
    //vtkSmartPointer<vtkMatrix4x4> orientTransform = vtkSmartPointer<vtkMatrix4x4>::New();

    //Read image
    //vtkSmartPointer<vtkTransform> transform2 = vtkSmartPointer<vtkTransform>::New(); //transform matrix
    //    transform2->Identity();
    //    transform2->PostMultiply();
    if(imageArg.isSet())
    {
        cout << ">> Blend: Reading Image" << endl;
        errorReading = false;
        success = reader->openImage(imageName.c_str(), img.data());

        ///Generate Views
        if(success)
        {
            //Generate the images
            img->setName(imageName.c_str());
            img->generateImage();

#ifdef ITK_REVIEW //Review only members
            if(labelArg.isSet())
            {
                img->overlayContour(labelName.c_str()); //Generates an RGB image, from VTK class
            }
#endif
//            img->flip(false, true, false);
        }
        else
            errorReading = true;

        if(errorReading)
        {
            cerr << "Error Reading the image file. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        ///Setup orientation transform matrices
        //Flip y to make same coordinates as image
        //orientTransform->DeepCopy(img->getTransformMatrix());
        //orientTransform->Invert();
        //transform2->Concatenate(orientTransform);
//        transform2->Concatenate(transform->GetMatrix());
        //cout << ">> Overlay: Transforming Actors" << endl;
    }

    //Checks
    if (!img->GetFloatImage())
    {
      QMessageBox msgBox;
      msgBox.setText("Mapping image not float type.");
      msgBox.setInformativeText("Image " + img->strippedBaseName() + " values were unexpected. Please Check image pixel values.");
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setDefaultButton(QMessageBox::Ok);
      msgBox.exec();
    }

    //mask the map
    if (maskArg.isSet())
        img->mask(maskName.c_str());

    //Read Mask
    /*QScopedPointer<milxQtImage> mask(new milxQtImage);  //hierarchy deletion
    if (maskArg.isSet())
    {
        cout << ">> Blend: Reading Mask" << endl;
        errorReading = false;
        success = reader->openImage(maskName.c_str(), mask.data());

        ///Generate Views
        if (success)
        {
          //Generate the images
          mask->setName(maskName.c_str());
          mask->generateImage();
          //            img->flip(false, true, false);
        }
        else
          errorReading = true;

        if (errorReading)
        {
          cerr << "Error Reading the mask file. Exiting." << endl;
          exit(EXIT_FAILURE);
        }
    }*/

    ///Display
    ///Setup image display
    QScopedPointer<milxQtImage> image(new milxQtImage);  //hierarchy deletion
    if (labelledImages)
      image->setData(labelledCollection[0]);
    else
      image->setData(collection[0]);
    image->setName(filenames[0].c_str());
    image->generateImage();

    double range[2];
    image->GetOutput()->GetScalarRange(range);
    if (minArg.isSet())
      range[0] = minValue;
    if (maxArg.isSet())
      range[1] = maxValue;

    ///Blend
    image->colourMapToGray();
    image->autoLevel(0.95);
    qApp->processEvents();

    //Correct the blending image
    if(gagCestArg.isSet())
    {
        cout << ">> Blend: GagCEST scaling applied" << endl;
        //CEST[%]=(RGB-2048)/2000*SF
        CorrectCESTMapping<floatImageType>(img->GetFloatImage());
        img->generateImage();
    }

    ///Colour maps
    cout << ">>> Blend: Setting Colourmap" << endl;
    img->colourMapToJet();
    if (vtkArg.isSet())
      img->colourMapToVTK();
    if (hsvArg.isSet())
      img->colourMapToHSV();
    if (rainbowArg.isSet())
      img->colourMapToRainbow();
    if (spectralArg.isSet())
      img->colourMapToSpectral();
    if (nihArg.isSet())
      img->colourMapToNIH();
    if (fireArg.isSet())
      img->colourMapToNIH_Fire();
    if (aalArg.isSet())
      img->colourMapToAAL();
    if (fsArg.isSet())
      img->colourMapToFS();
    if (hotArg.isSet())
      img->colourMapToHOT();
    if (coolArg.isSet())
      img->colourMapToCOOL();

    if (scalarBarArg.isSet())
    {
      if (minArg.isSet() || maxArg.isSet())
        img->enableScale(barName.c_str(), true, range[0], range[1]);
      else
        img->enableScale(barName.c_str(), true);
    }

    //Visualise Map
    //img->GetWindowLevel()->PassAlphaToOutputOn();
    img->autoLevel(0.95);
    img->enableScale("Values", true);
    qApp->processEvents();

    //Blend
    cout << ">> Blend: Begining ..." << endl;
    QPointer<milxQtImage> blendImage = imagesBlend(image.data(), img.data(), opacity);

    //Transfer scalar bar
    blendImage->AddActor(img->GetScalarBar());
    mainWindow.setCentralWidget(blendImage);

    //Hide Mapping
    img->hide();

    if(whiteArg.isSet())
        blendImage->background(true);

    //setup background
    double backgroundColours[3];
    blendImage->GetBackground(backgroundColours[0], backgroundColours[1], backgroundColours[2]);
    if(redArg.isSet())
        backgroundColours[0] = redValue;
    if(greenArg.isSet())
        backgroundColours[1] = greenValue;
    if(blueArg.isSet())
        backgroundColours[2] = blueValue;
    if(redArg.isSet() || blueArg.isSet() || greenArg.isSet())
        blendImage->SetBackground(backgroundColours[0], backgroundColours[1], backgroundColours[2]);

    //Apply transform to main image last to ensure sub images are correctly transform also
    //if(transformArg.isSet())
    //    image->SetTransform(transform);

    cout << ">> Blend: Rendering" << endl;
    if (onscreenArg.isSet())
    {
        cout << ">> Blend: Displaying Onscreen" << endl;
        if (axialArg.isSet())
            blendImage->viewToAxial();
        if (coronalArg.isSet())
            blendImage->viewToCoronal();
        if (saggitalArg.isSet())
            blendImage->viewToSagittal();
        blendImage->SetSize(windowHeight, windowWidth);
        mainWindow.resize(windowWidth, windowHeight);
        blendImage->show();
        mainWindow.show();
    }
    else
        cout << ">> Blend: Output to File" << endl;

    QScopedPointer<milxQtFile> writer(new milxQtFile); //Smart deletion
    if (!onscreenArg.isSet())
    {
        /*vtkNew<vtkImageCast> caster;
        caster->SetInputData(blendImage->GetOutput());
        caster->SetOutputScalarTypeToUnsignedChar();
        caster->Update();*/
        milx::PrintDebug("Output to be written is of type " + std::string(blendImage->GetOutput()->GetScalarTypeAsString()));
        milx::PrintDebug("Output has number of components of " + milx::NumberToString(blendImage->GetOutput()->GetNumberOfScalarComponents()));
        writer->saveImage(outName.c_str(), blendImage->GetOutput());
    }
    cout << ">> Complete" << endl;

    //blendImage->OffScreenRenderingOff(); //Required to prevent double-free*/
    if(!onscreenArg.isSet())
        return EXIT_SUCCESS;
    else
        return app.exec();
}

QPointer<milxQtImage> imagesBlend(milxQtImage *firstImg, milxQtImage *secondImg, float opacity)
{
  vtkSmartPointer<vtkImageMapToColors> filterColorsImage = vtkSmartPointer<vtkImageMapToColors>::New();
  filterColorsImage->SetLookupTable(firstImg->GetLookupTable());
#if VTK_MAJOR_VERSION <= 5
  filterColorsImage->SetInput(firstImg->GetOutput());
#else
  filterColorsImage->SetInputData(firstImg->GetOutput());
#endif
  filterColorsImage->PassAlphaToOutputOn();
  filterColorsImage->Update();
  vtkSmartPointer<vtkImageData> ucharData1 = filterColorsImage->GetOutput();

  int initialExtent[6];
  firstImg->GetOutput()->GetExtent(initialExtent);

  // Combine the images (blend takes multiple connections on the 0th input port)
  vtkSmartPointer<vtkImageBlend> blend = vtkSmartPointer<vtkImageBlend>::New();
  blend->SetOpacity(0, 1.0);
#if VTK_MAJOR_VERSION <= 5
  blend->AddInput(ucharData1);
#else
  blend->AddInputData(ucharData1);
#endif
  //        blend->SetBlendModeToCompound();
  blend->SetBlendModeToNormal();

  milx::PrintInfo("Blending with Image: " + secondImg->strippedName().toStdString() + " with Opacity: " + milx::NumberToString(opacity));
  //        vtkSmartPointer<vtkImageData> ucharData2 = secondImg->GetWindowLevel()->GetOutput();
  vtkSmartPointer<vtkImageMapToColors> filterColorsOverlay = vtkSmartPointer<vtkImageMapToColors>::New();
  filterColorsOverlay->SetLookupTable(secondImg->GetLookupTable());
#if VTK_MAJOR_VERSION <= 5
  filterColorsOverlay->SetInput(secondImg->GetOutput());
#else
  filterColorsOverlay->SetInputData(secondImg->GetOutput());
#endif
  filterColorsOverlay->PassAlphaToOutputOn();
  filterColorsOverlay->Update();
  vtkSmartPointer<vtkImageData> ucharData2 = filterColorsOverlay->GetOutput();

  if (!secondImg->GetLookupTable())
    milx::PrintWarning("Colourmap is not set. Please set a colour map to ensure proper blending.");

  int actualExtent[6];
  secondImg->GetOutput()->GetExtent(actualExtent);

  if (initialExtent[1] == actualExtent[1] && initialExtent[3] == actualExtent[3] && initialExtent[5] == actualExtent[5])
  {
#if VTK_MAJOR_VERSION <= 5
    blend->AddInput(ucharData2);
#else
    blend->AddInputData(ucharData2);
#endif
    blend->SetOpacity(1, opacity);
  }
  else
    milx::PrintError("Images are not the same size. Skipping.");

  milx::PrintInfo("Blending");
  blend->Update();
  milx::PrintDebug("Number of components: " + milx::NumberToString(blend->GetOutput()->GetNumberOfScalarComponents()));

  //Convert to RGB from RGBA for ITK medical image formats
  /*double tblRange[2] = { 0, 255 };
  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetRange(tblRange[0], tblRange[1]);
  lut->SetNumberOfTableValues(256);
  for (int i = 0; i <= VTK_UNSIGNED_CHAR_MAX; i++)
  {
    double doubleIndex = static_cast<double>(i);
    lut->SetTableValue(i, VTK_UNSIGNED_CHAR_MAX - doubleIndex, VTK_UNSIGNED_CHAR_MAX - doubleIndex, VTK_UNSIGNED_CHAR_MAX - doubleIndex, 1.0);
  }
  lut->Build();

  vtkSmartPointer<vtkImageMapToColors> rgbImage = vtkSmartPointer<vtkImageMapToColors>::New();
  rgbImage->SetLookupTable(lut);
#if VTK_MAJOR_VERSION <= 5
  rgbImage->SetInput(blend->GetOutput());
#else
  rgbImage->SetInputData(blend->GetOutput());
#endif
  rgbImage->PassAlphaToOutputOff();
  rgbImage->SetOutputFormatToRGB();
  rgbImage->Update();
  milx::PrintDebug("Blend result is of type " + std::string(rgbImage->GetOutput()->GetScalarTypeAsString()));
  milx::PrintDebug("Blend result has number of components of " + milx::NumberToString(rgbImage->GetOutput()->GetNumberOfScalarComponents()));
  */

  QPointer<milxQtImage> blendResult = new milxQtImage;
  blendResult->setName("Blended Images");
  blendResult->setData(blend->GetOutput());
  blendResult->generateImage();

  return blendResult;
}

template<typename TImage>
void CorrectCESTMapping(itk::SmartPointer<TImage> img, float scaleFactor)
{
  itk::ImageRegionIterator<TImage> imageIterator(img, img->GetLargestPossibleRegion());

  //CEST[%]=(RGB-2048)/2000*SF
  while (!imageIterator.IsAtEnd())
  {
    // Get the value of the current pixel
    float val = (imageIterator.Get() - 2048) / 2000 * scaleFactor;
    //std::cout << (int)val << std::endl;

    // Set the current pixel to white
    imageIterator.Set(val);

    ++imageIterator;
  }
}
