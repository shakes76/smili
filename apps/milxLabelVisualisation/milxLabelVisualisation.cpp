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
    \brief This program is specialised to display labelled images.
    \author Shekhar S. Chandra, 2014

    Linearise can be used to evenly space the label values across the selected colourmap. Use division and offset to 'place' the labels in this linear scale.
    This is useful when label values are very close together.

    Simple Example: milxLabelVisualisation --white -l seg.nii.gz -o cartilage.png

    Example: ./bin/milxLabelVisualisation --white -l bin/Anon14_20111028_001_012_t2_trufi3d_we_sag_iso_output.nii.gz  -o cartilage.png --nohuman --postprocess -s ./bin/Atlas_MRI_surface_R_fem_noscalars.vtk -t ./bin/Anon14_20111028_001_012_t2_trufi3d_we_sag_iso_output_Transform_initial_rigid_Atlas_000_to_Image.txt --inverse --above 160 --below 125 --loadviewfile ./bin/fem_atlas.cam
*/
//Qt
#include <QApplication>
#include <QMainWindow>

#include "milxQtImage.h"
#include "milxQtModel.h"
#include "milxQtPlot.h"
#include "milxQtFile.h"

#include <vtkWindowToImageFilter.h>

#include <tclap/CmdLine.h> //Command line parser library

using namespace TCLAP;

int main(int argc, char* argv[])
{
    QApplication app(argc,argv);
    QMainWindow mainWindow;//app takes ownership so need to be on stack, not heap

    ///Process Arguments
    CmdLine cmd("A visualisation tool for labelled images", ' ', milx::NumberToString(milx::Version));

    ///Mandatory
    ValueArg<std::string> labelledImageArg("l", "labels", "Labelled image to visualise", true, "seg.nii.gz", "Labels");
    ///Optional
    ValueArg<std::string> outputArg("o", "output", "Output Screenshot name", false, "screenshot.png", "Output");
    ValueArg<std::string> surfaceArg("s", "surface", "Surface to overlay with labelled image", false, "surface.vtk", "Surface");
    ValueArg<std::string> transformArg("t", "transform", "Transform (ITK Format) to apply to objects being rendered", false, "rigid.txt", "Transform");
    ValueArg<std::string> loadViewFileArg("", "loadviewfile", "Load saved view from file (use onscreen mode to save view files)", false, "camera.cam", "Load View File");
    ValueArg<float> aboveArg("", "above", "Add above value to thresholding of the labels. Controls max of scalar range of colourmap also.", false, 255, "Above");
    ValueArg<float> belowArg("", "below", "Add below value to thresholding of the labels. Controls min of scalar range of colourmap also.", false, 0, "Below");
    ValueArg<int> heightArg("y", "height", "Set the height of the window.", false, 600, "Height");
    ValueArg<int> widthArg("x", "width", "Set the width of the window.", false, 800, "Width");
    ValueArg<int> linearDivisionArg("", "lineardivision", "Evenly space the label values across N to get full use of colourmap.", false, 16, "Linear Division");
    ValueArg<int> linearOffsetArg("", "linearoffset", "Offset the even spacing of the label values across N to get full use of colourmap.", false, 1, "Linear Offset");
    ///Switches
    SwitchArg onscreenArg("", "onscreen", "Enable on screen rendering, i.e. display the rendering as an interactive window.", false);
    SwitchArg inverseArg("", "inverse", "Use the inverse transform when transform file is provided.", false);
    SwitchArg whiteArg("", "white", "Make background white rather than default gradient colour.", false);
    SwitchArg loadViewArg("", "loadview", "Load saved view (use smilx or onscreen render mode to view and save with Right Click->View->Save View", false);
    SwitchArg humanArg("", "nohuman", "Disable human orientation glyph.", false);
    SwitchArg processArg("", "postprocess", "Process the resulting surfaces from labels to reduce vertices and smooth etc.", false);
    SwitchArg volumeArg("", "volume", "Display using volume rendering instead of marching cude iso-surfaces.", false);
    SwitchArg lineariseArg("", "linearise", "Evenly space the label values to get full use of colourmap.", false);
    ///Colourmaps
    SwitchArg jetArg("", "jet", "Change colourmap of the scalars to the Jet map (default)", false);
    SwitchArg vtkArg("", "vtk", "Change colourmap of the scalars to blue-red (rainbow VTK) map", false);
    SwitchArg rainbowArg("", "rainbow", "Change colourmap of the scalars to the rainbow map", false);
    SwitchArg nihArg("", "NIH", "Change colourmap of the scalars to NIH", false);
    SwitchArg fireArg("", "NIH_FIRE", "Change colourmap of the scalars to NIH Fire", false);
    SwitchArg aalArg("", "AAL", "Change colourmap of the scalars to AAL", false);
    SwitchArg hotArg("", "HOT", "Change colourmap of the scalars to HOT", false);

    ///Add arguments
    cmd.add( labelledImageArg );
    cmd.add( outputArg );
    cmd.add( surfaceArg );
    cmd.add( transformArg );
    cmd.add( loadViewFileArg );
    cmd.add( aboveArg );
    cmd.add( belowArg );
    cmd.add( heightArg );
    cmd.add( widthArg );
    cmd.add( linearDivisionArg );
    cmd.add( linearOffsetArg );

    cmd.add( onscreenArg );
    cmd.add( inverseArg );
    cmd.add( whiteArg );
    cmd.add( loadViewArg );
    cmd.add( humanArg );
    cmd.add( processArg );
    cmd.add( volumeArg );
    cmd.add( lineariseArg );

    cmd.add( jetArg );
    cmd.add( vtkArg );
    cmd.add( rainbowArg );
    cmd.add( nihArg );
    cmd.add( fireArg );
    cmd.add( aalArg );
    cmd.add( hotArg );

    ///Parse the argv array.
    cmd.parse( argc, argv );

    ///Save argument values
    const std::string labelsName = labelledImageArg.getValue();
    const std::string screenName = outputArg.getValue();
    const std::string surfaceName = surfaceArg.getValue();
    const std::string transformName = transformArg.getValue();
    const std::string loadViewName = loadViewFileArg.getValue();
    const float aboveValue = aboveArg.getValue();
    const float belowValue = belowArg.getValue();
    const int windowHeight = heightArg.getValue();
    const int windowWidth = widthArg.getValue();
    int N = linearDivisionArg.getValue();
    int offset = linearOffsetArg.getValue();

    //Check arguments
    //Most of the checking is done by TCLAP
    if(inverseArg.isSet())
    {
        if(!transformArg.isSet())
        {
            cerr << "Error in arguments! Inverse argument needs to be used with the transform argument." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    //Read transform
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New(); //transform matrix
        transform->Identity();
        transform->PostMultiply();
    if(transformArg.isSet())
    {
        vtkSmartPointer<vtkMatrix4x4> matrix = milx::File::OpenITKTransform(transformName);
        transform->Concatenate(matrix);
        if(inverseArg.isSet())
            transform->Inverse();
    }

    std::cerr << "Reading Header for type info etc." << std::endl;
    std::string pixelType, componentType;
    size_t dimensions = 3;
    if( !milx::File::ReadImageInformation(labelsName, pixelType, componentType, dimensions) )
    {
        milx::PrintError("Failed Reading Image. Check the image type/file. Exiting.");
        exit(EXIT_FAILURE);
    }
    milx::PrintInfo("Pixel Type: " + pixelType);
    milx::PrintInfo("Component Type: " + componentType);
    milx::PrintInfo("Dimensions: " + milx::NumberToString(dimensions));

    //Reader object
    QScopedPointer<milxQtFile> reader(new milxQtFile);
    QScopedPointer<milxQtImage> labelledImage(new milxQtImage);  //smart deletion
    bool success = reader->openImage(labelsName.c_str(), labelledImage.data());
        labelledImage->setName(labelsName.c_str());
    if(aboveArg.isSet() || belowArg.isSet())
        labelledImage->threshold(0, belowValue, aboveValue);
        labelledImage->generateImage();

    ///Convert image to 8-bit image if not already so
    QPointer<milxQtImage> eightBitImage;  //smart deletion
    if(componentType == "unsigned_char" || componentType == "unsigned char")
    {
        eightBitImage = labelledImage.data();
        milx::PrintInfo("Detected labelled images.");
    }
    else
    {
        milx::PrintWarning("8-bit image not found and will be cast to 8-bit image for visualisation.");
        eightBitImage = new milxQtImage(labelledImage.data());
        eightBitImage->setData( milx::Image<floatImageType>::CastImage<charImageType>(labelledImage->GetFloatImage()) );
    }
    eightBitImage->setName(labelsName.c_str());
    eightBitImage->generateImage();

    std::vector<unsigned char> values = milx::Image<charImageType>::LabelValues(eightBitImage->GetCharImage());
    std::cout << ">> " << values.size() << " Labels present: " << std::endl;
    for(size_t j = 0; j < values.size(); j ++)
        std::cout << static_cast<unsigned>(values[j]) << ", ";
    std::cout << std::endl;

    if(!linearDivisionArg.isSet())
        N = values.size();

    ///Surface based
    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
        lookupTable->SetTableRange(0.0, values.size()+1);
        lookupTable->Build();

    std::vector< QSharedPointer<milxQtModel> > isoSurfaces;
    QSharedPointer<milxQtRenderWindow> model(new milxQtRenderWindow);
        model->setName("Label Surface Rendering");
        model->SetSize(windowHeight, windowWidth);
        model->generateRender();

if(volumeArg.isSet())
{
//    float minVolumeValue = numeric_limits<float>::max();
//    float maxVolumeValue = numeric_limits<float>::min();
    std::vector<float> colourValuesAsFloat, opacityValuesAsFloat;
    opacityValuesAsFloat.push_back(0.0); //background transparent
    colourValuesAsFloat.push_back(0.0); //background
    for(size_t j = 0; j < values.size(); j ++)
    {
        opacityValuesAsFloat.push_back(1.0);
        colourValuesAsFloat.push_back( static_cast<float>(values[j]) );
    }

    vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    #if VTK_MAJOR_VERSION <= 5
        imageReorient->SetInput(eightBitImage->GetOutput());
    #else
        imageReorient->SetInputData(eightBitImage->GetOutput());
    #endif
        imageReorient->SetFilteredAxis(1);
        imageReorient->FlipAboutOriginOn();
        imageReorient->Update(); //ITK image would have been flipped

    QScopedPointer<milxQtPlot> label(new milxQtPlot);  //smart deletion
        label->SetSource(imageReorient->GetOutput(), true);
        label->setOpacityTransferValues(opacityValuesAsFloat);
        label->setColourTransferValues(colourValuesAsFloat);
    //Colour maps
        label->colourMapToJet(belowValue, aboveValue); //default
    if(vtkArg.isSet())
        label->colourMapToVTK(belowValue, aboveValue);
    if(rainbowArg.isSet())
        label->colourMapToRainbow(belowValue, aboveValue);
    if(nihArg.isSet())
        label->colourMapToNIH(belowValue, aboveValue);
    if(fireArg.isSet())
        label->colourMapToNIH_Fire(belowValue, aboveValue);
    if(aalArg.isSet())
        label->colourMapToAAL(belowValue, aboveValue);
    if(hotArg.isSet())
        label->colourMapToHOT(belowValue, aboveValue);
        label->volumePlot(imageReorient->GetOutput(), true, true);

    vtkVolumeCollection *volumeCollection = label->GetVolumes();

    volumeCollection->InitTraversal();
    vtkVolume * volume = volumeCollection->GetNextVolume();
    if(transformArg.isSet())
        volume->SetUserTransform(transform);
    model->AddVolume(volume);
}
else
{
    const float linearScale = (aboveValue-belowValue)/N;
    for(size_t j = 0; j < values.size(); j ++)
    {
        if(lineariseArg.isSet())
            std::cout << ">> Iso surface label " << static_cast<unsigned>(values[j]) << " -> " << (j+offset)*linearScale << std::endl;
        else
            std::cout << ">> Iso surface label " << static_cast<unsigned>(values[j]) << std::endl;

        QScopedPointer<milxQtImage> label(new milxQtImage);  //smart deletion
            label->SetInput(eightBitImage->GetCharImage());
            label->binaryThreshold(1.0, values[j], values[j]);
            label->generateImage();

        QSharedPointer<milxQtModel> mdl(new milxQtModel);
            mdl->generateIsoSurface(label->GetOutput(), 0, 0.5);
        if(processArg.isSet())
        {
            mdl->quadricDecimate(0.5);
            mdl->smoothSinc(15);
        }
        if(transformArg.isSet())
            mdl->SetTransform(transform);
            mdl->generateModel();

        //Colour maps
        model->colourMapToJet(belowValue, aboveValue); //default
    if(vtkArg.isSet())
        model->colourMapToVTK(belowValue, aboveValue);
    if(rainbowArg.isSet())
        model->colourMapToRainbow(belowValue, aboveValue);
    if(nihArg.isSet())
        model->colourMapToNIH(belowValue, aboveValue);
    if(fireArg.isSet())
        model->colourMapToNIH_Fire(belowValue, aboveValue);
    if(aalArg.isSet())
        model->colourMapToAAL(belowValue, aboveValue);
    if(hotArg.isSet())
        model->colourMapToHOT(belowValue, aboveValue);

        ///Get Colour
        double colour[3];
        if(lineariseArg.isSet())
            model->GetLookupTable()->GetColor((j+offset)*linearScale, colour); //!< Pull colour for data
        else
            model->GetLookupTable()->GetColor(values[j], colour); //!< Pull colour for data
        mdl->changeColour(colour[0], colour[1], colour[2]);

        model->AddActor(mdl->GetActor());
        isoSurfaces.push_back(mdl); //prevent deletion
    }
}
    model->generateRender();

    ///Overlay surface if requested
    QScopedPointer<milxQtModel> surface(new milxQtModel); //smart deletion
    if(surfaceArg.isSet())
    {
        success = reader->openModel(surfaceName.c_str(), surface.data());
        if(!success)
        {
            milx::PrintError("Failed Reading surface. Check the surface type/file. Exiting.");
            exit(EXIT_FAILURE);
        }

        cout << ">> Overlaying surface" << std::endl;
        surface->setName(surfaceName.c_str());
        surface->generateModel();

        model->AddActor(surface->GetActor());
    }
    model->generateRender();

    mainWindow.setCentralWidget(model.data());
    mainWindow.resize(windowWidth, windowHeight);
    if(humanArg.isSet())
        model->disableOrient();
    if(whiteArg.isSet())
        model->background(true);
    if(loadViewArg.isSet())
        model->loadView();
    if(loadViewFileArg.isSet())
        model->loadView(loadViewName.c_str());
    cout << ">> Overlay: Rendering" << std::endl;
    if(!onscreenArg.isSet())
        model->OffScreenRenderingOn();
    else
        mainWindow.show();

    //Update view
    model->GetRenderWindow()->Render();
    qApp->processEvents();

    ///Take screenshot
    vtkSmartPointer<vtkWindowToImageFilter> windowToImage = vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImage->SetInput(model->GetRenderWindow());
        windowToImage->Update();

    QScopedPointer<milxQtFile> writer(new milxQtFile); //Smart deletion
    model->GetRenderWindow()->Render();
    writer->saveImage(screenName.c_str(), windowToImage->GetOutput());
    cout << ">> Complete" << std::endl;

    model->OffScreenRenderingOff(); //Required to prevent double-free
    if(!onscreenArg.isSet())
        return EXIT_SUCCESS;
    else
        return app.exec();
}
