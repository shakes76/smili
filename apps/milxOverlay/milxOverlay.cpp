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
    \brief This program overlays image slices in the model view for comparison.
    \author Shekhar S. Chandra, 2014

    Usage:
    Basic mesh display with scalar bar
    ./bin/milxOverlay 'Atlas_MRI_surface_R_combined.vtk' --scalarbar "Scalars" --coronal --onscreen --jet

    Overlay surface (masked scalars) with image slice and labelled overlay and vector field.
    ./milxOverlay case_dess_9332085_boneSurfaceWithContrastMap_atlas_0.vtk -i case_dess_9332085_preprocess.nii.gz --scalarmask case_dess_9332085_boneSurfaceRelaxed_atlas_0.vtk --vectors case_dess_9332085_boneSurfaceWithThicknessMap_atlas_0.vtk -o knee.png --onscreen --label case_dess_9332085_boneSegmentation_0003_atlas_0.nii.gz
*/
//Qt
#include <QApplication>

#include "milxQtModel.h"
#include "milxQtImage.h"
#include "milxQtFile.h"

#include <vtkMath.h>
#include <vtkWindowToImageFilter.h>

#include <tclap/CmdLine.h> //Command line parser library

using namespace TCLAP;

int main(int argc, char* argv[])
{
    QApplication app(argc,argv);
    QMainWindow mainWindow;//app takes ownership so need to be on stack, not heap

    //---------------------------
    ///Program Info
    milx::PrintInfo("--------------------------------------------------------");
    milx::PrintInfo("SMILI Overlay Tool for Medical Imaging");
    milx::PrintInfo("(c) Copyright Chandra et al., 2015.");
    milx::PrintInfo("University of Queensland, Australia.");
    milx::PrintInfo("Australian e-Health Research Centre, CSIRO, Australia.");
    milx::PrintInfo("SMILI Version: " + milx::NumberToString(milx::Version));
    milx::PrintInfo("milxQt Version: " + milx::NumberToString(milxQtVersion));
    milx::PrintInfo("Application Version: 1.03");
    milx::PrintInfo("--------------------------------------------------------\n");

    //---------------------------
    ///Process Arguments
    CmdLine cmd("An overlay tool for models/images", ' ', milx::NumberToString(milx::Version));

    ///Optional
    ValueArg<std::string> outputArg("o", "output", "Output Screenshot name", false, "screenshot.png", "Output");
    ValueArg<std::string> imageArg("i", "image", "Load and generate image slices as planes (axial, coronal & saggital)", false, "image.nii.gz", "Image");
    ValueArg<std::string> vectorsArg("v", "vectors", "Generate and display vector field from model. If no vectors are present then scalars and normals will be used", false, "vectors.vtk", "Vectors");
    ValueArg<std::string> maskArg("s", "scalarmask", "Clip the main mesh based on binary scalar values in model", false, "mask.vtk", "Scalar Mask");
    ValueArg<std::string> setScalarsArg("", "setscalars", "Set active scalars of model as array name given", false, "Stats", "Set Scalars");
    ValueArg<std::string> labelArg("l", "label", "Overlay binary image/label on the image", false, "labelling.nii.gz", "Label");
    ValueArg<std::string> transformArg("t", "transform", "Transform (ITK Format) to apply to objects being rendered", false, "rigid.txt", "Transform");
    ValueArg<std::string> isoArg("", "isosurface", "Generate Iso-surface (surface from label image) from the image", false, "data.nii.gz", "Isosurface");
    ValueArg<std::string> loadScalarsArg("", "loadscalars", "Load scalars for the mesh from another mesh", false, "scalars.vtk", "Scalars");
    ValueArg<std::string> loadViewFileArg("", "loadviewfile", "Load saved view from file (use onscreen mode to save view files)", false, "camera.cam", "Load View File");
    ValueArg<std::string> scalarBarArg("", "scalarbar", "Enable the scalar bar for display with given name.", false, "Scalars", "Scalar Bar");
    ValueArg<float> opacityArg("", "opacity", "Set the opacity of the models", false, 0.5, "Opacity");
    ValueArg<float> isoValueArg("", "isovalue", "Set the label or iso-surface value for option --isosurface.", false, 0.5, "Isovalue");
    ValueArg<float> minArg("", "min", "Set the minimum value for scalars on model.", false, 0.0, "Minimum Scalar");
    ValueArg<float> maxArg("", "max", "Set the maximum value for scalars on model.", false, 1.0, "Maximum Scalar");
    ValueArg<float> redArg("", "redbackground", "Set the redness value (0-1) for the background.", false, 1.0, "Red Component");
    ValueArg<float> greenArg("", "greenbackground", "Set the greenness value (0-1) for the background.", false, 1.0, "Green Component");
    ValueArg<float> blueArg("", "bluebackground", "Set the blueness value (0-1) for the background.", false, 1.0, "Blue Component");
    ValueArg<float> zoomArg("", "zoom", "Zoom in on the current view by given factor.", false, 2.0, "Zoom");
    ValueArg<int> heightArg("y", "height", "Set the height of the window.", false, 600, "Height");
    ValueArg<int> widthArg("x", "width", "Set the width of the window.", false, 800, "Width");
    ValueArg<int> axialSliceArg("", "axialslice", "Set the axial slice of the image.", false, 0, "Axial Slice");
    ValueArg<int> coronalSliceArg("", "coronalslice", "Set the coronal slice of the image.", false, 0, "Coronal Slice");
    ValueArg<int> saggitalSliceArg("", "saggitalslice", "Set the saggital slice of the image.", false, 0, "Saggital Slice");
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
    SwitchArg loadViewArg("", "loadview", "Load saved view (use smilx or onscreen render mode to view and save with Right Click->View->Save View", false);
    SwitchArg onscreenArg("", "onscreen", "Enable on screen rendering, i.e. display the rendering as an interactive window.", false);
    SwitchArg whiteArg("", "white", "Make background white rather than default gradient colour.", false);
    SwitchArg inverseArg("", "inverse", "Use the inverse transform when transform file is provided.", false);
    SwitchArg noaxialArg("", "noaxial", "Do not show axial slice when image is given.", false);
    SwitchArg nocoronalArg("", "nocoronal", "Do not show coronal slice when image is given.", false);
    SwitchArg nosaggitalArg("", "nosaggital", "Do not show saggital slice when image is given.", false);
    SwitchArg axialArg("", "axial", "Show axial view based on position of first surface.", false);
    SwitchArg coronalArg("", "coronal", "Show coronal view based on position of first surface.", false);
    SwitchArg saggitalArg("", "saggital", "Show saggital view based on position of first surface.", false);
    SwitchArg removeScalarsArg("", "removescalars", "Remove the scalars of the models.", false);
    SwitchArg clampScalarsArg("", "clampscalars", "Clamp the scalars of the models to background colour (outside colourmap) from min/max values provided.", false);
    SwitchArg clampScalarsBelowArg("", "clampbelow", "Clamp the scalars of the models to background colour (outside colourmap) below min values provided.", false);
    SwitchArg outlineArg("", "outline", "Display outline box for the model. Model outline overided by image if present.", false);
    SwitchArg slicesArg("", "noslices", "Prevent any slices from being shown. Other elements from image option etc. will still be displayed.", false);
    SwitchArg autoColourArg("", "autocolour", "Colour each additional model automatically based on order.", false);
    SwitchArg wireframeArg("", "wireframe", "Display initial surface as a wireframe model.", false);
    SwitchArg wireframesArg("", "wireframes", "Display all surfaces as a wireframe models.", false);
    SwitchArg humanArg("", "nohuman", "Disable human orientation glyph.", false);
    SwitchArg specularArg("", "nospecular", "Disable specular lighting for the rendered view.", false);

    ///Mandatory
    UnlabeledMultiArg<std::string> multinames("surfaces", "Surfaces to overlay", true, "Surfaces");

    ///XOR args
    std::vector<Arg*> xorlist;
    xorlist.push_back(&axialArg);
    xorlist.push_back(&coronalArg);
    xorlist.push_back(&saggitalArg);
    cmd.xorAdd(xorlist);

    ///Add argumnets
    cmd.add( multinames );
    cmd.add( outputArg );
    cmd.add( imageArg );
    cmd.add( vectorsArg );
    cmd.add( maskArg );
#ifdef ITK_REVIEW //Review only members
    cmd.add( labelArg );
#endif
    cmd.add( transformArg );
    cmd.add( isoArg );
    cmd.add( minArg );
    cmd.add( maxArg );
    cmd.add( redArg );
    cmd.add( greenArg );
    cmd.add( blueArg );
    cmd.add( zoomArg );
    cmd.add( opacityArg );
    cmd.add( isoValueArg );
    cmd.add( loadScalarsArg );
    cmd.add( setScalarsArg );
    cmd.add( heightArg );
    cmd.add( widthArg );
    cmd.add( axialSliceArg );
    cmd.add( coronalSliceArg );
    cmd.add( saggitalSliceArg );

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
    cmd.add( loadViewArg );
    cmd.add( loadViewFileArg );
    cmd.add( scalarBarArg );
    cmd.add( onscreenArg );
    cmd.add( whiteArg );
    cmd.add( inverseArg );
    cmd.add( noaxialArg );
    cmd.add( nocoronalArg );
    cmd.add( nosaggitalArg );
    cmd.add( removeScalarsArg );
    cmd.add( clampScalarsArg );
    cmd.add( clampScalarsBelowArg );
    cmd.add( outlineArg );
    cmd.add( slicesArg );
    cmd.add( autoColourArg );
    cmd.add( wireframeArg );
    cmd.add( wireframesArg );
    cmd.add( humanArg );
    cmd.add( specularArg );

    ///Parse the argv array.
    cmd.parse( argc, argv );

    ///Get the value parsed by each arg.
    //Filenames of surfaces
    std::vector<std::string> filenames = multinames.getValue();
    const std::string screenName = outputArg.getValue();
    const std::string imageName = imageArg.getValue();
    const std::string vectorsName = vectorsArg.getValue();
    const std::string scalarMaskName = maskArg.getValue();
    const std::string labelName = labelArg.getValue();
    const std::string transformName = transformArg.getValue();
    const std::string isoName = isoArg.getValue();
    const std::string loadScalarsName = loadScalarsArg.getValue();
    const std::string loadViewName = loadViewFileArg.getValue();
    const std::string arrayName = setScalarsArg.getValue();
    const std::string barName = scalarBarArg.getValue();
    const float opacity = opacityArg.getValue();
    const float isoValue = isoValueArg.getValue();
    const float minValue = minArg.getValue();
    const float maxValue = maxArg.getValue();
    const float redValue = redArg.getValue();
    const float greenValue = greenArg.getValue();
    const float blueValue = blueArg.getValue();
    const float zoomFactor = zoomArg.getValue();
    const int windowHeight = heightArg.getValue();
    const int windowWidth = widthArg.getValue();
    const int axialSliceNumber = axialSliceArg.getValue();
    const int coronalSliceNumber = coronalSliceArg.getValue();
    const int saggitalSliceNumber = saggitalSliceArg.getValue();
//    const float scaleFactor = scaleArg.getValue();

    //Check arguments
    //Most of the checking is done by TCLAP
    if(inverseArg.isSet())
    {
        if(!transformArg.isSet())
        {
            cerr << "Error in arguments! Inverse argument needs to be used with the transform argument." << endl;
            exit(EXIT_FAILURE);
        }
    }
    if(noaxialArg.isSet() || nosaggitalArg.isSet() || nocoronalArg.isSet() || axialSliceArg.isSet() || coronalSliceArg.isSet() || saggitalSliceArg.isSet())
    {
        if(!imageArg.isSet())
        {
            cerr << "Error in arguments! View/Slice arguments need to be used with the image argument." << endl;
            exit(EXIT_FAILURE);
        }
    }
    if(isoValueArg.isSet())
    {
        if(!isoArg.isSet())
        {
            cerr << "Error in arguments! Isovalue argument needs to be used with the isosurface argument." << endl;
            exit(EXIT_FAILURE);
        }
    }

    ///Read
    //Reader object
    QScopedPointer<milxQtFile> reader(new milxQtFile);
    bool errorReading = false;
    bool success = false;

    //Read model
    cout << ">> Overlay: Reading Models" << endl;
    vtkSmartPointer<vtkPolyDataCollection> collection;
    success = milx::File::OpenModelCollection(filenames, collection);
    if(!success) //Error printed inside
    {
        cerr << "Error reading models!" << endl;
        exit(EXIT_FAILURE);
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

    //Load Models
    QList< QSharedPointer<milxQtModel> > models; //Use shared pointer here, more overhead but allows assignment etc.
    QSharedPointer<milxQtModel> model;
    const size_t n = collection->GetNumberOfItems();

    if(n < 1)
    {
        cerr << "At least one model must be provided!" << endl;
        exit(EXIT_FAILURE);
    }

    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    if(autoColourArg.isSet())
    {
        lookupTable->SetTableRange(0.0, n+1);
        lookupTable->Build();
    }

    double range[2];
    collection->InitTraversal();
    for(size_t j = 0; j < n; j ++)
    {
        double colour[3];

        QSharedPointer<milxQtModel> mdl(new milxQtModel);
          mdl->setName(filenames[j].c_str());
          mdl->SetInput(collection->GetNextItem()); //!< create model

        if(removeScalarsArg.isSet() && mdl->GetScalars())
            mdl->RemoveScalars();
        if(loadScalarsArg.isSet())
            mdl->loadScalars(loadScalarsName.c_str());

        mdl->generateModel();

        //Colouring
        if(autoColourArg.isSet() && j > 0)
        {
            ///Get Colour
            lookupTable->GetColor(j, colour); //!< Pull colour for data
            mdl->changeColour(colour[0], colour[1], colour[2]);
        }
        if(setScalarsArg.isSet())
            mdl->showArray(arrayName.c_str());

        ///Setup scalar ranges
        if(mdl->GetScalars())
        {
            mdl->GetScalarRange(range);
            if(clampScalarsBelowArg.isSet())
                mdl->thresholdScalars(minValue, range[1], vtkMath::Nan());
            if(clampScalarsArg.isSet())
                mdl->thresholdScalars(minValue, maxValue, vtkMath::Nan());
            if(minArg.isSet())
                range[0] = minValue;
            if(maxArg.isSet())
                range[1] = maxValue;
        }

        if( opacityArg.isSet() && (!autoColourArg.isSet() || j == 0) ) //needs to be after generate Model
            mdl->SetOpacity(opacity);
        if(wireframesArg.isSet())
            mdl->generateWireframe();

        models.append(mdl);
    }
    model = models[0];
    mainWindow.setCentralWidget(model.data());
    std::cerr << "Done" << std::endl;

    //Read vectors model
    QScopedPointer<milxQtModel> modelVectors(new milxQtModel); //smart deletion
    if(vectorsArg.isSet())
    {
      success = reader->openModel(vectorsName.c_str(), modelVectors.data());
      if(success)
      {
          cout << ">> Applying Vectors" << endl;
          modelVectors->setName(vectorsName.c_str());
          modelVectors->generateModel();
          modelVectors->generateVectorField();
//          if(transformArg.isSet())
//              modelVectors->SetTransform(transform);
      }
      else
          errorReading = true;
    }
    //Read mask model
    QScopedPointer<milxQtModel> modelMask(new milxQtModel); //smart deletion
    if(maskArg.isSet())
    {
      success = reader->openModel(scalarMaskName.c_str(), modelMask.data());
      if(success)
      {
          cout << ">> Overlay: Applying Mask" << endl;
          modelMask->setName(scalarMaskName.c_str());
          modelMask->generateModel();

          //Clipping and scalars stuff
          modelMask->GetScalars()->SetName("Mask");
          model->GetScalars()->SetName("Scalars");
          model->GetOutput()->GetPointData()->AddArray(modelMask->GetScalars());
          model->showArray("Mask");
          model->threshold(1.0, 1.0);
          model->showArray("Scalars");

          //Reset scalar range for display
          model->resetScalarRange();

          if(model->GetNumberOfPoints() == 0)
          {
              cerr << "Error using scalar mask. Model no longer has any points!" << endl;
              exit(EXIT_FAILURE);
          }
      }
      else
          errorReading = true;
    }
    //Read iso surface
    QScopedPointer<milxQtModel> modelIso(new milxQtModel); //smart deletion
    QScopedPointer<milxQtImage> imgIso(new milxQtImage);  //smart deletion
    if(isoArg.isSet())
    {
        cout << ">> Overlay: Applying Isosurface" << endl;
        success = reader->openImage(isoName.c_str(), imgIso.data());
            imgIso->setName(isoName.c_str());
            imgIso->generateImage();
            imgIso->threshold(0, isoValue-1.0, isoValue);// fix thresholding here

        if(!success)
            errorReading = true;

//        modelIso->generatePolyDataFromImage(imgIso->GetOutput());
        modelIso->generateIsoSurface(imgIso->GetOutput(), 0, isoValue-1.0);// fix thresholding here
        modelIso->generateModel();
        modelIso->changeColour(0,1,0);
//        if(transformArg.isSet())
//              modelIso->SetTransform(transform);
    }

    if(errorReading)
    {
        cerr << "Error Reading one or more of the input files. Exiting." << endl;
        exit(EXIT_FAILURE);
    }

    //Read Image
    QScopedPointer<milxQtImage> img(new milxQtImage);  //hierarchy deletion
    QScopedPointer<milxQtImage> imgCoronal(new milxQtImage);  //hierarchy deletion
    QScopedPointer<milxQtImage> imgSagittal(new milxQtImage);  //hierarchy deletion
    vtkSmartPointer<vtkMatrix4x4> orientTransform = vtkSmartPointer<vtkMatrix4x4>::New();

    //Read image
    vtkSmartPointer<vtkTransform> transform2 = vtkSmartPointer<vtkTransform>::New(); //transform matrix
        transform2->Identity();
        transform2->PostMultiply();
    if(imageArg.isSet())
    {
        cout << ">> Overlay: Reading Image" << endl;
        errorReading = false;
        success = reader->openImage(imageName.c_str(), img.data());

        ///Generate Views
        if(success)
        {
            //Generate the images
            //Axial default view atm
            img->setName(imageName.c_str());
            img->generateImage();
            if(axialSliceArg.isSet())
                img->setSlice(axialSliceNumber);
    //        if(indexNIH > 0)
    //            img->colourMapToNIH();
    //        if(indexNIH_FIRE > 0)
    //            img->colourMapToNIH_Fire();
    //        if(indexAAL > 0)
    //            img->colourMapToAAL();

#ifdef ITK_REVIEW //Review only members
            if(labelArg.isSet())
            {
                img->overlayContour(labelName.c_str()); //Generates an RGB image, from VTK class
            }
#endif

//            img->flip(false, true, false);

            //Share to create other views
			imgCoronal->setDisplayData(img.data());
			imgSagittal->setDisplayData(img.data());

            //ATM in API, fastest coding way is to load image for each view desired
            //Thus need to copy image, memory usage is high here
            //Setup coronal view
            if(!nocoronalArg.isSet())
            {
                imgCoronal->setName(imageName.c_str());
                imgCoronal->generateImage();
                imgCoronal->viewToCoronal();
                if(coronalSliceArg.isSet())
                    imgCoronal->setSlice(coronalSliceNumber);
            }

            //Setup sagittal view
            if(!nosaggitalArg.isSet())
            {
                imgSagittal->setName(imageName.c_str());
                imgSagittal->generateImage();
                imgSagittal->viewToSagittal();
                if(saggitalSliceArg.isSet())
                    imgSagittal->setSlice(saggitalSliceNumber);
            }
        }
        else
            errorReading = true;

        if(errorReading)
        {
            cerr << "Error Reading the image file. Exiting." << endl;
            exit(EXIT_FAILURE);
        }

        ///Setup orientation transform matrices
        //Flip y to make same coordinates as model
        orientTransform->DeepCopy(img->getTransformMatrix());
        orientTransform->Invert();
        transform2->Concatenate(orientTransform);
//        transform2->Concatenate(transform->GetMatrix());
        cout << ">> Overlay: Transforming Actors" << endl;
    }

    ///Display
    model->setWindowTitle("MILX-Overlay");
    model->SetSize(windowHeight, windowWidth);
    mainWindow.resize(windowWidth, windowHeight);
    //Add models
    foreach(QSharedPointer<milxQtModel> mdl, models)
    {
        model->addModelActor(mdl->GetActor());
    }
    //Add image slices
    if(imageArg.isSet() && !slicesArg.isSet())
    {
        if(!noaxialArg.isSet())
            model->addImageActor(img->GetImageActor(), transform2->GetMatrix());
        if(!nocoronalArg.isSet())
            model->addImageActor(imgCoronal->GetImageActor(), transform2->GetMatrix());
        if(!nosaggitalArg.isSet())
            model->addImageActor(imgSagittal->GetImageActor(), transform2->GetMatrix());

        ///Here we use the Qt signals and slots directly as it was found that the VTK-Qt connector caused problems
        ///with the image actors.
        QObject::connect(img.data(), SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), model.data(), SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        QObject::connect(imgCoronal.data(), SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), model.data(), SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        QObject::connect(imgSagittal.data(), SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), model.data(), SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
    }

    //Colour maps
    cout << ">>> Overlay: Setting Colourmap" << endl;
    if(jetArg.isSet())
        model->colourMapToJet();
    if(vtkArg.isSet())
        model->colourMapToVTK();
    if(hsvArg.isSet())
      model->colourMapToHSV();
    if(rainbowArg.isSet())
        model->colourMapToRainbow();
    if(spectralArg.isSet())
      model->colourMapToSpectral();
    if(nihArg.isSet())
        model->colourMapToNIH();
    if(fireArg.isSet())
        model->colourMapToNIH_Fire();
    if(aalArg.isSet())
        model->colourMapToAAL();
    if(fsArg.isSet())
      model->colourMapToFS();
    if(hotArg.isSet())
        model->colourMapToHOT();
    if(coolArg.isSet())
      model->colourMapToCOOL();
    if(minArg.isSet() || maxArg.isSet()) //set the range if another colour map is used
        model->SetScalarRange(range);
    if(scalarBarArg.isSet())
    {
        if(minArg.isSet() || maxArg.isSet())
            model->enableScale(barName.c_str(), true, range[0], range[1]);
        else
            model->enableScale(barName.c_str(), true);
    }

    if(vectorsArg.isSet())
        model->AddActor(modelVectors->GetActor());
    if(isoArg.isSet())
    {
//        if(onscreenArg.isSet())
//            modelIso->show();
        model->AddActor(modelIso->GetActor());
    }
    if(whiteArg.isSet())
        model->background(true);
    if(outlineArg.isSet())
    {
        if(imageArg.isSet())
            model->enableOutline(img->GetOutput());
        else
            model->enableOutline();
    }

    //setup background
    double backgroundColours[3];
    model->GetBackground(backgroundColours[0], backgroundColours[1], backgroundColours[2]);
    if(humanArg.isSet())
        model->disableOrient();
    if(redArg.isSet())
        backgroundColours[0] = redValue;
    if(greenArg.isSet())
        backgroundColours[1] = greenValue;
    if(blueArg.isSet())
        backgroundColours[2] = blueValue;
    if(redArg.isSet() || blueArg.isSet() || greenArg.isSet())
        model->SetBackground(backgroundColours[0], backgroundColours[1], backgroundColours[2]);

    //Apply transform to main model last to ensure sub models are correctly transform also
    if(wireframeArg.isSet())
        model->generateWireframe();
    if(transformArg.isSet())
        model->SetTransform(transform);
    if(axialArg.isSet())
        model->viewToAxial();
    if(coronalArg.isSet())
        model->viewToCoronal();
    if(saggitalArg.isSet())
        model->viewToSagittal();
    if(loadViewArg.isSet())
        model->loadView();
    if(loadViewFileArg.isSet())
        model->loadView(loadViewName.c_str());
    if(specularArg.isSet())
        model->disableSpecularDisplay();
    //Zoom
    if(zoomArg.isSet())
    {
      vtkCamera *camera = model->GetRenderer()->GetActiveCamera();
      camera->Zoom(zoomFactor);
    }

    cout << ">> Overlay: Rendering" << endl;
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
    cout << ">> Complete" << endl;

    model->OffScreenRenderingOff(); //Required to prevent double-free
    if(!onscreenArg.isSet())
        return EXIT_SUCCESS;
    else
        return app.exec();
}
