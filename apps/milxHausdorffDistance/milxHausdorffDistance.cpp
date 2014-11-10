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
    \brief This program computes the Hausdorff distance of surfaces given labels or surfaces.
    \author Shekhar S. Chandra, 2014

    Usage:

*/
//ITK
#include <itkImage.h>
//VTK
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
//SMILI
#include "milxDeformableModel.h"
#include "milxImage.h"
#include "milxFile.h"

#include <tclap/CmdLine.h> //Command line parser library

using namespace TCLAP;

typedef float FloatPixelType;
typedef itk::Image<FloatPixelType, 3> FloatImageType;
typedef unsigned char InputPixelType;
typedef itk::Image<InputPixelType, 3> LabelImageType;
typedef LabelImageType::SizeValueType SizeValueType;

int main(int argc, char* argv[])
{
    //---------------------------
    ///Program Info
    milx::PrintInfo("--------------------------------------------------------");
    milx::PrintInfo("MILX-SMILI Hausdorff Distance Tool for Medical Imaging");
    milx::PrintInfo("(c) Copyright CSIRO, 2012.");
    milx::PrintInfo("Australian e-Health Research Centre, CSIRO.");
    milx::PrintInfo("SMILI Version: " + milx::NumberToString(milx::Version));
    milx::PrintInfo("Application Version: 1.00");
    milx::PrintInfo("Processors to be used: " + milx::NumberToString(milx::NumberOfProcessors()/2));
    milx::PrintInfo("--------------------------------------------------------\n");

    //---------------------------
    ///Process Arguments
    CmdLine cmd("A Hausdorff Distance tool for models", ' ', milx::NumberToString(milx::Version));

    ///Optional
    ValueArg<std::string> outputArg("o", "output", "Output model name", false, "hausdorff.vtk", "Output");
    ValueArg<std::string> prefixArg("p", "prefix", "Output prefix for multiple output", false, "hausdorff_", "Output Prefix");
    ValueArg<std::string> labelArg("l", "label", "Compute the distances from the labelled image to the surface(s) provided.", true, "labelling.nii.gz", "Label");
    ValueArg<float> labelValueArg("", "labelvalue", "Set the label value for option --label.", true, 255, "Label Value");
    ValueArg<size_t> caseArg("c", "case", "Set the case ID being done. Used to name extra output.", false, 0, "Case");
    ///Switches
    SwitchArg symmetricArg("s", "symmetric", "Compute forward and backward distances. This is required to get Hausdorff distance.", false);

    ///Mandatory
    UnlabeledMultiArg<std::string> multinames("surfaces", "Surfaces to compute the distances with.", true, "Surfaces");

    ///Add argumnets
    cmd.add( multinames );
    cmd.add( outputArg );
    cmd.add( prefixArg );
    cmd.add( labelArg );
    cmd.add( labelValueArg );
    cmd.add( symmetricArg );
    cmd.add( caseArg );

    ///Parse the argv array.
    cmd.parse( argc, argv );

    ///Get the value parsed by each arg.
    //Filenames of surfaces
    std::vector<std::string> filenames = multinames.getValue();
    const std::string outputName = outputArg.getValue();
    const std::string prefixName = prefixArg.getValue();
    const std::string labelName = labelArg.getValue();
    const float labelValue = labelValueArg.getValue();
    const size_t caseID = caseArg.getValue();

    ///Setup ITK Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(milx::NumberOfProcessors()/2);

    //Check arguments
    //Most of the checking is done by TCLAP
    if(labelValueArg.isSet())
    {
        if(!labelArg.isSet())
        {
            cerr << "Error in arguments! Label argument needs to be used with the label value argument." << endl;
            exit(EXIT_FAILURE);
        }
    }

    ///Read
    bool success = false;

    //Load Models
    cout << ">> Hausdorff Distance: Reading Models" << endl;
    vtkSmartPointer<vtkPolyDataCollection> collection;
    success = milx::File::OpenModelCollection(filenames, collection);
    if(!success) //Error printed inside
    {
        cerr << "Error reading models!" << endl;
        exit(EXIT_FAILURE);
    }
    const size_t n = collection->GetNumberOfItems();
    if(n < 1)
    {
        cerr << "At least one model must be provided!" << endl;
        exit(EXIT_FAILURE);
    }
    collection->InitTraversal();
    vtkSmartPointer<vtkPolyData> surface = collection->GetNextItem();
    double bounds[6];
    surface->GetBounds(bounds);
    std::cerr << "Done" << std::endl;

    //Read labels
    itk::SmartPointer<LabelImageType> labelledImage, resizedLabelledImage, thresholdedImage;  //smart deletion
    itk::SmartPointer<FloatImageType> distanceMap;
    const bool binary = false, signedDistance = true, insideDistance = false, squaredDistance = false;
    const unsigned char objectValue = 255;
    if(labelArg.isSet())
    {
        cout << ">> Hausdorff Distance: Using Labelling" << endl;
        cout << "Loading... " << endl;
        success = milx::File::OpenImage<LabelImageType>(labelName, labelledImage);

        cout << "Thresholding... " << endl;
        //~ thresholdedImage = milx::Image<LabelImageType>::BinaryThresholdImage<LabelImageType>(resizedLabelledImage, 0, objectValue, labelValue, labelValue);
        thresholdedImage = milx::Image<LabelImageType>::BinaryThresholdImage<LabelImageType>(labelledImage, 0, objectValue, labelValue, labelValue);

        cout << "Computing Distance Map of Label... " << endl;
        distanceMap = milx::Image<LabelImageType>::DistanceMap<FloatImageType>(thresholdedImage, binary, signedDistance, insideDistance, squaredDistance);

//        success = milx::File::SaveImage<FloatImageType>(prefixName + "_dmap.nii.gz", distanceMap);
    }
    if(!success)
    {
        cerr << "Error Reading one or more of the input files. Exiting." << endl;
        exit(EXIT_FAILURE);
    }

    //Read labels and generate iso surface
    //Clip mesh to ensure same FoV
    cout << ">> Hausdorff Distance: Computing Absolute Surface Distances... " << endl;
    vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
    milx::DeformableModel model(surface);
    model.RemoveScalars(); //remove because mark will not set as 1.0, causing problems with meshes having scalars already
    model.MarkSurfaceInsideImage<FloatImageType>(model.GetOutput(), distanceMap, weights); //sets weights as scalars in here
    model.Clip(1.0, 1.0);

    //Save distances as scalars on mesh
    const bool absValues = true;
    vtkSmartPointer<vtkFloatArray> scalars = model.SurfaceScalarsFromImage<FloatImageType>(model.GetOutput(), distanceMap, absValues);
    scalars->SetName("Surface Distance Errors");
    model.SetScalars(scalars);

    //Restore correspondence to copy scalars over to full (unclipped) mesh, since the mesh is clipped to ensure same FoV
    cout << "Regions of Surface Outside the Image are marked with -1." << endl;
    weights->FillComponent(0, -1.0); //resized within MarkSurfaceInsideImage function
    for(int j = 0; j < model.GetOutput()->GetNumberOfPoints(); j ++)
    {
      vtkIdType index = surface->FindPoint(model.GetOutput()->GetPoint(j));
      weights->SetTuple1(index, scalars->GetTuple1(j));
    }

    milx::DeformableModel symmetricModel;
    if(symmetricArg.isSet())
    {
      if(labelArg.isSet())
      {
        //Voxelise surface
        LabelImageType::SpacingType labelSpacing = labelledImage->GetSpacing();
        double bounds[6];
        model.GetOutput()->GetBounds(bounds);
        vtkSmartPointer<vtkImageData> voxelisedModel = model.Voxelise(objectValue, labelSpacing.GetDataPointer(), bounds);
        LabelImageType::Pointer modelLabel = milx::Image<LabelImageType>::ConvertVTKImageToITKImage(voxelisedModel);
        cout << "Computing Distance Map of Surface... " << endl;
        itk::SmartPointer<FloatImageType> modelDistanceMap = milx::Image<LabelImageType>::DistanceMap<FloatImageType>(modelLabel, binary, signedDistance, insideDistance, squaredDistance);

        //~ //Debug, write distance map
        //~ milx::File::SaveImage(prefixName + "_model_distance_map.nii.gz", modelDistanceMap);

        //Isosurface label
        cout << "Generating Iso Surface... " << endl;
        vtkSmartPointer<vtkImageData> labelledImageVTK = vtkSmartPointer<vtkImageData>::New();
        labelledImageVTK->DeepCopy( milx::Image<LabelImageType>::ConvertITKImageToVTKImage(thresholdedImage) ); //no orientation kept here
        milx::PrintDebug("Apply orientation to VTK Image form of Labelled Image");
        vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New(); //not used but required by ApplyOrientationToVTKImage
        vtkSmartPointer<vtkImageData> restoredLabelledImageVTK = vtkSmartPointer<vtkImageData>::New();
        restoredLabelledImageVTK->DeepCopy( milx::Image<LabelImageType>::ApplyOrientationToVTKImage(labelledImageVTK, labelledImage, matrix, true) );

        //Write result
        //~ cout << "Saving VTK Image Result... " << endl;
        //~ milx::File::SaveImage<LabelImageType>(prefixName + "_vtkimage.vti", restoredLabelledImageVTK);

        //Iso surface
        milx::DeformableModel isoSurface;
        milx::PrintDebug("VTK Marching cubes");
        isoSurface.IsoSurface(restoredLabelledImageVTK, objectValue-1);

        //Write result
//        cout << "Saving Iso Surface Result... " << endl;
//        milx::File::SaveModel(prefixName + milx::NumberToString(caseID) + "_isosurface.vtp", isoSurface.GetOutput());

        cout << ">> Hausdorff Distance: Computing Backward Absolute Surface Distances... " << endl;
        vtkSmartPointer<vtkFloatArray> weights2 = vtkSmartPointer<vtkFloatArray>::New();
        symmetricModel.SetInput(isoSurface.GetOutput());
        symmetricModel.RemoveScalars(); //remove because mark will not set as 1.0, causing problems with meshes having scalars already
        symmetricModel.MarkSurfaceInsideImage<FloatImageType>(symmetricModel.GetOutput(), modelDistanceMap, weights2); //sets weights2 as scalars in here
        isoSurface.SetScalars(weights2);
        //~ cout << "Saving Clipped Iso Surface Result... " << endl;
        //~ milx::File::SaveModel(prefixName + "_isosurface_clipped.vtp", symmetricModel.GetOutput());
        symmetricModel.Clip(1.0, 1.0);

        //Save distances as scalars on mesh
        const bool absValues = true;
        vtkSmartPointer<vtkFloatArray> scalars2 = model.SurfaceScalarsFromImage<FloatImageType>(symmetricModel.GetOutput(), modelDistanceMap, absValues);
        scalars2->SetName("Surface Distance Errors");
        symmetricModel.SetScalars(scalars2);

        //Restore correspondence to copy scalars over to full (unclipped) mesh, since the mesh is clipped to ensure same FoV
        weights2->FillComponent(0, -1.0);
        for(int j = 0; j < symmetricModel.GetOutput()->GetNumberOfPoints(); j ++)
        {
          vtkIdType index = isoSurface.GetOutput()->FindPoint(symmetricModel.GetOutput()->GetPoint(j));
          weights2->SetTuple1(index, scalars2->GetTuple1(j));
        }

        ///Compute Hausdorff
        double range1[2], range2[2];
        surface->GetScalarRange(range1);
        isoSurface.GetOutput()->GetScalarRange(range2);
        std::cout << "Forward Max Distance: " << range1[1] << std::endl;
        std::cout << "Backward Max Distance: " << range2[1] << std::endl;
        std::cout << "Hausdorff Distance: " << milx::Maximum<double>(range1[1], range2[1]) << std::endl;

        //Write result
        cout << "Saving Backward Result... " << endl;
        milx::File::SaveModel(prefixName + milx::NumberToString(caseID) + "_backward.vtp", isoSurface.GetOutput());
      }
    }

    //Copy result over as scalars
    surface->GetPointData()->SetScalars(weights);

    //Write result
    cout << "Saving Result... " << endl;
    milx::File::SaveModel(outputName, surface);

    cout << ">> Hausdorff Distance: Operation Complete" << endl;
    return EXIT_SUCCESS;
}
