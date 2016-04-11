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
#include <vector>
#include <tclap/CmdLine.h> //Command line parser library

// vtk includes
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkMultiThreader.h>
#include <vtkPolyDataCollection.h>
#include <vtkTransformCollection.h>
// SMILI
#include <milxFile.h>
#include <milxModel.h>

using namespace TCLAP;

typedef std::vector< std::string >::iterator stringiterator;

//Supported operations
enum operations {none = 0, convert, duplicate, cat, split, scale, decimate, smooth, laplacian, thresholdscalars, clip, flip, diffscalars, copyscalars, diffscalarspairs, statscalars, removescalars, mse, procrustes, icp};

/**
  \file milxModelApp.cxx
  \brief A swiss-army knife of model/surface processing. This application implements single/batch processing of models with various processing options found in the SMILI library.
  \ingroup Applications

  You are able to smooth, align, decimate etc. a number of surfaces easily via the commandline with this application (use --help). This application uses the milxModel class to do all processing.

  Usage:
  Concatenate surfaces
  \verbatim
  $ milxModelApp *.vtk --cat -o initial.vtk
  \endverbatim
  
  Split surfaces based on knowledge of known components of that surface
  \verbatim
  milxModelApp --split combined_new_weights.vtk --component Atlas_MRI_surface_R_ace.vtk --component Atlas_MRI_surface_R_fem.vtk -p split_ 
  \endverbatim

  Copy scalars from atlas mesh to other meshes in directory and output to another directory with same names
  \verbatim
  milxModelApp --scalarcopy focus_atlases/focus_bladder_atlas.vtk results/bladder/asm_bladder_*.vtk -p results_scalars/bladder/
  \endverbatim

  Clip meshes in directory keeping only parts with value of 1 and output to another directory with same names
  \verbatim
  milxModelApp --clip 1 results_scalars/bladder/asm_bladder_*.vtk -p results_clipped/bladder/
  \endverbatim

  Compute scalars stats (mean, variance etc. per point) of meshes in directory and output single mesh (of first mesh) with stats as arrays in mesh
  \verbatim
  milxModelApp --scalarstats Hausdorff/bladder/bladder__*.vtk -o bladder_stats.vtk
  \endverbatim
*/
int main(int argc, char *argv[])
{
  //---------------------------
  ///Program Info
  milx::PrintInfo("--------------------------------------------------------");
  milx::PrintInfo("SMILI Model Tool for Models/Surfaces/Meshes.");
  milx::PrintInfo("(c) Copyright Chandra et al., 2015.");
  milx::PrintInfo("Version: " + milx::NumberToString(milx::Version));
  milx::PrintInfo("University of Queensland, Australia.");
  milx::PrintInfo("Australian e-Health Research Centre, CSIRO, Australia.");
  milx::PrintInfo("--------------------------------------------------------\n");

  //---------------------------
  ///Process Arguments
  CmdLine cmd("A diagnostic tool for models/surface operations", ' ', milx::NumberToString(milx::Version));

  ///Optional
  ValueArg<size_t> threadsArg("", "threads", "Set he number of global threads to use.", false, milx::NumberOfProcessors(), "Threads");
  ValueArg<std::string> outputArg("o", "output", "Output Surface", false, "result.vtk", "Output");
  ValueArg<std::string> prefixArg("p", "prefix", "Output prefix for multiple output", false, "surface_", "Output Prefix");
  ValueArg<std::string> outputFormatArg("", "outputformat", "Specify the default output format for multiple outputs (vtk, vtp, ply, stl, default: same as input)", false, "vtk", "Output format");
//  MultiArg<std::string> altmultinames("s", "altsurfaces", "Another set of surfaces to be used in operation", false, "AltSurfaces");
  ValueArg<float> decimateArg("d", "decimate", "Decimate all the meshes provided using the Quadric Decimate algorithm with decimation factor.", false, 0.5, "Decimate");
  ValueArg<float> smoothArg("", "smooth", "Smooth all the meshes provided using the Windowed Sinc Algorithm with iterations.", false, 18, "Smooth");
  ValueArg<float> laplacianArg("", "laplacian", "Smooth all the meshes provided using the Laplacian Algorithm with iterations.", false, 18, "Laplacian");
  ValueArg<float> scaleArg("s", "scale", "Scale the coordinates of the points by scale.", false, 0.9, "Scale");
  ValueArg<float> thresholdAboveArg("", "thresholdabove", "Thresold scalars above value.", false, 0.0, "Above");
  ValueArg<float> thresholdBelowArg("", "thresholdbelow", "Thresold scalars below value.", false, 0.0, "Below");
  ValueArg<float> clipArg("", "clip", "Clip model based on scalars value (keeping only parts with value).", false, 1.0, "Clip");
  MultiArg<std::string> componentsArg("", "component", "Surface is a component of the surfaces.", false, "Component");
  ///Clamped Optional
  std::vector<size_t> axesAllowed;
  axesAllowed.push_back(0);
  axesAllowed.push_back(1);
  axesAllowed.push_back(2);
  ValuesConstraint<size_t> allowedAxesVals( axesAllowed );
  ValueArg<size_t> flipArg("f", "flip", "Flip the meshes in the axis provided (0: x-axis, 1: y-axis, 2: z-axis).", false, 0, &allowedAxesVals);
  ///Switches
  SwitchArg partitionArg("", "partitionList", "Partition the list of surfaces provided into two for operating on one vs the other", false);
  ///XOR Switches
  SwitchArg duplicateArg("", "duplicate", "Simply open and save the input file(s). Supports single (-o) or multiple (-p) input.", false);
  SwitchArg convertArg("c", "convert", "Convert the model from the current format to the one given at output (from file extension). Supports single (-o) or multiple (-p) input.", false);
  SwitchArg concatenateArg("", "cat", "Concatenate N surfaces into single mesh with filename provided.", false);
  SwitchArg colourConcatenateArg("", "colourcat", "Concatenate N surfaces into single mesh with filename provided and colour them.", false);
  SwitchArg splitArg("", "split", "Split each surface given components.", false);
  SwitchArg diffScalarArg("", "scalardiff", "Compute the differences in Scalars.", false);
  SwitchArg statsScalarArg("", "scalarstats", "Compute statistics of scalars (mean, variance etc. per point) output mesh with stats as arrays.", false);
  SwitchArg removeScalarArg("", "scalarremove", "Remove the scalars.", false);
  SwitchArg copyScalarArg("", "scalarcopy", "Copy the Scalars from first mesh to all others while removing existing ones.", false);
  SwitchArg mseArg("", "mse", "Mean Squared Error of Points in models.", false);
  SwitchArg procrustesArg("", "procrustes", "Similarity alignment of surfaces assuming points have correspondence. Use --rigid if rigid alignment is required.", false);
  SwitchArg rigidArg("", "rigid", "Rigid alignment of surfaces assuming points have correspondence. Use with procrustes or other registration arguments.", false);
  SwitchArg icpArg("", "icp", "Iterative Closest Points alignment of surfaces assuming points don't have correspondence. Last surface in list is used as the reference 'fixed' surface.", false);
  SwitchArg saveTransformsArg("", "savetransforms", "Save the transformation matrix after surface alignment. For use with --icp", false);

  ///Mandatory
  UnlabeledMultiArg<std::string> multinames("surfaces", "Surfaces to operate on", true, "Surfaces");

  ///Add argumnets
  cmd.add( threadsArg );
  cmd.add( multinames );
//  cmd.add( altmultinames );
  cmd.add( outputArg );
  cmd.add( prefixArg );
  cmd.add( outputFormatArg );
  cmd.add( componentsArg );
  cmd.add( rigidArg );
  cmd.add( partitionArg );
  cmd.add( saveTransformsArg );
//  cmd.add( thresholdAboveArg );
//  cmd.add( thresholdBelowArg );
  ///XOR args
  std::vector<Arg*> xorlist;
  xorlist.push_back(&duplicateArg);
  xorlist.push_back(&convertArg);
  xorlist.push_back(&concatenateArg);
  xorlist.push_back(&colourConcatenateArg);
  xorlist.push_back(&splitArg);
  xorlist.push_back(&scaleArg);
  xorlist.push_back(&decimateArg);
  xorlist.push_back(&smoothArg);
  xorlist.push_back(&laplacianArg);
  xorlist.push_back(&diffScalarArg);
  xorlist.push_back(&statsScalarArg);
  xorlist.push_back(&removeScalarArg);
  xorlist.push_back(&copyScalarArg);
  xorlist.push_back(&mseArg);
  xorlist.push_back(&procrustesArg);
  xorlist.push_back(&icpArg);
  xorlist.push_back(&thresholdAboveArg);
  xorlist.push_back(&thresholdBelowArg);
  xorlist.push_back(&clipArg);
  xorlist.push_back(&flipArg);
  cmd.xorAdd(xorlist);

  ///Parse the argv array.
  cmd.parse( argc, argv );

  ///Get the value parsed by each arg.
  const size_t threads = threadsArg.getValue();
  //Filenames of surfaces
  std::vector<std::string> filenames = multinames.getValue();
//  std::vector<std::string> altfilenames = altmultinames.getValue();
  std::string outputName = outputArg.getValue();
  const std::string prefixName = prefixArg.getValue();
  const float decimateFactor = decimateArg.getValue();
  const float scaleFactor = scaleArg.getValue();
  const float smoothIterations = smoothArg.getValue();
  const float laplacianIterations = laplacianArg.getValue();
  float thresholdAbove = thresholdAboveArg.getValue();
  float thresholdBelow = thresholdBelowArg.getValue();
  float clipValue = clipArg.getValue();
  size_t flipAxis = flipArg.getValue();
  std::vector<std::string> componentNames = componentsArg.getValue();

  ///Setup ITK Threads
  vtkMultiThreader::SetGlobalDefaultNumberOfThreads(threads);
  milx::PrintInfo("Threads to use: " + milx::NumberToString(threads));

  ///Display operation
  operations operation = none;
  bool componentsValid = false;
  if(duplicateArg.isSet())
  {
    //Duplicate files
    milx::PrintInfo("Duplicating surfaces: ");
    operation = duplicate;
  }
  if(convertArg.isSet())
  {
    //Info
    milx::PrintInfo("Converting surfaces: ");
    operation = convert;
  }
  if(concatenateArg.isSet() || colourConcatenateArg.isSet())
  {
    //Cat
    milx::PrintInfo("Concatenating surfaces: ");
    operation = cat;
  }
  if(mseArg.isSet())
  {
    //Cat
    milx::PrintInfo("MSE of surfaces: ");
    operation = mse;
  }
  if(componentsArg.isSet())
  {
    if(!splitArg.isSet())
    {
      milx::PrintError("Components Argument Error: One or more of the operations");
      milx::PrintError("required to use components is not set (such as split)");
      milx::PrintError("Operations needing components are split.");
      exit(EXIT_FAILURE);
    }
    milx::PrintInfo("Using components: ");
    for(stringiterator name = componentNames.begin(); name != componentNames.end(); name ++)
      std::cout << *name << ", ";
    std::cout << std::endl;
  }
  if(splitArg.isSet())
  {
    //Split
    if(!componentsArg.isSet())
    {
      milx::PrintError("Split Argument Error: Components of the given surfaces must be set.");
      milx::PrintError("These surfaces provide info for the splitting process.");
      milx::PrintError("Set components via the component argument as many times as needed.");
      exit(EXIT_FAILURE);
    }
    milx::PrintInfo("For splitting surfaces: ");
    operation = split;
    componentsValid = true;
  }
  if(scaleArg.isSet())
  {
    //Scale
    milx::PrintInfo("Scaling coordinates of each point in the surfaces: ");
    operation = scale;
  }
  if(smoothArg.isSet())
  {
    //Smooth
    milx::PrintInfo("Smoothing points in the surfaces: ");
    operation = smooth;
  }
  if(laplacianArg.isSet())
  {
    //Smooth Laplacian
    milx::PrintInfo("Laplacian smoothing points in the surfaces: ");
    operation = laplacian;
  }
  if(decimateArg.isSet())
  {
    //Smooth
    milx::PrintInfo("Decimating number of points in the surfaces: ");
    operation = decimate;
  }
  if(thresholdAboveArg.isSet() || thresholdBelowArg.isSet())
  {
    if(!prefixArg.isSet())
    {
      milx::PrintError("Threshold Argument Error: Output Prefix (-p) must be provided.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }
    if(thresholdAboveArg.isSet() && !thresholdBelowArg.isSet())
    {
      thresholdBelow = -std::numeric_limits<float>::max(); //min negative value
    }
    else if(!thresholdAboveArg.isSet() && thresholdBelowArg.isSet())
    {
      thresholdAbove = std::numeric_limits<float>::max();
    }
    std::cout << "Threshold Above Value: " << thresholdAbove << std::endl;
    std::cout << "Threshold Below Value: " << thresholdBelow << std::endl;
    //threshold
    milx::PrintInfo("Thresholding Scalars of surfaces: ");
    operation = thresholdscalars;
  }
  if(clipArg.isSet())
  {
    //Scale
    milx::PrintInfo("Clipping surfaces: ");
    operation = clip;
  }
  if(diffScalarArg.isSet())
  {
    //Scale
    milx::PrintInfo("Differencing Scalars of surfaces: ");
    operation = diffscalars;
  }
  if(flipArg.isSet())
  {
    //Scale
    milx::PrintInfo("Flipping surfaces: ");
    operation = flip;
  }
  if(statsScalarArg.isSet())
  {
    //Scale
    milx::PrintInfo("Statistics of Scalars in surfaces: ");
    operation = statscalars;
  }
  if(removeScalarArg.isSet())
  {
    //Remove scalars
    milx::PrintInfo("Removing Scalars of surfaces: ");
    operation = removescalars;
  }
  if(copyScalarArg.isSet())
  {
    if(filenames.size() == 1)
    {
      milx::PrintError("Argument Error: Need another filename to copy from mesh A to mesh B.");
      exit(EXIT_FAILURE);
    }
    else if(filenames.size() == 2 && !outputArg.isSet())
    {
      milx::PrintError("Argument Error: Use Output (-o) for copying from mesh A to mesh B.");
      milx::PrintError("Re-run with the output name set.");
      exit(EXIT_FAILURE);
    }
    else if(filenames.size() > 2 && !prefixArg.isSet())
    {
      milx::PrintError("Argument Error: Output Prefix (-p) must be provided for more than 2 meshes.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }

    //Remove scalars
    milx::PrintInfo("Copying Scalars of surfaces: ");
    operation = copyscalars;
  }
  if(rigidArg.isSet() && (!procrustesArg.isSet() && !icpArg.isSet()))
  {
    //Scale
    milx::PrintError("Rigid option can only be used with the Procrustes/ICP option");
    exit(EXIT_FAILURE);
  }
  if(procrustesArg.isSet())
  {
    if(!prefixArg.isSet())
    {
      milx::PrintError("Procrustes Argument Error: Output Prefix (-p) must be provided.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }
    //align
    milx::PrintInfo("Computing Procrustes Alignment of surfaces: ");
    operation = procrustes;
  }
  if(icpArg.isSet())
  {
    if(!prefixArg.isSet())
    {
      milx::PrintError("ICP Argument Error: Output Prefix (-p) must be provided.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }
    //align
    milx::PrintInfo("Computing ICP Alignment of surfaces: ");
    operation = icp;
  }

  ///Partition list into two if option given
  std::vector<std::string> altfilenames;
  if(partitionArg.isSet())
  {
    for(size_t j = filenames.size()/2; j < filenames.size(); j ++)
      altfilenames.push_back( filenames[j] );

    for(size_t j = 0; j < altfilenames.size(); j ++)
      filenames.pop_back();
  }

  std::cout << "Total Surfaces: " << filenames.size() << std::endl;
  for(stringiterator name = filenames.begin(); name != filenames.end(); name ++)
    std::cout << *name << ", ";
  std::cout << std::endl;

  if(partitionArg.isSet())
  {
    if(!diffScalarArg.isSet())
    {
      milx::PrintError("Partitioning list of surfaces is only supported with the following operations:");
      milx::PrintError(diffScalarArg.getName());
      exit(EXIT_FAILURE);
    }
    if(!prefixArg.isSet())
    {
      milx::PrintError("Partition Surfaces List Argument Error: Output Prefix (-p) must be provided.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }

    //alt
    milx::PrintInfo("Operation will be applied in conjuction with surfaces: ");
    std::cout << "Total Surfaces: " << altfilenames.size() << std::endl;
    for(stringiterator name = altfilenames.begin(); name != altfilenames.end(); name ++)
      std::cout << *name << ", ";
    std::cout << std::endl;

    if(diffScalarArg.isSet())
      operation = diffscalarspairs;
  }

  //---------------------------
  ///Open files
  std::cerr << "Reading Surfaces... ";
  vtkSmartPointer<vtkPolyDataCollection> collection;
  if( !milx::File::OpenModelCollection(filenames, collection) ) //Error printed inside
    exit(EXIT_FAILURE);
  std::cerr << "Done" << std::endl;
  std::cout << "Read " << collection->GetNumberOfItems() << " models" << std::endl;

  vtkSmartPointer<vtkPolyDataCollection> altcollection;
  if(partitionArg.isSet())
  {
    std::cerr << "Reading Alterate Surfaces... ";
    if( !milx::File::OpenModelCollection(altfilenames, altcollection) ) //Error printed inside
      exit(EXIT_FAILURE);
    std::cerr << "Done" << std::endl;
  }

  vtkSmartPointer<vtkPolyDataCollection> components;
  if(componentsArg.isSet() && componentsValid)
  {
    std::cerr << "Reading Components... ";
    if( !milx::File::OpenModelCollection(componentNames, components) ) //Error printed inside
      exit(EXIT_FAILURE);
    std::cerr << "Done" << std::endl;
  }

  if(colourConcatenateArg.isSet())
  {
    size_t count = 0;
    collection->InitTraversal();
    for(int j = 0; j < collection->GetNumberOfItems(); j ++)
    {
      vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
      const size_t numPoints = mesh->GetNumberOfPoints();

      vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
        weights->SetNumberOfComponents(1);
        weights->SetNumberOfTuples(numPoints);
        weights->FillComponent(0, count);

      mesh->GetPointData()->SetScalars(weights);
      count ++;
    }
  }

  //---------------------------
  ///Geometric transforms potentially associated with the 'collection' item (e.g. with --icp).
  vtkSmartPointer<vtkTransformCollection> transformCollection = vtkSmartPointer<vtkTransformCollection>::New();

  //---------------------------
  ///Operate
  milx::Model Model;
  bool outputRequired = true;
  bool multiOutputRequired = false;
  bool standardOperation = false;
  vtkSmartPointer<vtkPolyDataCollection> resultCollections;

  if (collection->GetNumberOfItems() == 1)
  {
    outputRequired = true;
    multiOutputRequired = false;
  }
  else
  {
    outputRequired = false;
    multiOutputRequired = true;
  }

  std::cerr << "Applying... ";
  switch(operation)
  {
    case duplicate: //Tested 04/2016
      standardOperation = true;
      break;

    case convert: //Tested 04/2016
      standardOperation = true; //requires user to set output format
      break;

    case cat: //Tested 04/2016
      Model.ConcatenateCollection(collection);
      outputRequired = true;
      multiOutputRequired = false;
      break;

    case mse: //--------------------------------
      milx::PrintInfo( "MSE: " + milx::NumberToString(Model.MeanSquaredErrorCollection(collection)) );
      outputRequired = false; //write mean
      break;

    case split: //--------------------------------
      {// scope for local variables
      std::vector< vtkSmartPointer<vtkPolyDataCollection> > splitCollections;

      //Split surfaces
      Model.SplitCollection(collection, components, splitCollections);

      //write
      typedef std::vector< vtkSmartPointer<vtkPolyDataCollection> >::iterator iterator;
      stringiterator name = filenames.begin();
      for(iterator surfaces = splitCollections.begin();
        surfaces != splitCollections.end();
        surfaces ++, name ++)
      {
        std::vector<std::string> splitNames;

        std::string splitName = prefixName + milx::File::GetBaseName(*name) + "_";
        std::string ext = milx::File::GetFileExtension(*name);
        for(stringiterator componentname = componentNames.begin(); componentname != componentNames.end(); componentname ++)
        {
          splitNames.push_back(splitName + milx::File::GetBaseName(*componentname) + "." + ext);
          milx::PrintInfo("Wrote " + splitName + milx::File::GetBaseName(*componentname) + "." + ext);
        }
        if( !milx::File::SaveModelCollection(splitNames, *surfaces) )
          exit(EXIT_FAILURE);
      }
      }// scope for local variables

      outputRequired = false;
      multiOutputRequired = false;
      break;

    case scale: //Tested 04/2016
      Model.ScaleCollection(collection, scaleFactor);

      standardOperation = true;
      break;

    case smooth: //Tested 04/2016
      Model.SmoothCollection(collection, smoothIterations);

      standardOperation = true;
      break;

    case laplacian: //Tested 04/2016
      Model.LaplacianCollection(collection, laplacianIterations);

      standardOperation = true;
      break;

    case decimate: //Tested 04/2016
      Model.DecimateCollection(collection, decimateFactor);

      standardOperation = true;
      break;

    case thresholdscalars:
      Model.ScalarThresholdCollection(collection, thresholdAbove, thresholdBelow);

      standardOperation = true;
      break;

    case clip:
      Model.ClipCollection(collection, clipValue, clipValue);

      standardOperation = true;
      break;

    case flip:
      {
        bool xAxis = false, yAxis = false, zAxis = false;
        if(flipAxis == 2)
          zAxis = true;
        else if(flipAxis == 1)
          yAxis = true;
        else
          xAxis = true;
        Model.FlipCollection(collection, xAxis, yAxis, zAxis);

        standardOperation = true;
        break;
      }

    case diffscalars:
      Model.ScalarDifferenceCollection(collection);
      break;

    case diffscalarspairs:
      resultCollections = vtkSmartPointer<vtkPolyDataCollection>::New();
      Model.ScalarDifferenceCollection(collection, altcollection, resultCollections);

      standardOperation = true;
      collection = resultCollections;
      break;

    case statscalars:
      Model.ScalarStatisticsCollection(collection);
      break;

    case removescalars:
      Model.ScalarRemoveCollection(collection);
    
      standardOperation = true;
      break;

    case copyscalars:
      Model.ScalarCopyCollection(collection);

      if(filenames.size() == 2)
      {
        outputRequired = true;
        multiOutputRequired = false;
        collection->InitTraversal();
        collection->GetNextItem();
        Model.SetInput(collection->GetNextItem());
      }
      else
        standardOperation = true;
      break;

    case procrustes:
    {
      collection = Model.ProcrustesAlignCollection(collection, rigidArg.isSet()); ///\todo leak here?

      std::string ext = milx::File::GetFileExtension(filenames[0]);
      outputName = prefixName + "_mean." + ext;
      outputRequired = true; //write mean
      multiOutputRequired = true; //write collection
      break;
    } //scope for ext

    case icp:
      collection = Model.IterativeClosestPointsAlignCollection(collection, rigidArg.isSet(), transformCollection);

      standardOperation = true; //write collection
      break;

    case none: //--------------------------------
      break;
  }
  std::cerr << "Done" << std::endl;

  if (collection->GetNumberOfItems() == 1 && standardOperation)
  {
    collection->InitTraversal();
    Model.Result() = collection->GetNextItem();
  }

  //---------------------------
  ///Write result
  if(outputRequired)
  {
    milx::PrintInfo("Writing result to " + outputName);
    milx::File::SaveModel(outputName, Model.Result());
  }
  if(multiOutputRequired)
  {
    if(collection->GetNumberOfItems() > 0)
    {
      const int N = collection->GetNumberOfItems();

      milx::PrintInfo("Writing results with prefix " + prefixName);
      std::vector<std::string> names;

      int n = 0;
      for(stringiterator filename = filenames.begin(); filename != filenames.end() && n < N; filename ++, n ++)
      {
        std::string ext;
        if (outputFormatArg.isSet()) {
          ext = outputFormatArg.getValue();
        } else {
          ext = milx::File::GetFileExtension(*filename);
        }

        names.push_back(prefixName + milx::File::GetBaseName(*filename) + "." + ext);
        milx::PrintInfo("Will be writing " + prefixName + milx::File::GetBaseName(*filename) + "." + ext);
      }

      milx::File::SaveModelCollection(names, collection);

      // Saving the vtkTransforms
      if (saveTransformsArg.isSet() && transformCollection->GetNumberOfItems() == collection->GetNumberOfItems())
      {
        const std::string ext = "trsf.gz";
        std::vector<std::string> tnames;
        n = 0;
        for(stringiterator name = names.begin(); name != names.end() && n < N; name ++, n++)
        {
          tnames.push_back(milx::File::StripFileExtension(*name) + "." + ext);
          milx::PrintInfo("Will be writing " + tnames.back());
        }
        milx::File::SaveTransformCollection(tnames, transformCollection, false);
      }
    }
    else
    {
      milx::PrintError("Result was empty. Skipping Output.");
    }
  }

  milx::PrintInfo("Operation Complete");
  return EXIT_SUCCESS;
}
