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

// SMILI
#include <milxFile.h>
#include <milxImage.h>

using namespace TCLAP;

typedef std::vector< std::string >::iterator stringiterator;

//Supported operations
enum operations {none = 0, info, convert, labelinfo, rescale, invert, relabel, smooth, median, gradmag, laplacian, distancemap, threshold, Otsu, crop, mask, resample, match, checker, add, diff, mean, merge, cast, flip};

//Image stuff
typedef unsigned char charPixelType;
typedef itk::Image<charPixelType, milx::imgDimension> charImageType;
typedef float floatPixelType;
typedef itk::Image<floatPixelType, milx::imgDimension> floatImageType;
typedef itk::VectorImage<floatPixelType, milx::imgDimension> vectorImageType;

/**
  \file milxImageApp.cxx
  \brief A swiss-army knife of image processing. This application implements single/batch processing of images with various processing options found in the SMILI library.
  \ingroup Applications

  You are able to smooth, mask, threshold etc. a number of images easily via the commandline with this application (use --help). This application uses the milxImage class to do all processing.

  In the following examples, to process all nifti images in MRIs directory (you can just provide a single image filename to do just one):

  Otsu multiple threshold:
  \verbatim
  milxImageApp --Otsu 128 features/median/*.nii.gz --labels 4 -p features/Otsu/Otsu_128bins_
  \endverbatim

  Rescale image intensities:
  \verbatim
  milxImageApp --rescale MRIs/*.nii.gz --above 255 --below 0 -p features/rescale/rescale_
  \endverbatim

  Labelling info:
  \verbatim
  milxImageApp --info --labels 0 Labels/*.nii.gz
  \endverbatim
  The unique label values present in each image is output. A concatenated output is provided at the end for analysis (histograms etc.)

  Checkerboards Atlas_MRI_Mean2_R_Preprocessed_x2.nii.gz to all Nifti images in current directory and writes to files with prefix check_:
  \verbatim
  milxImageApp --checkerboard SyngoHipData/Atlas_MRI_Mean2_R_Preprocessed_x2.nii.gz *.nii.gz -p check_
  \endverbatim
*/
int main(int argc, char *argv[])
{
  //---------------------------
  ///Program Info
  milx::PrintInfo("--------------------------------------------------------");
  milx::PrintInfo("MILX-SMILI Image Diagnostic Tool for Images/Volumes.");
  milx::PrintInfo("(c) Copyright Shekhar Chandra et al., 2013.");
  milx::PrintInfo("Version: " + milx::NumberToString(milx::Version));
  milx::PrintInfo("Australian e-Health Research Centre, CSIRO.");
  milx::PrintInfo("--------------------------------------------------------\n");

  //---------------------------
  ///Process Arguments
  CmdLine cmd("A diagnostic tool for image/volume operations", ' ', milx::NumberToString(milx::Version));

  ///Optional
  ValueArg<std::string> outputArg("o", "output", "Output Image", false, "result.nii.gz", "Output");
  ValueArg<std::string> prefixArg("p", "prefix", "Output prefix for multiple output.", false, "img_", "Output Prefix");
  //~ ValueArg<float> decimateArg("d", "decimate", "Decimate all the meshes provided using the Quadric Decimate algorithm", false, 0.5, "Decimate");
  ValueArg<float> smoothArg("s", "smooth", "Smooth the images using the Gradient Anisotropic Diffusion algorithm given timestep (use negative value to auto select).", false, -1.0, "Smooth");
  ValueArg<size_t> medianArg("", "median", "Smooth the images using the median filtering given the number of pixels in neighbourhood.", false, 1, "Median");
  ValueArg<size_t> mergeArg("", "merge", "Merge the labels in labelled images, 0:keep, 1:aggregate, 2:pack, 3:strict.", false, 0, "Merge");
  ValueArg<std::string> cropArg("", "crop", "Masks and crops the images using the mask image name provided.", false, "mask.nii.gz", "Crop");
  ValueArg<std::string> maskArg("m", "mask", "Masks the images using the mask image name provided.", false, "mask.nii.gz", "Mask");
  ValueArg<std::string> resampleArg("", "resample", "Resample the images to the reference image provided assuming they are both already in the correct space.", false, "reference.nii.gz", "Resample");
  ValueArg<std::string> matchArg("", "match", "Match the histograms of the images to the reference image provided.", false, "reference.nii.gz", "Match");
  ValueArg<std::string> checkerArg("", "checkerboard", "Checkerboard all of the images to the reference image provided.", false, "reference.nii.gz", "Checkerboard");
  //~ ValueArg<float> scaleArg("", "scale", "Scale the coordinates of the point", false, 0.9, "Scale");
  ValueArg<float> aboveArg("", "above", "Add above value to operation (such as to thresholding).", false, 0.0, "Above");
  ValueArg<float> belowArg("", "below", "Add below value to operation (such as to thresholding).", false, 255.0, "Below");
  ValueArg<float> insideArg("", "inside", "Add inside value to operation (such as to thresholding).", false, 1.0, "Inside");
  ValueArg<size_t> iterationsArg("i", "iterations", "Number of iterations to use in the operation (such as smoothing).", false, 5, "Labels");
  ValueArg<size_t> labelsArg("", "labels", "Number of labels to use in the operation (such as Otsu thresholding) else print labelling info.", false, 3, "Labels");
  ValueArg<size_t> OtsuArg("", "Otsu", "Otsu multiple threshold with the number of bins to use.", false, 128, "Otsu");
  ValueArg<size_t> paddingArg("", "padding", "Number of pixels to pad in the operation in question (such as crop).", false, 1, "Padding");
  ValueArg<size_t> flipArg("f", "flip", "Flip the image about origin at axis indicated (0: x-axis, 1:y-axis, 2:z-axis).", false, 0, "Flip");
  ///Switches
  ///XOR Switches
  SwitchArg infoArg("", "info", "Report the image information(s).", false);
  SwitchArg labelInfoArg("", "labelinfo", "Report the label information for a labelled image. Same as --info --labels arguments.", false);
  SwitchArg convertArg("c", "convert", "Convert the image from current format to the one given at output.", false);
  SwitchArg gradMagArg("g", "gradmag", "Compute the Gradient magnitude (show edges) of the image(s).", false);
  SwitchArg laplacianArg("l", "laplacian", "Compute the Laplacian of the image(s).", false);
  SwitchArg distancemapArg("d", "distancemap", "Compute the Signed Distance map of the image(s).", false);
  SwitchArg thresholdArg("t", "threshold", "Compute the Threshold of the image(s).", false);
  SwitchArg binaryArg("b", "binary", "Compute the binary version of the operation(s).", false);
  SwitchArg rescaleArg("r", "rescale", "Re-scale the intensities of the image to a given range (set inconjuction with the above and below arguments).", false);
  SwitchArg invertArg("", "invert", "Invert the intensities of the images.", false);
  SwitchArg relabelArg("", "relabel", "Relabel the labels of the labelled image consecutatively based on connectivity.", false);
  SwitchArg addArg("", "add", "Add/Sum the provided images together (pixel-wise).", false);
  SwitchArg diffArg("", "diff", "Difference the provided images from the first image (pixel-wise).", false);
  SwitchArg meanArg("", "mean", "Average the provided images together (pixel-wise).", false);
  SwitchArg castArg("", "cast", "Cast the images from 8-bit to floating-point type (or vice-versa) depending on the input type.", false);
  //~ SwitchArg mseArg("", "mse", "Mean Squared Error of Points in models", false);
  //~ SwitchArg rigidArg("", "rigid", "Rigid alignment of surfaces assuming points have correspondence.", false);
  //~ SwitchArg saveTransformsArg("", "savetransforms", "Save the transformation matrix after surface alignment. For use with --icp", false);

  ///Mandatory
  UnlabeledMultiArg<std::string> multinames("images", "Images to operate on (pixel type is auto detected from the first image)", true, "Images");

  ///Add argumnets
  cmd.add( multinames );
  cmd.add( outputArg );
  cmd.add( prefixArg );
  //~ cmd.add( rigidArg );
  //~ cmd.add( saveTransformsArg );
  cmd.add( aboveArg );
  cmd.add( belowArg );
  cmd.add( insideArg );
  cmd.add( binaryArg );
  cmd.add( iterationsArg );
  cmd.add( labelsArg );
  cmd.add( paddingArg );
  ///XOR args
  std::vector<Arg*> xorlist;
  xorlist.push_back(&infoArg);
  xorlist.push_back(&convertArg);
  //~ xorlist.push_back(&scaleArg);
  //~ xorlist.push_back(&decimateArg);
  xorlist.push_back(&smoothArg);
  xorlist.push_back(&medianArg);
  xorlist.push_back(&gradMagArg);
  xorlist.push_back(&laplacianArg);
  xorlist.push_back(&distancemapArg);
  xorlist.push_back(&thresholdArg);
  xorlist.push_back(&OtsuArg);
  xorlist.push_back(&rescaleArg);
  xorlist.push_back(&invertArg);
  xorlist.push_back(&relabelArg);
  xorlist.push_back(&maskArg);
  xorlist.push_back(&resampleArg);
  xorlist.push_back(&matchArg);
  xorlist.push_back(&checkerArg);
  xorlist.push_back(&addArg);
  xorlist.push_back(&diffArg);
  xorlist.push_back(&meanArg);
  xorlist.push_back(&castArg);
  xorlist.push_back(&flipArg);
  //~ xorlist.push_back(&mseArg);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
  xorlist.push_back(&labelInfoArg);
  xorlist.push_back(&mergeArg);
  xorlist.push_back(&cropArg);
#endif // ITK_REVIEW
  cmd.xorAdd(xorlist);

  ///Parse the argv array.
  cmd.parse( argc, argv );

  ///Get the value parsed by each arg.
  //Filenames of surfaces
  std::vector<std::string> filenames = multinames.getValue();
  std::string outputName = outputArg.getValue();
  const std::string prefixName = prefixArg.getValue();
  const std::string resampleName = resampleArg.getValue();
  const std::string matchName = matchArg.getValue();
  const std::string checkerName = checkerArg.getValue();
  //~ const float decimateFactor = decimateArg.getValue();
  //~ const float scaleFactor = scaleArg.getValue();
  const float smoothTimestep = smoothArg.getValue();
  const size_t radius = medianArg.getValue();
  const size_t mergeType = mergeArg.getValue();
  float aboveValue = aboveArg.getValue();
  float belowValue = belowArg.getValue();
  float insideValue = insideArg.getValue();
  size_t iterations = iterationsArg.getValue(); //number of labels
  size_t labelsValue = labelsArg.getValue(); //number of labels
  size_t OtsuValue = OtsuArg.getValue(); //number of bins
  size_t paddingValue = paddingArg.getValue(); //number of bins
  size_t flipAxis = flipArg.getValue(); //number of bins

  std::string maskName;
  if(cropArg.isSet())
    maskName = cropArg.getValue();
  else
    maskName = maskArg.getValue();

  ///Display operation
  operations operation = none;
  bool outputRequired = false;
  bool multiOutputRequired = true;
  if(infoArg.isSet())
  {
    outputRequired = false;
    multiOutputRequired = false;

    //Info
    operation = info;
  }
  if(labelInfoArg.isSet())
  {
    outputRequired = false;
    multiOutputRequired = false;

    //Info
    operation = labelinfo;
  }
  if(convertArg.isSet())
  {
    if(!outputArg.isSet())
    {
      milx::PrintError("Argument Error: Use Output (-o) for conversion.");
      milx::PrintError("Re-run with the output name set.");
      exit(EXIT_FAILURE);
    }
    if(filenames.size() != 1)
    {
      milx::PrintError("Argument Error: Only one output is supported for conversion.");
      milx::PrintError("Re-run with the correct number of inputs.");
      exit(EXIT_FAILURE);
    }

    outputRequired = true;
    multiOutputRequired = false;

    //Info
    operation = convert;
  }
  if( smoothArg.isSet() || medianArg.isSet() || gradMagArg.isSet() || laplacianArg.isSet() || distancemapArg.isSet() || cropArg.isSet() || maskArg.isSet()
     || resampleArg.isSet() || matchArg.isSet() || checkerArg.isSet() || thresholdArg.isSet() || OtsuArg.isSet() || rescaleArg.isSet() || invertArg.isSet() || relabelArg.isSet() || addArg.isSet()
     || diffArg.isSet() || meanArg.isSet() || mergeArg.isSet() || castArg.isSet() || flipArg.isSet() )
  {
    ///Check if output argument given and only doing one image
    if(filenames.size() == 1 && (addArg.isSet() || diffArg.isSet() || meanArg.isSet() || mergeArg.isSet()))
    {
      milx::PrintError("Argument Error: Cannot use arithmetic operations on single input.");
      milx::PrintError("Re-run with more inputs set.");
      exit(EXIT_FAILURE);
    }
    else if(filenames.size() == 1 && prefixArg.isSet())
    {
      milx::PrintError("Argument Error: Cannot use prefix argument on single input.");
      milx::PrintError("Re-run with output (-o) option set instead.");
      exit(EXIT_FAILURE);
    }
    else if(filenames.size() == 1)
    {
      outputRequired = true;
      multiOutputRequired = false;
    }
    else if(filenames.size() > 1 && (addArg.isSet() || diffArg.isSet() || meanArg.isSet() || mergeArg.isSet()))
    {
      outputRequired = true;
      multiOutputRequired = false;
    }
    else if(outputArg.isSet() && filenames.size() > 1)
    {
      milx::PrintError("Argument Error: Use Output Prefix (-p) for multiple inputs.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }
    else if(!prefixArg.isSet() && filenames.size() > 1)
    {
      milx::PrintError("Argument Error: Output Prefix (-p) must be provided for multiple inputs.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }

    //Operations allowed
    if(smoothArg.isSet())
      operation = smooth;
    if(medianArg.isSet())
      operation = median;
    if(gradMagArg.isSet())
      operation = gradmag;
    if(laplacianArg.isSet())
      operation = laplacian;
    if(distancemapArg.isSet())
      operation = distancemap;
    if(thresholdArg.isSet())
      operation = threshold;
    if(cropArg.isSet())
      operation = crop;
    if(maskArg.isSet())
      operation = mask;
    if(resampleArg.isSet())
      operation = resample;
    if(matchArg.isSet())
      operation = match;
    if(checkerArg.isSet())
      operation = checker;
    if(OtsuArg.isSet())
      operation = Otsu;
    if(rescaleArg.isSet())
      operation = rescale;
    if(invertArg.isSet())
      operation = invert;
    if(relabelArg.isSet())
      operation = relabel;
    if(addArg.isSet())
      operation = add;
    if(diffArg.isSet())
      operation = diff;
    if(meanArg.isSet())
      operation = mean;
    if(mergeArg.isSet())
      operation = merge;
    if(castArg.isSet())
      operation = cast;
    if(flipArg.isSet())
      operation = flip;
  }
  if(aboveArg.isSet() || belowArg.isSet())
  {
    if(!thresholdArg.isSet() && !rescaleArg.isSet())
    {
      milx::PrintError("Argument Error: Another argument (such as threshold) must be provided.");
      milx::PrintError("Re-run with one of these flags set.");
      exit(EXIT_FAILURE);
    }
  }
  if(labelsArg.isSet())
  {
    if(!OtsuArg.isSet() && infoArg.isSet())
    {
    #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
      milx::PrintError("Argument Warning: Another argument (such as Otsu) not provided. Printing label info instead.");
      operation = labelinfo;
    #else
      milx::PrintError("Argument Warning: Another argument (such as Otsu) not provided and ITK version not sufficient for printing label info. Exiting.");
      exit(EXIT_FAILURE);
    #endif
    }
  }
  if(iterationsArg.isSet())
  {
    if(!smoothArg.isSet())
    {
      milx::PrintError("Argument Error: Another argument (such as smoothing) must be provided.");
      milx::PrintError("Re-run with one of these flags set.");
      exit(EXIT_FAILURE);
    }
  }
  if(paddingArg.isSet())
  {
    if(!cropArg.isSet())
    {
      milx::PrintError("Argument Error: Another argument (such as crop) must be provided.");
      milx::PrintError("Re-run with one of these flags set.");
      exit(EXIT_FAILURE);
    }
  }
  if(flipArg.isSet())
  {
    if(flipAxis < 0 || flipAxis > 2)
    {
      milx::PrintError("Argument Error: Incorrect value provided for axis flipping.");
      milx::PrintError("Re-run with value correctly set.");
      exit(EXIT_FAILURE);
    }
  }
  if(thresholdArg.isSet() || rescaleArg.isSet())
  {
    if( (!aboveArg.isSet() || !belowArg.isSet()) && rescaleArg.isSet() )
    {
      milx::PrintError("Argument Error: Above and below argument(s) must be provided.");
      milx::PrintError("Re-run with one of these flags set.");
      exit(EXIT_FAILURE);
    }
    else
    {
      if( (!aboveArg.isSet() && !belowArg.isSet()) && !binaryArg.isSet() )
      {
        milx::PrintError("Argument Error: Either above or below argument(s) must be provided.");
        milx::PrintError("Re-run with one of these flags set.");
        exit(EXIT_FAILURE);
      }
      if( (!aboveArg.isSet() || !belowArg.isSet() || !insideArg.isSet()) && binaryArg.isSet() )
      {
        milx::PrintError("Argument Error: The inside, above and below argument(s) must be provided.");
        milx::PrintError("Re-run with one of these flags set.");
        exit(EXIT_FAILURE);
      }
    }
  }

  std::cout << "Total Images: " << filenames.size() << std::endl;
  if(filenames.empty())
  {
    milx::PrintError("No image file names provided. Exiting.");
    exit(EXIT_FAILURE);
  }
  for(stringiterator name = filenames.begin(); name != filenames.end(); name ++)
    std::cout << *name << ", ";
  std::cout << std::endl;

  //---------------------------
  ///Open files
  std::cerr << "Reading Header for type info etc." << std::endl;
  std::string pixelType, componentType;
  size_t dimensions = 3;
  if( !milx::File::ReadImageInformation(filenames[0], pixelType, componentType, dimensions) )
  {
    milx::PrintError("Failed Reading First Image. Check the image type/file. Exiting.");
    exit(EXIT_FAILURE);
  }
  if(infoArg.isSet() || labelInfoArg.isSet())
  {
    milx::PrintInfo("Pixel Type: " + pixelType);
    milx::PrintInfo("Component Type: " + componentType);
    milx::PrintInfo("Dimensions: " + milx::NumberToString(dimensions));
  }

  std::cerr << "Reading Images... ";
  std::vector< itk::SmartPointer<vectorImageType> > vectorCollection;
  std::vector< itk::SmartPointer<charImageType> > labelledCollection;
  std::vector< itk::SmartPointer<floatImageType> > collection;
  bool labelledImages = false, vectorImages = false;
  if(pixelType == "vector" || dimensions > 3) ///\todo handle 4D images here properly
  {
    milx::PrintInfo("Detected vector images.");
    if( !milx::File::OpenImages<vectorImageType>(filenames, vectorCollection) ) //Error NOT printed inside
    {
      milx::PrintError("Failed Reading Labelled Images. Exiting.");
      exit(EXIT_FAILURE);
    }
    vectorImages = true;
    std::cout << "Read " << vectorCollection.size() << " vector images" << std::endl;

  #if ITK_MAJOR_VERSION > 3
    if(!infoArg.isSet() && !maskArg.isSet() && !addArg.isSet() && !diffArg.isSet() && !meanArg.isSet())
  #else
    if(!infoArg.isSet() && !addArg.isSet() && !diffArg.isSet() && !meanArg.isSet())
  #endif
    {
      milx::PrintError("Input Error: Operation provided not supported for vector images yet.");
      exit(EXIT_FAILURE);
    }
  }
  else if(componentType == "unsigned_char" || componentType == "unsigned char")
  {
    milx::PrintInfo("Detected labelled images.");
    if( !milx::File::OpenImages<charImageType>(filenames, labelledCollection) ) //Error NOT printed inside
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
    if( !milx::File::OpenImages<floatImageType>(filenames, collection) ) //Error NOT printed inside
    {
      milx::PrintError("Failed Reading Images. Exiting.");
      exit(EXIT_FAILURE);
    }
    std::cout << "Read " << collection.size() << " images" << std::endl;
  }

  charImageType::Pointer maskImage;
  if(cropArg.isSet() || maskArg.isSet())
  {
    if(!milx::File::OpenImage<charImageType>(maskName, maskImage))
    {
      milx::PrintError("Failed Reading Mask. Exiting.");
      exit(EXIT_FAILURE);
    }
    std::cout << "Read " << maskName << " as mask" << std::endl;
  }
  std::cerr << "Done" << std::endl;

  //---------------------------
  ///Operate
  vectorImageType::Pointer vecResult;
  charImageType::Pointer charResult;
  floatImageType::Pointer floatResult;
  std::cerr << "Applying... ";
  if(vectorImages)
  {
    switch(operation)
    {
      case info:
        milx::Image<vectorImageType>::InformationCollection(vectorCollection);
        break;

      case convert:

        break;

    #if ITK_MAJOR_VERSION > 3
      case mask:
        milx::Image<vectorImageType>::MaskCollection<charImageType>(vectorCollection, maskImage);
        break;
    #endif

      case add:
        vecResult = milx::Image<vectorImageType>::AddCollection(vectorCollection);
        break;

      case diff:
        vecResult = milx::Image<vectorImageType>::DifferenceCollection(vectorCollection);
        break;

      case mean:
        vecResult = milx::Image<vectorImageType>::AverageVectorCollection(vectorCollection, vectorCollection[0]->GetNumberOfComponentsPerPixel());
        break;

      case none: //--------------------------------
        break;

      default:
        milx::PrintError("Operation not supported. Exiting");
        exit(EXIT_FAILURE);
        break;
    }
  }
  else if(labelledImages)
  {
    switch(operation)
    {
      case info:
        milx::Image<charImageType>::InformationCollection(labelledCollection);
        break;

      case convert:

        break;

    #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
      case labelinfo:
      {
        milx::PrintInfo("\n>> Labels present: ");
        milx::PrintInfo("-----------------------------------------------");
        std::vector<unsigned char> allValues;
        for(size_t k = 0; k < labelledCollection.size(); k ++)
        {
            milx::PrintInfo("Image " + milx::NumberToString(k));
            milx::PrintInfo(milx::File::ExtractFilename(filenames[k]) + ": \n");
            std::vector<unsigned char> values = milx::Image<charImageType>::LabelValues(labelledCollection[k]);
            for(size_t j = 0; j < values.size(); j ++)
            {
                allValues.push_back(values[j]);
                std::cout << static_cast<unsigned>(values[j]) << ", ";
            }
            std::cout << std::endl;
        }
        milx::PrintInfo("\n>> All Labels present (as single list, n = " + milx::NumberToString(allValues.size()) + "): ");
        for(size_t j = 0; j < allValues.size(); j ++)
            std::cout << static_cast<unsigned>(allValues[j]) << ", ";
        std::cout << std::endl;
        break;
      }
    #endif

      case rescale:
        milx::Image<charImageType>::RescaleIntensityCollection(labelledCollection, belowValue, aboveValue);
        break;

      case invert:
        milx::Image<charImageType>::InvertIntensityCollection(labelledCollection);
        break;

      case relabel:
        milx::Image<charImageType>::RelabelCollection(labelledCollection);
        break;

      case smooth:
        {
          collection = milx::Image<charImageType>::AnisotropicDiffusionCollection<floatImageType>(labelledCollection, iterations, smoothTimestep);
          labelledImages = false;
          break;
        }

      case median:
        milx::Image<charImageType>::MedianCollection(labelledCollection, radius);
        break;

      case gradmag:
        milx::Image<charImageType>::GradientMagnitudeCollection(labelledCollection);
        break;

      case laplacian:
        milx::Image<charImageType>::LaplacianCollection(labelledCollection);
        break;

      case distancemap:
        {
          const bool binary = false, signedDistance = true, insideDistance = false, squaredDistance = false;

          collection = milx::Image<charImageType>::DistanceMapCollection<floatImageType>(labelledCollection, binary, signedDistance, insideDistance, squaredDistance);

          labelledImages = false;
          break;
        }

      case threshold:
        if(aboveArg.isSet() && !belowArg.isSet())
        {
          milx::PrintInfo("Thresholding Above " + milx::NumberToString(aboveValue));
          milx::Image<charImageType>::ThresholdAboveCollection(labelledCollection, 0.0, aboveValue);
        }
        else if(!aboveArg.isSet() && belowArg.isSet())
        {
          milx::PrintInfo("Thresholding Below " + milx::NumberToString(belowValue));
          milx::Image<charImageType>::ThresholdBelowCollection(labelledCollection, 0.0, belowValue);
        }
        else if(binaryArg.isSet())
        {
          milx::PrintInfo("Binary Thresholding Above " + milx::NumberToString(aboveValue) + " and Below " + milx::NumberToString(belowValue));
          labelledCollection = milx::Image<charImageType>::BinaryThresholdCollection<charImageType>(labelledCollection, 0.0, insideValue, belowValue, aboveValue);
        }
        else
        {
          milx::PrintInfo("Thresholding Above " + milx::NumberToString(aboveValue) + " and Below " + milx::NumberToString(belowValue));
          milx::Image<charImageType>::ThresholdCollection(labelledCollection, 0.0, belowValue, aboveValue);
        }

        break;

      case Otsu:
        milx::PrintInfo("Otsu Multiple Thresholding with " + milx::NumberToString(OtsuValue) + " and " + milx::NumberToString(labelsValue) + " labels");
        labelledCollection = milx::Image<charImageType>::OtsuMultipleThresholdCollection<charImageType>(labelledCollection, OtsuValue, labelsValue);
        break;

      case crop:
      #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
        milx::Image<charImageType>::MaskAndCropCollection<charImageType>(labelledCollection, maskImage, paddingValue);
      #endif
        break;

      case mask:
        milx::Image<charImageType>::MaskCollection<charImageType>(labelledCollection, maskImage);
        break;

      case resample:
        {
          charImageType::Pointer referenceImage;
          if(resampleArg.isSet())
          {
            if(!milx::File::OpenImage<charImageType>(resampleName, referenceImage))
            {
              milx::PrintError("Failed Reading Reference Image. Exiting.");
              exit(EXIT_FAILURE);
            }
            std::cout << "Read " << resampleName << " as the reference image for resampling" << std::endl;
          }
          milx::Image<charImageType>::ResampleLabelCollection<charImageType>(labelledCollection, referenceImage);
          break;
        }

      case match:
        {
          charImageType::Pointer referenceImage;
          if(matchArg.isSet())
          {
            if(!milx::File::OpenImage<charImageType>(matchName, referenceImage))
            {
              milx::PrintError("Failed Reading Reference Image. Exiting.");
              exit(EXIT_FAILURE);
            }
            std::cout << "Read " << matchName << " as the reference image for matching" << std::endl;
          }
          milx::Image<charImageType>::MatchHistogramCollection(labelledCollection, referenceImage);
          break;
        }

      case checker:
        {
          charImageType::Pointer referenceImage;
          if(checkerArg.isSet())
          {
            if(!milx::File::OpenImage<charImageType>(checkerName, referenceImage))
            {
              milx::PrintError("Failed Reading Reference Image. Exiting.");
              exit(EXIT_FAILURE);
            }
            std::cout << "Read " << checkerName << " as the reference image for checkerboarding" << std::endl;
          }
          milx::Image<charImageType>::CheckerboardCollection(labelledCollection, referenceImage);
          break;
        }

      case add:
        charResult = milx::Image<charImageType>::AddCollection(labelledCollection);
        break;

      case diff:
        charResult = milx::Image<charImageType>::DifferenceCollection(labelledCollection);
        break;

      case mean:
        floatResult = milx::Image<charImageType>::AverageCollection<floatImageType>(labelledCollection);
        labelledImages = false;
        break;


      case merge:
      #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
        charResult = milx::Image<charImageType>::MergeLabelledImages(labelledCollection, mergeType);
      #endif
        break;

      case cast:
        collection = milx::Image<charImageType>::CastCollection<floatImageType>(labelledCollection);
        labelledCollection.clear();
        labelledImages = false;
        break;

      case flip:
        if(flipAxis == 0)
          milx::Image<charImageType>::FlipCollection(labelledCollection, true, false, false, true);
        else if(flipAxis == 1)
          milx::Image<charImageType>::FlipCollection(labelledCollection, false, true, false, true);
        else
          milx::Image<charImageType>::FlipCollection(labelledCollection, false, false, true, true);
        break;

      case none: //--------------------------------
        break;

      default:
        milx::PrintError("Operation not supported. Exiting");
        exit(EXIT_FAILURE);
        break;
    }
    std::cerr << "Done" << std::endl;
  }
  else
  {
    switch(operation)
    {
      case info:
        milx::Image<floatImageType>::InformationCollection(collection);
        break;

      case convert:

        break;

    #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
      case labelinfo:
      {
        milx::PrintWarning("Float images found and will be cast to 8-bit images for labelling info.");
        labelledCollection = milx::Image<floatImageType>::CastCollection<charImageType>(collection);
        collection.clear();

        milx::PrintInfo("\n>> Labels present: ");
        milx::PrintInfo("-----------------------------------------------");
        std::vector<unsigned char> allValues;
        for(size_t k = 0; k < labelledCollection.size(); k ++)
        {
            milx::PrintInfo("Image " + milx::NumberToString(k));
            milx::PrintInfo(milx::File::ExtractFilename(filenames[k]) + ": \n");
            std::vector<unsigned char> values = milx::Image<charImageType>::LabelValues(labelledCollection[k]);
            for(size_t j = 0; j < values.size(); j ++)
            {
                allValues.push_back(values[j]);
                std::cout << static_cast<unsigned>(values[j]) << ", ";
            }
            std::cout << std::endl;
        }
        milx::PrintInfo("\n>> All Labels present (as single list, n = " + milx::NumberToString(allValues.size()) + "): ");
        for(size_t j = 0; j < allValues.size(); j ++)
            std::cout << static_cast<unsigned>(allValues[j]) << ", ";
        std::cout << std::endl;
        break;
      }
    #endif

      case rescale:
        milx::Image<floatImageType>::RescaleIntensityCollection(collection, belowValue, aboveValue);
        break;

      case invert:
        milx::Image<floatImageType>::InvertIntensityCollection(collection);
        break;

      case relabel:
        milx::PrintError("Failed Relabelling. Only applicable to Labelled images. Exiting.");
        exit(EXIT_FAILURE);

      case smooth:
        collection = milx::Image<floatImageType>::AnisotropicDiffusionCollection<floatImageType>(collection, iterations, smoothTimestep);
        break;

      case median:
        milx::Image<floatImageType>::MedianCollection(collection, radius);
        break;

      case gradmag:
        milx::Image<floatImageType>::GradientMagnitudeCollection(collection);
        break;

      case laplacian:
        milx::Image<floatImageType>::LaplacianCollection(collection);
        break;

      case distancemap:
        {
          const bool binary = false, signedDistance = true, insideDistance = false, squaredDistance = false;

          collection = milx::Image<floatImageType>::DistanceMapCollection<floatImageType>(collection, binary, signedDistance, insideDistance, squaredDistance);

          break;
        }

      case threshold:
        if(aboveArg.isSet() && !belowArg.isSet())
        {
          milx::PrintInfo("Thresholding Above " + milx::NumberToString(aboveValue));
          milx::Image<floatImageType>::ThresholdAboveCollection(collection, 0.0, aboveValue);
        }
        else if(!aboveArg.isSet() && belowArg.isSet())
        {
          milx::PrintInfo("Thresholding Below " + milx::NumberToString(belowValue));
          milx::Image<floatImageType>::ThresholdBelowCollection(collection, 0.0, belowValue);
        }
        else if(binaryArg.isSet())
        {
          milx::PrintInfo("Binary Thresholding Above " + milx::NumberToString(aboveValue) + " and Below " + milx::NumberToString(belowValue));
          labelledCollection = milx::Image<floatImageType>::BinaryThresholdCollection<charImageType>(collection, 0.0, insideValue, belowValue, aboveValue);
          labelledImages = true;
        }
        else
        {
          milx::PrintInfo("Thresholding Above " + milx::NumberToString(aboveValue) + " and Below " + milx::NumberToString(belowValue));
          milx::Image<floatImageType>::ThresholdCollection(collection, 0.0, belowValue, aboveValue);
        }

        break;

      case Otsu:
        milx::PrintInfo("Otsu Multiple Thresholding with " + milx::NumberToString(OtsuValue) + " and " + milx::NumberToString(labelsValue) + " labels");
        labelledCollection = milx::Image<floatImageType>::OtsuMultipleThresholdCollection<charImageType>(collection, OtsuValue, labelsValue);
        labelledImages = true;

        break;

      case crop:
      #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
        milx::Image<floatImageType>::MaskAndCropCollection<charImageType>(collection, maskImage, paddingValue);
      #endif

        break;

      case mask:
        milx::Image<floatImageType>::MaskCollection<charImageType>(collection, maskImage);
        break;

      case resample:
        {
          floatImageType::Pointer referenceImage;
          if(resampleArg.isSet())
          {
            if(!milx::File::OpenImage<floatImageType>(resampleName, referenceImage))
            {
              milx::PrintError("Failed Reading Reference Image. Exiting.");
              exit(EXIT_FAILURE);
            }
            std::cout << "Read " << resampleName << " as the reference image for resampling" << std::endl;
          }
          milx::Image<floatImageType>::ResampleCollection<floatImageType>(collection, referenceImage);
          break;
        }

      case match:
        {
          floatImageType::Pointer referenceImage;
          if(matchArg.isSet())
          {
            if(!milx::File::OpenImage<floatImageType>(matchName, referenceImage))
            {
              milx::PrintError("Failed Reading Reference Image. Exiting.");
              exit(EXIT_FAILURE);
            }
            std::cout << "Read " << matchName << " as the reference image for matching" << std::endl;
          }
          milx::Image<floatImageType>::MatchHistogramCollection(collection, referenceImage);
          break;
        }

      case checker:
        {
          floatImageType::Pointer referenceImage;
          if(checkerArg.isSet())
          {
            if(!milx::File::OpenImage<floatImageType>(checkerName, referenceImage))
            {
              milx::PrintError("Failed Reading Reference Image. Exiting.");
              exit(EXIT_FAILURE);
            }
            std::cout << "Read " << checkerName << " as the reference image for checkerboarding" << std::endl;
          }
          milx::Image<floatImageType>::CheckerboardCollection(collection, referenceImage);
          break;
        }

      case add:
        floatResult = milx::Image<floatImageType>::AddCollection(collection);
        break;

      case diff:
        floatResult = milx::Image<floatImageType>::DifferenceCollection(collection);
        break;

      case mean:
        floatResult = milx::Image<floatImageType>::AverageCollection<floatImageType>(collection);
        break;

      case merge:
      #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
        milx::PrintWarning("Float images found and will be cast to 8-bit images for merge.");
        labelledCollection = milx::Image<floatImageType>::CastCollection<charImageType>(collection);
        collection.clear();
        charResult = milx::Image<charImageType>::MergeLabelledImages(labelledCollection, mergeType);
        labelledImages = true;
      #endif
        break;

      case cast:
        labelledCollection = milx::Image<floatImageType>::CastCollection<charImageType>(collection);
        collection.clear();
        labelledImages = true;
        break;

      case flip:
        if(flipAxis == 0)
          milx::Image<floatImageType>::FlipCollection(collection, true, false, false, true);
        else if(flipAxis == 1)
          milx::Image<floatImageType>::FlipCollection(collection, false, true, false, true);
        else
          milx::Image<floatImageType>::FlipCollection(collection, false, false, true, true);
        break;

      case none: //--------------------------------
        break;

      default:
        milx::PrintError("Operation not supported. Exiting");
        exit(EXIT_FAILURE);
        break;
    }
    std::cerr << "Done" << std::endl;
  }

  if(outputArg.isSet() && filenames.size() == 1)
  {
    if(vectorImages)
      vecResult = vectorCollection[0];
    if(labelledImages)
      charResult = labelledCollection[0];
    else
      floatResult = collection[0];
  }

  ///Write result, separate as some members convert to/from labelled images
  if(vectorImages)
  {
    if(outputRequired)
    {
      milx::PrintInfo("Writing vector result to " + outputName);
      milx::File::SaveImage<vectorImageType>(outputName, vecResult);
    }
    if(multiOutputRequired)
    {
      if(!vectorCollection.empty())
      {
        const int N = vectorCollection.size();

        milx::PrintInfo("Writing vector results with prefix " + prefixName);
        std::vector<std::string> names;

        int n = 0;
        for(stringiterator filename = filenames.begin(); filename != filenames.end() && n < N; filename ++, n ++)
        {
          std::string ext = milx::File::GetFileExtension(*filename);
          names.push_back(prefixName + milx::File::GetBaseName(*filename) + "." + ext);
          milx::PrintInfo("Will be writing " + prefixName + milx::File::GetBaseName(*filename) + "." + ext);
        }

        milx::File::SaveImages<vectorImageType>(names, vectorCollection);
      }
      else
      {
        milx::PrintError("Result was empty. Skipping Output.");
      }
    }
  }
  else if(labelledImages)
  {
    if(outputRequired)
    {
      milx::PrintInfo("Writing labelled result to " + outputName);
      milx::File::SaveImage<charImageType>(outputName, charResult);
    }
    if(multiOutputRequired)
    {
      if(!labelledCollection.empty())
      {
        const int N = labelledCollection.size();

        milx::PrintInfo("Writing labelled results with prefix " + prefixName);
        std::vector<std::string> names;

        int n = 0;
        for(stringiterator filename = filenames.begin(); filename != filenames.end() && n < N; filename ++, n ++)
        {
          std::string ext = milx::File::GetFileExtension(*filename);
          names.push_back(prefixName + milx::File::GetBaseName(*filename) + "." + ext);
          milx::PrintInfo("Will be writing " + prefixName + milx::File::GetBaseName(*filename) + "." + ext);
        }

        milx::File::SaveImages<charImageType>(names, labelledCollection);
      }
      else
      {
        milx::PrintError("Result was empty. Skipping Output.");
      }
    }
  }
  else
  {
    if(outputRequired)
    {
      milx::PrintInfo("Writing result to " + outputName);
      milx::File::SaveImage<floatImageType>(outputName, floatResult);
    }
    if(multiOutputRequired)
    {
      if(!collection.empty())
      {
        const int N = collection.size();

        milx::PrintInfo("Writing results with prefix " + prefixName);
        std::vector<std::string> names;

        int n = 0;
        for(stringiterator filename = filenames.begin(); filename != filenames.end() && n < N; filename ++, n ++)
        {
          std::string ext = milx::File::GetFileExtension(*filename);
          names.push_back(prefixName + milx::File::GetBaseName(*filename) + "." + ext);
          milx::PrintInfo("Will be writing " + prefixName + milx::File::GetBaseName(*filename) + "." + ext);
        }

        milx::File::SaveImages<floatImageType>(names, collection);
      }
      else
      {
        milx::PrintError("Result was empty. Skipping Output.");
      }
    }
  }

  //---------------------------
  milx::PrintInfo("Operation Complete");
  return EXIT_SUCCESS;
}
