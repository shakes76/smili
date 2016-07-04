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
#include <locale> //isspace
#include <algorithm> //remove_if
#include <tclap/CmdLine.h> //Command line parser library

// SMILI
#include <milxFile.h>
#include <milxImage.h>

using namespace TCLAP;

typedef std::vector< std::string >::iterator stringiterator;
typedef std::vector< std::vector< std::string > >::iterator listiterator;
typedef std::vector< std::pair<std::string, std::string> >::iterator tagiterator;

//Supported operations
enum operations {none = 0, info, convert, print, csv};

//Image typedefs
typedef unsigned char charPixelType;
typedef itk::Image<charPixelType, milx::imgDimension> charImageType;
typedef short shortPixelType;
typedef itk::Image<shortPixelType, milx::imgDimension> shortImageType;
typedef unsigned short ushortPixelType;
typedef itk::Image<ushortPixelType, milx::imgDimension> ushortImageType;
typedef int intPixelType;
typedef itk::Image<intPixelType, milx::imgDimension> intImageType;
typedef unsigned int uintPixelType;
typedef itk::Image<intPixelType, milx::imgDimension> uintImageType;
typedef float floatPixelType;
typedef itk::Image<floatPixelType, milx::imgDimension> floatImageType;
typedef itk::VectorImage<floatPixelType, milx::imgDimension> vectorImageType;

/**
  \file milxDICOMApp.cxx
  \brief A swiss-army knife of DICOM images. This application implements conversion etc. of DICOM images with the SMILI library.
  \ingroup Applications

  You can convert DICOM series, print their tags and get other info.
*/
int main(int argc, char *argv[])
{
  //---------------------------
  ///Program Info
  milx::PrintInfo("--------------------------------------------------------");
  milx::PrintInfo("SMILI DICOM Image Tool for Images/Volumes.");
  milx::PrintInfo("(c) Copyright Chandra et al., 2015.");
  milx::PrintInfo("Version: " + milx::NumberToString(milx::Version));
  milx::PrintInfo("Supported Pixel/Voxel Types: Float, Short, UShort, Int, UInt, UChar, Vector.");
  milx::PrintInfo("Use suffix argument to set output format for conversions. Default is Nifti (.nii.gz)");
  milx::PrintInfo("University of Queensland, Australia.");
  milx::PrintInfo("Australian e-Health Research Centre, CSIRO, Australia.");
  milx::PrintInfo("--------------------------------------------------------\n");

  //---------------------------
  ///Process Arguments
  CmdLine cmd("A diagnostic tool for DICOM series operations", ' ', milx::NumberToString(milx::Version));

  ///Optional
  ValueArg<std::string> outputArg("o", "output", "Output Image", false, "result.nii.gz", "Output");
  ValueArg<std::string> prefixArg("p", "prefix", "Output prefix for multiple output.", false, "img_", "Output Prefix");
  ValueArg<std::string> suffixArg("s", "suffix", "Output suffix for output format and additional text.", false, ".nii.gz", "Output Suffix");
  ValueArg<std::string> exportArg("e", "export-tags", "Output DICOM tags as text pairs to file.", false, "tags.csv", "Export Tags");
  //~ ValueArg<float> decimateArg("d", "decimate", "Decimate all the meshes provided using the Quadric Decimate algorithm", false, 0.5, "Decimate");
  //ValueArg<size_t> mergeArg("", "merge", "Merge the labels in labelled images, 0:keep, 1:aggregate, 2:pack, 3:strict.", false, 0, "Merge");
  ///Switches
  ///XOR Switches
  SwitchArg infoArg("", "info", "Report the image information(s) for each image in the DICOM series.", false);
  SwitchArg convertArg("c", "convert", "Convert the DICOM series to image volumes given at output.", false);
  SwitchArg printArg("t", "print-tags", "Output DICOM tags as text pairs in terminal.", false);
  SwitchArg recursiveArg("r", "recursive", "Recursively parse series directory provided.", false);
  SwitchArg instanceArg("", "instance", "Add Instance ID to filename if found.", false);
  SwitchArg echoArg("", "echo", "Add Echo ID to filename if found.", false);
  SwitchArg acquisitionArg("", "acquisition", "Add Acquisition ID to filename if found.", true);

  ///Mandatory
  UnlabeledMultiArg<std::string> multinames("series", "DICOM Image series to operate on", true, "Series");

  ///Add argumnets
  cmd.add( multinames );
  cmd.add( outputArg );
  cmd.add( prefixArg );
  cmd.add( suffixArg );
  cmd.add( recursiveArg );
  cmd.add( instanceArg );
  cmd.add( echoArg );
  cmd.add( acquisitionArg );
  ///XOR args
  std::vector<Arg*> xorlist;
  xorlist.push_back(&infoArg);
  xorlist.push_back(&convertArg);
  xorlist.push_back(&printArg);
  xorlist.push_back(&exportArg);
  //~ xorlist.push_back(&mseArg);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members

#endif // ITK_REVIEW
  cmd.xorAdd(xorlist);

  ///Parse the argv array.
  cmd.parse( argc, argv );

  ///Get the value parsed by each arg.
  //Filenames of surfaces
  std::vector<std::string> filenames = multinames.getValue();
  std::string outputName = outputArg.getValue();
  const std::string prefixName = prefixArg.getValue();
  const std::string suffixName = suffixArg.getValue();
  const std::string exportName = exportArg.getValue();
  const bool recurse = recursiveArg.isSet();

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
  if(convertArg.isSet())
  {
    if(outputArg.isSet())
    {
      milx::PrintError("Argument Error: Use Prefix (-p) for conversion.");
      milx::PrintError("Re-run with the prefix name set.");
      exit(EXIT_FAILURE);
    }

    outputRequired = true;
    multiOutputRequired = false;

    //Convert
    operation = convert;
  }
  if(exportArg.isSet())
  {
    //Export Tags
    operation = csv;
  }
  if(printArg.isSet())
  {
      //Print Tags
      operation = print;
  }

  std::cout << "Total Folders: " << filenames.size() << std::endl;
  if(filenames.empty())
  {
    milx::PrintError("No folders provided. Exiting.");
    exit(EXIT_FAILURE);
  }
  for(stringiterator name = filenames.begin(); name != filenames.end(); name ++)
    std::cout << *name << ", ";
  std::cout << std::endl;

  //---------------------------
  ///Get UIDs and filenames
  std::cerr << "Reading series UIDs" << std::endl;
  std::vector< std::vector<std::string> > UIDList;
  std::vector<std::string> validFilenames;
  for (stringiterator name = filenames.begin(); name != filenames.end(); name++)
  {
    std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(*name, recurse);
    milx::PrintInfo("Found " + milx::NumberToString(UIDs.size()) + " in " + *name);
    if (!UIDs.empty())
    {
      UIDList.push_back(UIDs);
      validFilenames.push_back(*name);
    }
  }
  std::cerr << "Done" << std::endl;

  if (UIDList.empty())
  {
    milx::PrintError("Found no series found in input directory. You may need to provide internal directory paths explicitly.");
    return EXIT_FAILURE;
  }

  stringiterator dir;
  listiterator list;
  for (dir = validFilenames.begin(), list = UIDList.begin(); list != UIDList.end(); list++, dir++)
  {
    std::cerr << "Applying operation to " << *dir << std::endl;
    for (stringiterator name = list->begin(); name != list->end(); name++)
    {
      std::cerr << "Processing UID: " << *name << std::endl;
      const std::vector<std::string> seriesFilenames = milx::File::GetDICOMSeriesFilenames(*dir, *name, recurse);

      //Read Header
      size_t dimensions = 3;
      std::string pixelType, componentType;
      if (!milx::File::ReadImageInformation(seriesFilenames.front(), pixelType, componentType, dimensions))
      {
        milx::PrintError("Failed Reading First Image. Check the image type/file. Exiting.");
        exit(EXIT_FAILURE);
      }
      if (infoArg.isSet())
      {
        milx::PrintInfo("UID: " + *name);
        milx::PrintInfo("Pixel Type: " + pixelType);
        milx::PrintInfo("Component Type: " + componentType);
        milx::PrintInfo("Dimensions: " + milx::NumberToString(dimensions));
      }

      //Open image using relevant type
      std::string caseID;
      std::string echoID = "";
      std::string seriesID = "";
      std::string acqID = "";
      std::string instanceID = "";
      itk::SmartPointer<vectorImageType> vectorImage;
      itk::SmartPointer<charImageType> labelledImage;
      itk::SmartPointer<shortImageType> shortImage;
      itk::SmartPointer<ushortImageType> ushortImage;
      itk::SmartPointer<intImageType> intImage;
      itk::SmartPointer<uintImageType> uintImage;
      itk::SmartPointer<floatImageType> floatImage;
      std::vector< std::pair<std::string, std::string> > tags;
      bool labelledImages = false, shortImages = false, ushortImages = false, integerImages = false, uintegerImages = false, vectorImages = false;
      if (pixelType == "vector" || dimensions > 3) ///\todo handle 4D images here properly
      {
        milx::PrintInfo("Detected vector images.");
        if (!milx::File::OpenDICOMSeriesAndTags<vectorImageType>(*dir, vectorImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
        {
          milx::PrintError("Failed Reading Vector Images. Exiting.");
          exit(EXIT_FAILURE);
        }
        vectorImages = true;
      }
      else if (componentType == "unsigned_char" || componentType == "unsigned char")
      {
        milx::PrintInfo("Detected labelled images.");
        if (!milx::File::OpenDICOMSeriesAndTags<charImageType>(*dir, labelledImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
        {
          milx::PrintError("Failed Reading Labelled Images. Exiting.");
          exit(EXIT_FAILURE);
        }
        labelledImages = true;
      }
      else if (componentType == "short" || componentType == "int16")
      {
          milx::PrintInfo("Detected short images.");
          if (!milx::File::OpenDICOMSeriesAndTags<shortImageType>(*dir, shortImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
          {
              milx::PrintError("Failed Reading Short Images. Exiting.");
              exit(EXIT_FAILURE);
          }
          shortImages = true;
      }
      else if (componentType == "unsigned_short" || componentType == "unsigned short")
      {
          milx::PrintInfo("Detected unsigned short images.");
          if (!milx::File::OpenDICOMSeriesAndTags<ushortImageType>(*dir, ushortImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
          {
              milx::PrintError("Failed Reading Unsigned Short Images. Exiting.");
              exit(EXIT_FAILURE);
          }
          ushortImages = true;
      }
      else if (componentType == "int" || componentType == "signed" || componentType == "int32" || componentType == "int64")
      {
        milx::PrintInfo("Detected integer images.");
        if (!milx::File::OpenDICOMSeriesAndTags<intImageType>(*dir, intImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
        {
          milx::PrintError("Failed Reading Integer Images. Exiting.");
          exit(EXIT_FAILURE);
        }
        integerImages = true;
      }
      else if (componentType == "unsigned_int" || componentType == "unsigned int" || componentType == "unsigned")
      {
        milx::PrintInfo("Detected unsigned int images.");
        if (!milx::File::OpenDICOMSeriesAndTags<uintImageType>(*dir, uintImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
        {
          milx::PrintError("Failed Reading Unsigned Integer Images. Exiting.");
          exit(EXIT_FAILURE);
        }
        uintegerImages = true;
      }
      else
      {
        milx::PrintInfo("Detected floating point images.");
        if (!milx::File::OpenDICOMSeriesAndTags<floatImageType>(*dir, floatImage, tags, *name, caseID, echoID, seriesID, acqID, instanceID, recurse)) //Error NOT printed inside
        {
          milx::PrintError("Failed Reading Images. Exiting.");
          exit(EXIT_FAILURE);
        }
      }

      //Remove illegal characters from names
      name->erase(std::remove(name->begin(),name->end(),' '),name->end());
      caseID.erase(std::remove(caseID.begin(),caseID.end(),' '),caseID.end());
      if(!echoID.empty())
        echoID.erase(std::remove(echoID.begin(),echoID.end(),' '),echoID.end());
      if(!seriesID.empty())
        seriesID.erase(std::remove(seriesID.begin(),seriesID.end(),' '),seriesID.end());
      if(!acqID.empty())
        acqID.erase(std::remove(acqID.begin(),acqID.end(),' '),acqID.end());
      if(!instanceID.empty())
        instanceID.erase(std::remove(instanceID.begin(),instanceID.end(),' '),instanceID.end());

      //Create directories
      std::string path = caseID + "/" + *name;
      if(!echoID.empty())
        path += "_" + echoID;
      if (prefixArg.isSet())
        path = prefixName + path;
      milx::PrintDebug("Making Directory " + path);
      milx::File::MakeDirectory(path);

      //create filename if needed
      std::string filename = path + "/" + caseID + "_" + *name;
      if(!echoID.empty() && echoArg.isSet())
        filename += "_" + echoID;
      if(!acqID.empty() && acquisitionArg.isSet())
        filename += "_Acq_" + acqID;
      if(!seriesID.empty())
        filename += "_Series_" + seriesID;
      if(!instanceID.empty() && instanceArg.isSet())
        filename += "_Instance_" + instanceID;
      filename += suffixName;
      milx::PrintDebug("Output Filename would be " + filename);

      if(operation == print) //independent operation from pixel type
      {
        std::cout << "Tag      \t| Value" << std::endl;
        for (tagiterator tag = tags.begin(); tag != tags.end(); tag++)
          std::cout << tag->first << "\t| " << tag->second << std::endl;
      }
      else if(operation == csv) //independent operation from pixel type
      {
          ofstream outFile(exportName.c_str(), ios::app);
          //outFile << "Tag,Value" << std::endl;
          for (tagiterator tag = tags.begin(); tag != tags.end(); tag++)
              outFile << tag->first << "," << tag->second << std::endl;
          outFile.close();
      }

      //std::cout << "Series: " << *name << std::endl;
      //std::cout << "Case ID: " << caseID << std::endl;
      if (vectorImages)
      {
        switch (operation)
        {
        case info:
          milx::Image<vectorImageType>::Information(vectorImage);
          break;

        case convert:
          milx::File::SaveImage<vectorImageType>(filename, vectorImage);
          break;

        case print:
          break;
        case csv:
          break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      else if (labelledImages)
      {
        switch (operation)
        {
        case info:
          milx::Image<charImageType>::Information(labelledImage);
          break;

        case convert:
          milx::File::SaveImage<charImageType>(filename, labelledImage);
          break;

        case print:
          break;
        case csv:
          break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      else if (shortImages)
      {
        switch (operation)
        {
        case info:
          milx::Image<shortImageType>::Information(shortImage);
          break;

        case convert:
          milx::File::SaveImage<shortImageType>(filename, shortImage);
          break;

        case print:
          break;
        case csv:
          break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      else if (ushortImages)
      {
        switch (operation)
        {
        case info:
          milx::Image<ushortImageType>::Information(ushortImage);
          break;

        case convert:
          milx::File::SaveImage<ushortImageType>(filename, ushortImage);
          break;

        case print:
          break;
        case csv:
          break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      else if (integerImages)
      {
        switch (operation)
        {
        case info:
          milx::Image<intImageType>::Information(intImage);
          break;

        case convert:
          milx::File::SaveImage<intImageType>(filename, intImage);
          break;

        case print:
          break;
        case csv:
            break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      else if (uintegerImages)
      {
        switch (operation)
        {
        case info:
          milx::Image<uintImageType>::Information(uintImage);
          break;

        case convert:
          milx::File::SaveImage<uintImageType>(filename, uintImage);
          break;

        case print:
          break;
        case csv:
            break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      else
      {
        switch (operation)
        {
        case info:
          milx::Image<floatImageType>::Information(floatImage);
          break;

        case convert:
          milx::File::SaveImage<floatImageType>(filename, floatImage);
          break;

        case print:
          break;
        case csv:
            break;

        case none: //--------------------------------
          break;

        default:
          milx::PrintError("Operation not supported. Exiting");
          exit(EXIT_FAILURE);
          break;
        }
      }
      std::cerr << "Done" << std::endl;
    }
    std::cerr << "Done" << std::endl;
  }

  //---------------------------
  milx::PrintInfo("Operation Complete");
  return EXIT_SUCCESS;
}
