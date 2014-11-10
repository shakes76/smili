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
#include "milxQtFile.h"

#include <zlib.h>

#include <QFileInfo>
//VTK
#include <vtkImageFlip.h>
#include <vtkImageCast.h>
#include <vtkOBJReader.h>
#include <vtkOBJExporter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPNMReader.h>
///ITK Imaging
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkPNGImageIOFactory.h>
#include <itkJPEGImageIOFactory.h>
#include <itkBMPImageIOFactory.h>
#include <itkNrrdImageIOFactory.h>
#include <itkRawImageIO.h>

#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"

milxQtFile::milxQtFile(QObject *theParent) : QFile(theParent)
{
    //Allocate critical variables
    observeProgress = itkEventQtObserver::New();
    Connector = vtkSmartPointer<vtkEventQtSlotConnect>::New();

    name = "Empty";
    dataPixelType = "scalar";
    dataComponentType = "float";
    dataDimensions = milx::imgDimension;
    dataComponents = 1;

    ///Connect Image Object for progress updates
    linkProgressEventOf(milx::ProgressUpdates->GetUpdateObject()); //link itk observer propagator
}

milxQtFile::~milxQtFile()
{
    //dtor
}

bool milxQtFile::openImage(const QString filename, vtkImageData* data)
{
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();
    bool integerFormat = false, vtkFormat = false, medical = true;
    int bounds[6];

    const QString charStr = "unsigned char";
    const QString typeStr = data->GetScalarTypeAsString();

    if(extension == "png" || extension == "jpg" || extension == "jpeg" || extension == "bmp")
    {
        integerFormat = true;
        medical = false;
    }
    else if(typeStr == charStr)
    {
        integerFormat = true;
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char, 3>::New() );
    }
    else if(extension == "vti")
    {
        vtkFormat = true;
    }
    else
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<float, 3>::New() );

    ///Flip - ITK and VTK do not have the same orientation
    vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();

    if(integerFormat)
    {
        if(!medical)
        {
          //Add some default image types
          itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char,2>::New() );
          itk::ObjectFactoryBase::RegisterFactory( itk::PNGImageIOFactory::New() );
          itk::ObjectFactoryBase::RegisterFactory( itk::JPEGImageIOFactory::New() );
          itk::ObjectFactoryBase::RegisterFactory( itk::BMPImageIOFactory::New() );
        }

        charImageType::Pointer charImg = milx::File::ReadImageUsingITK<charImageType>(filename.toStdString());

        if(!charImg)
            return false;

        ///Export to VTK and flip
    #if VTK_MAJOR_VERSION <=5
        imageReorient->SetInput(milx::Image<charImageType>::ConvertITKImageToVTKImage(charImg));
    #else
        imageReorient->SetInputData(milx::Image<charImageType>::ConvertITKImageToVTKImage(charImg));
    #endif // VTK_MAJOR_VERSION
    }
    else if(vtkFormat)
    {
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();

        if(reader->CanReadFile(filename.toStdString().c_str()))
        {
            reader->SetFileName(filename.toStdString().c_str());
            reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(reader);
            reader->Update();

            if(!errorObserver->ReportsFailure())
            {
                data->DeepCopy(reader->GetOutput());
                data->GetExtent(bounds);
            }
            else
            {
                cerr << "VTI Reader Encountered the following error." << endl;
                cerr << errorObserver->GetMessage() << endl;
                return false;
            }
        }
        else
        {
            cerr << "Could not load VTI file!" << endl;
            return false;
        }
    }
    else
    {
        floatImageType::Pointer floatImg = milx::File::ReadImageUsingITK<floatImageType>(filename.toStdString());

        if(!floatImg)
            return false;

        ///Export to VTK and flip
    #if VTK_MAJOR_VERSION <=5
        imageReorient->SetInput(milx::Image<floatImageType>::ConvertITKImageToVTKImage(floatImg));
    #else
        imageReorient->SetInputData(milx::Image<floatImageType>::ConvertITKImageToVTKImage(floatImg));
    #endif // VTK_MAJOR_VERSION
    }

    if(!vtkFormat)
    {
        imageReorient->SetFilteredAxis(1);
        imageReorient->FlipAboutOriginOn();

        linkProgressEventOf(imageReorient);
        imageReorient->Update();
        data->DeepCopy(imageReorient->GetOutput());
    }

    return true;
}

QString milxQtFile::supportedImageFormats()
{
  QString exts = "";
  vector<string> extensions = milx::File::GetSupportedImageFileExtensions();

  for(size_t j = 0; j < extensions.size(); j ++)
    exts += "*" + QString(extensions[j].c_str()) + " ";

  return exts;
}

bool milxQtFile::isIntegerFormat(const QString filename, bool &errorEncountered)
{
  std::string pixelType, componentType;
  errorEncountered = false;

  //Check type of medical image
  if(!milx::File::ReadImageInformation(filename.toStdString(), pixelType, componentType, dataDimensions))
  {
      cerr << "Failed reading header of image. File may not be an image. Exiting" << endl;
      errorEncountered = true;
      return false;
  }
  dataPixelType = pixelType.c_str();
  dataComponentType = componentType.c_str();

  if(componentType == "unsigned_char")
      return true;

  return false;
}

bool milxQtFile::isFieldFormat(const QString filename, bool &errorEncountered)
{
  std::string pixelType, componentType;
  errorEncountered = false;

  //Check type of medical image
  if(!milx::File::ReadImageInformation(filename.toStdString(), pixelType, componentType, dataDimensions))
  {
      cerr << "Failed reading header of image. File may not be an image. Exiting" << endl;
      errorEncountered = true;
      return false;
  }
  dataPixelType = pixelType.c_str();
  dataComponentType = componentType.c_str();

  if(pixelType == "vector")
      return true;

  return false;
}

bool milxQtFile::openImage(const QString filename, milxQtImage* data)
{
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();
    bool integerFormat = false, vtkFormat = false, deformField = false, rgbImage = false, pnmImage = false;

    if(extension == "png" || extension == "jpg" || extension == "jpeg" || extension == "bmp")
    {
        //Add some default image types
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char,2>::New() );
        itk::ObjectFactoryBase::RegisterFactory( itk::PNGImageIOFactory::New() );
        itk::ObjectFactoryBase::RegisterFactory( itk::JPEGImageIOFactory::New() );
        itk::ObjectFactoryBase::RegisterFactory( itk::BMPImageIOFactory::New() );
    }
    else if(extension == "vti")
    {
        vtkFormat = true;
    }
    else if(extension == "pbm" || extension == "pgm" || extension == "ppm")
    {
        pnmImage = true;
    }

    if(!vtkFormat && !pnmImage)
    {
        cerr << "Trying to read image header ..." << endl;
        //Check type of medical image
        std::string pixelType, componentType;
        if(!milx::File::ReadImageInformation(filename.toStdString(), pixelType, componentType, dataDimensions))
        {
            cerr << "Failed reading header of image. File may not be an image. Exiting" << endl;
            return false;
        }
        dataPixelType = pixelType.c_str();
        dataComponentType = componentType.c_str();

        if(componentType == "unsigned_char" && pixelType == "scalar")
            integerFormat = true;
        if(componentType == "unsigned_short") //PNG etc. maybe 16-bit integers
            integerFormat = false;
        if(pixelType == "vector")
            deformField = true;
        if( (pixelType == "rgb" || pixelType == "rgba") && componentType == "unsigned_char" )
            rgbImage = true;

        data->setActualNumberOfDimensions(dataDimensions);
    }

//    cerr << "Open Image" << endl;
    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    if(integerFormat)
    {
        charImageType::Pointer charImg;

        if( !milx::File::OpenImage<charImageType>(filename.toStdString(), charImg) )
            return false;

        data->SetInput(charImg, true);
    }
    else if(vtkFormat)
    {
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();

        if(reader->CanReadFile(filename.toStdString().c_str()))
        {
            reader->SetFileName(filename.toStdString().c_str());
            reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(reader);
            reader->Update();

            if(!errorObserver->ReportsFailure())
                data->SetInput(reader->GetOutput());
            else
            {
                cerr << "VTI Reader Encountered the following error." << endl;
                cerr << errorObserver->GetMessage() << endl;
                return false;
            }
        }
        else
        {
            cerr << "Could not load VTI file!" << endl;
            return false;
        }
    }
    else if(pnmImage)
    {
        vtkSmartPointer<vtkPNMReader> reader = vtkSmartPointer<vtkPNMReader>::New();

        if(reader->CanReadFile(filename.toStdString().c_str()))
        {
            reader->SetFileName(filename.toStdString().c_str());
            reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(reader);
            reader->Update();

            if(!errorObserver->ReportsFailure())
            {
                cout << "Image Description: " << reader->GetDescriptiveName() << endl;
                data->SetInput(reader->GetOutput());
            }
            else
            {
                cerr << "PNM Reader Encountered the following error." << endl;
                cerr << errorObserver->GetMessage() << endl;
                return false;
            }
        }
        else
        {
            cerr << "Could not load PNM file!" << endl;
            return false;
        }
    }
    else if(deformField) //Read vector images
    {
        itk::ObjectFactoryBase::RegisterFactory( itk::NrrdImageIOFactory::New() );

        vectorImageType::Pointer vecImg;

        if( !milx::File::OpenImage<vectorImageType>(filename.toStdString(), vecImg) )
            return false;

        dataComponents = vecImg->GetNumberOfComponentsPerPixel();

        data->SetInput(vecImg, true);
    }
#if (ITK_VERSION_MAJOR > 3)
    else if(dataDimensions > 3)
    {
        vectorImageType::Pointer vecImg;

        typedef itk::Image<floatPixelType, 4> float4DImageType;
        if( !milx::File::OpenAsVectorImage<float4DImageType, 3>(filename.toStdString(), vecImg) )
            return false;

        dataComponents = vecImg->GetNumberOfComponentsPerPixel();

        data->SetInput(vecImg, true);
    }
#endif
    else if(rgbImage)
    {
        rgbImageType::Pointer rgbImg;

        if( !milx::File::OpenImage<rgbImageType>(filename.toStdString(), rgbImg) )
            return false;

        data->SetInput(rgbImg, true);
    }
    else
    {
        floatImageType::Pointer floatImg;

        if( !milx::File::OpenImage<floatImageType>(filename.toStdString(), floatImg) )
            return false;

        data->SetInput(floatImg, true);
    }

    return true;
}

bool milxQtFile::openImageSeries(milxQtImage* data, QString directoryPath)
{
  bool success = true;
  QPointer<QFileDialog> fileOpener = new QFileDialog;
  QSettings settings("Shekhar Chandra", "milxQt");

  ///If filenames list is empty ask for them
  if(directoryPath.isEmpty())
    {
      QString path = settings.value("recentPath").toString();
      directoryPath = fileOpener->getExistingDirectory(NULL, tr("Open DICOM Directory"),
                                           path,
                                           QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    }

  if(directoryPath.isEmpty())
    return !success;

  ///Get UIDs
  std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(directoryPath.toStdString());

  std::string seriesName;
  if(UIDs.size() > 1) //ask user which one to load
  {
      QLabel *uidLabel = new QLabel;
      uidLabel->setText("<b>Select Series</b>");
      QLabel *uidInfoLabel = new QLabel;
      uidInfoLabel->setText("Multiple DICOM Series Found. Please choose one:");
      QComboBox *uidComboBox = new QComboBox;
      for(size_t j = 0; j < UIDs.size(); j ++)
          uidComboBox->insertItem(j, UIDs[j].c_str());
      QPushButton *uidButton = new QPushButton;
      uidButton->setText("OK");

      QVBoxLayout *uidLayout = new QVBoxLayout;
      uidLayout->addWidget(uidLabel);
      uidLayout->addWidget(uidInfoLabel);
      uidLayout->addWidget(uidComboBox);
      uidLayout->addWidget(uidButton);

      QDialog *uidWidget = new QDialog;
      uidWidget->setWindowTitle("Select Series to Load");
      uidWidget->setLayout(uidLayout);
      uidWidget->setModal(true);
      QObject::connect(uidButton, SIGNAL(clicked()), uidWidget, SLOT(accept()));
      uidWidget->exec();

      seriesName = uidComboBox->currentText().toStdString();
      delete uidWidget;
  }
  else if(UIDs.empty())
  {
      cerr << "Error. No DICOM series was found in directory" << endl;
      return false;
  }
  else
      seriesName = UIDs.begin()->c_str();

  cout << "Reading series as float images" << endl;
  floatImageType::Pointer floatImg;
  milx::File::OpenDICOMSeries<floatImageType>(directoryPath.toStdString(), floatImg, seriesName);
  data->SetInput(floatImg, false);
  data->setName(seriesName.c_str());
  cout << "Completed Reading Series: " << seriesName << endl;

  //save path
  QFileInfo fi(directoryPath);
  settings.setValue("recentPath",fi.absolutePath());

  return success;
}

bool milxQtFile::saveImage(const QString filename, vtkImageData* data)
{
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();
    bool integerFormat = false, medical = true, vtkFormat = false, success = false;
    int bounds[6];

    const QString charStr = "unsigned char";
    const QString typeStr = data->GetScalarTypeAsString();

    if(extension == "png" || extension == "jpg" || extension == "jpeg" || extension == "bmp")
    {
        integerFormat = true;
        medical = false;
    }
    else if(typeStr == charStr)
    {
        integerFormat = true;
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char, 3>::New() );
    }
    else
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<float, 3>::New() );

    data->GetExtent(bounds);

    if(extension == "vti")
        vtkFormat = true;

    ///Flip - ITK and VTK do not have the same orientation
    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    #if VTK_MAJOR_VERSION <=5
        imageReorient->SetInputConnection(data->GetProducerPort());
    #else
        imageReorient->SetInputData(data);
    #endif
        imageReorient->SetFilteredAxis(1);
        imageReorient->FlipAboutOriginOn();
        linkProgressEventOf(imageReorient);
        imageReorient->Update();

    if(integerFormat)
    {
        if(!medical)
        {
            //Add some default image types
            itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char, 2>::New() );
            itk::ObjectFactoryBase::RegisterFactory( itk::PNGImageIOFactory::New() );
            itk::ObjectFactoryBase::RegisterFactory( itk::JPEGImageIOFactory::New() );
            itk::ObjectFactoryBase::RegisterFactory( itk::BMPImageIOFactory::New() );
//
            ///Export to ITK
            typedef itk::Image<rgbPixelType, 2> rgbImageSliceType;
            rgbImageSliceType::Pointer charImg = milx::Image<rgbImageSliceType>::ConvertVTKImageToITKImage(imageReorient->GetOutput());
            success = milx::File::WriteImageUsingITK<rgbImageSliceType>(filename.toStdString(), charImg);
        }
        else
        {
            ///Export to ITK
            rgbImageType::Pointer charImg = milx::Image<rgbImageType>::ConvertVTKImageToITKImage(imageReorient->GetOutput());
            success = milx::File::WriteImageUsingITK<rgbImageType>(filename.toStdString(), charImg);
        }
    }
    else if(vtkFormat)
    {
        vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data);
        #else
            writer->SetInputData(data);
        #endif
            writer->SetFileName(filename.toStdString().c_str());
            writer->SetDataModeToBinary();
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();

            if(errorObserver->ReportsFailure())
            {
                cerr << "VTI Writer Encountered the following error." << endl;
                cerr << errorObserver->GetMessage() << endl;
            }
            else
                success = true;
    }
    else
    {
        ///Export to ITK
        floatImageType::Pointer floatImg = milx::Image<floatImageType>::ConvertVTKImageToITKImage(imageReorient->GetOutput());

        success = milx::File::WriteImageUsingITK<floatImageType>(filename.toStdString(), floatImg);
    }

    return success;
}

bool milxQtFile::saveImage(const QString filename, milxQtImage* data)
{
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();
    bool integerFormat = false, medical = true, rgbFormat = false, vtkFormat = false, success = false;

    if(extension == "png" || extension == "jpg" || extension == "jpeg" || extension == "bmp")
    {
        integerFormat = true;
        medical = false;
    }
    else if(data->is8BitImage())
    {
        integerFormat = true;
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char,3>::New() );
    }
    else if(data->isRGBImage())
    {
        rgbFormat = true;
    }
    else
        itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<float,3>::New() );

    if(extension == "vti")
        vtkFormat = true;
    else if(data->isVTKImage())
    {
        vtkFormat = false;
        vtkSmartPointer<vtkImageCast> castVTKImage = vtkSmartPointer<vtkImageCast>::New();
            castVTKImage->SetOutputScalarTypeToFloat();
        #if VTK_MAJOR_VERSION <=5
            castVTKImage->SetInput(data->GetOutput());
        #else
            castVTKImage->SetInputData(data->GetOutput());
        #endif // VTK_MAJOR_VERSION
            castVTKImage->Update();

        vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
        #if VTK_MAJOR_VERSION <=5
            imageReorient->SetInput(castVTKImage->GetOutput());
        #else
            imageReorient->SetInputData(castVTKImage->GetOutput());
        #endif
            imageReorient->SetFilteredAxis(1);
            imageReorient->FlipAboutOriginOn();
            linkProgressEventOf(imageReorient);
            imageReorient->Update();

        cout << "Converted VTK Image to ITK Image since saving requested medical image format" << endl;
        floatImageType::Pointer ITKImage = milx::Image<floatImageType>::ConvertVTKImageToITKImage(imageReorient->GetOutput());
        data->SetInput(ITKImage);
        data->generateImage();
    }

    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    if(vtkFormat)
    {
        vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data->GetOutput());
        #else
            writer->SetInputData(data->GetOutput());
        #endif // VTK_MAJOR_VERSION
            writer->SetFileName(filename.toStdString().c_str());
            writer->SetDataModeToBinary();
            linkProgressEventOf(writer);
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            writer->Write();

            if(errorObserver->ReportsFailure())
            {
                cerr << "VTI Writer Encountered the following error." << endl;
                cerr << errorObserver->GetMessage() << endl;
            }
            else
                success = true;
    }
    else if(integerFormat)
    {
        if(!medical)
        {
            //Add some default image types
            itk::ObjectFactoryBase::RegisterFactory( itk::RawImageIOFactory<unsigned char, 2>::New() );
            itk::ObjectFactoryBase::RegisterFactory( itk::PNGImageIOFactory::New() );
            itk::ObjectFactoryBase::RegisterFactory( itk::JPEGImageIOFactory::New() );
            itk::ObjectFactoryBase::RegisterFactory( itk::BMPImageIOFactory::New() );
        }

        success = milx::File::WriteImageUsingITK<charImageType>(filename.toStdString(), data->GetCharImage());
    }
    else if(rgbFormat)
    {
        success = milx::File::WriteImageUsingITK<rgbImageType>(filename.toStdString(), data->GetRGBImage());
    }
    else if(data->isVectorImage())
    {
        success = milx::File::WriteImageUsingITK<vectorImageType>(filename.toStdString(), data->GetVectorImage());
    }
    else
    {
        success = milx::File::WriteImageUsingITK<floatImageType>(filename.toStdString(), data->GetFloatImage());
    }

    return success;
}

bool milxQtFile::openModel(const QString filename, vtkPolyData* data)
{
    bool legacy = false, wavefront = false, stanford = false, stereoLith = false;
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();

    if(extension == "vtk")
        legacy = true; //Load legacy VTK file
    else if(extension == "obj")
        wavefront = true;
    else if(extension == "ply")
        stanford = true;
    else if(extension == "stl")
        stereoLith = true;

    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    if(legacy)
    {
        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->DeepCopy(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else if(wavefront)
    {
        vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->DeepCopy(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else if(stanford)
    {
        vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->DeepCopy(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else if(stereoLith)
    {
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->DeepCopy(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();
        if(!errorObserver->ReportsFailure())
            data->DeepCopy(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }

    return true;
}

bool milxQtFile::openModel(const QString filename, milxQtModel* data)
{
    bool legacy = false, wavefront = false, stanford = false, stereoLith = false;
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();

    if(extension == "vtk")
        legacy = true; //Load legacy VTK file
    else if(extension == "obj")
        wavefront = true;
    else if(extension == "ply")
        stanford = true;
    else if(extension == "stl")
        stereoLith = true;

    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    if(legacy)
    {
        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->SetInput(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else if(wavefront)
    {
        vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->SetInput(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else if(stanford)
    {
        vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->SetInput(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else if(stereoLith)
    {
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->SetInput(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        linkProgressEventOf(reader);
        reader->Update();

        if(!errorObserver->ReportsFailure())
            data->SetInput(reader->GetOutput());
        else
        {
            cerr << "Reader Encountered the following error." << endl;
            cerr << errorObserver->GetMessage() << endl;
            return false;
        }
    }

    return true;
}

bool milxQtFile::openModelCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{
    bool success = true;
    QPointer<QFileDialog> fileOpener = new QFileDialog;
    QSettings settings("Shekhar Chandra", "milxQt");

    ///If filenames list is empty ask for them
    if(filenames.isEmpty())
    {
        QString path = settings.value("recentPath").toString();
        filenames = fileOpener->getOpenFileNames(NULL,
                    tr("Select File(s) to Open"),
                    path,
                    tr(openModelExts.c_str()));
    }

    if(filenames.size() == 0)
        return !success;

    QFileInfo fi(filenames[0]);
    settings.setValue("recentPath", fi.absolutePath()); //!< Remember path if successful.

    for(int j = 0; j < filenames.size(); j ++)
    {
        vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();

        success = openModel(filenames[j], data); ///Open supported model file

        if(!success)
        {
            cerr << "Encountered Error in Reading model. Aborting Collection Read." << endl;
            break;
        }
        else
            cout << "Opened " << filenames[j].toStdString() << " into collection." << endl;

        collection->AddItem(data);

        qApp->processEvents();
    }

    return success;
}

bool milxQtFile::saveModel(const QString filename, vtkPolyData* data, const bool binary)
{
    bool legacy = false, stanford = false, stereoLith = false;
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();

    if(extension == "vtk")
        legacy = true; //Save legacy VTK file
    else if(extension == "ply")
        stanford = true; //Save Stanford Poly file
    else if(extension == "stl")
        stereoLith = true; //STL files

    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    if(legacy)
    {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data);
        #else
            writer->SetInputData(data);
        #endif // VTK_MAJOR_VERSION
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else if(stanford)
    {
        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data);
        #else
            writer->SetInputData(data);
        #endif // VTK_MAJOR_VERSION
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else if(stereoLith)
    {
        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data);
        #else
            writer->SetInputData(data);
        #endif // VTK_MAJOR_VERSION
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data);
        #else
            writer->SetInputData(data);
        #endif // VTK_MAJOR_VERSION
            if(binary)
                writer->SetDataModeToBinary();
            writer->SetCompressorTypeToZLib();
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }

    if(errorObserver->ReportsFailure())
    {
        cerr << "Writer Encountered the following error." << endl;
        cerr << errorObserver->GetMessage() << endl;
        return false;
    }

    return true;
}

bool milxQtFile::saveModel(const QString filename, milxQtModel* data, const bool binary)
{
    bool legacy = false, stanford = false, wavefront = false, stereoLith = false;
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();

    if(extension == "vtk")
        legacy = true; //Save legacy VTK file
    else if(extension == "ply")
        stanford = true; //Save Stanford Poly file
    else if(extension == "obj")
        wavefront = true;
    else if(extension == "stl")
        stereoLith = true; //STL files

    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    if(legacy)
    {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data->GetOutput());
        #else
            writer->SetInputData(data->GetOutput());
        #endif // VTK_MAJOR_VERSION
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else if(stanford)
    {
        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data->GetOutput());
        #else
            writer->SetInputData(data->GetOutput());
        #endif // VTK_MAJOR_VERSION
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else if(wavefront)
    {
        vtkSmartPointer<vtkOBJExporter> writer = vtkSmartPointer<vtkOBJExporter>::New();
            QString namePrefix = fileInfo.path() + "/" + fileInfo.baseName();
            cout << "Exporting with prefix " << namePrefix.toStdString().c_str() << endl;
            writer->SetFilePrefix(namePrefix.toStdString().c_str());
            data->disableOrient();
            writer->SetInput(data->GetRenderWindow());
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else if(stereoLith)
    {
        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data->GetOutput());
        #else
            writer->SetInputData(data->GetOutput());
        #endif // VTK_MAJOR_VERSION
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
        #if VTK_MAJOR_VERSION <=5
            writer->SetInput(data->GetOutput());
        #else
            writer->SetInputData(data->GetOutput());
        #endif // VTK_MAJOR_VERSION
            if(binary)
                writer->SetDataModeToBinary();
            writer->SetCompressorTypeToZLib();
            writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            linkProgressEventOf(writer);
            writer->Write();
    }

    if(errorObserver->ReportsFailure())
    {
        cerr << "Writer Encountered the following error." << endl;
        cerr << errorObserver->GetMessage() << endl;
        return false;
    }

    return true;
}

bool milxQtFile::saveScalarsOfModel(const QString filename, milxQtModel* data)
{
    bool csvOutput = false;
    QFileInfo fileInfo(filename);
    QString extension = fileInfo.suffix().toLower();

    if(extension == "csv" || extension == "txt")
        csvOutput = true; //Save csv file

    if(csvOutput)
    {
        QFile::setFileName(filename);
        if(!QFile::open(QIODevice::WriteOnly | QIODevice::Text))
            return false;

        QTextStream outFile(this);
        vtkSmartPointer<vtkDataArray> scalars = data->GetOutput()->GetPointData()->GetScalars();

//        for(int k = 0; k < scalars->GetNumberOfComponents(); k ++)
//        {
            for(int j = 0; j < scalars->GetNumberOfTuples(); j ++)
            {
                outFile << scalars->GetTuple1(j);
                if(j < scalars->GetNumberOfTuples()-1)
                   outFile << ", ";
            }
//            outFile << "\n";
//        }

        QFile::close();
    }
    else
    {

    }

    return true;
}

void milxQtFile::linkProgressEventOf(vtkObject * obj)
{
    Connector->Connect(obj,
                       vtkCommand::ProgressEvent,
                       this,
                       SLOT( updateQtEvents() ),
                       NULL, 1.0); //High Priority
}
