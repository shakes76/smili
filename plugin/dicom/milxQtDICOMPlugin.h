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
#ifndef MILXQTDICOMPLUGIN_H
#define MILXQTDICOMPLUGIN_H

#include <QThread>
#include <QMenu>
#include <QDockWidget>
#include <QWizardPage>
#include <QGroupBox>

//ITK
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkPolygonSpatialObject.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#if (ITK_VERSION_MAJOR > 3) //Review only members
    #include <gdcmTypes.h>
    #include <gdcmAttribute.h>
    #include <gdcmReader.h>
#endif // (ITK_VERSION_MAJOR > 3)
//VTK
#include <vtkPolyDataCollection.h>
//milxSMILI
#include "milxFile.h"
//milxQt
#include "milxQtAliases.h"
#include "milxQtPluginInterface.h"
#include "milxQtImage.h"
#include "milxQtMain.h"
#include "milxQtManager.h"

/**
    \class milxQtDICOMPlugin
    \brief The interface for the DICOM plugin for milxQt
    \author 

    
*/
class MILXQT_PLUGIN_EXPORT milxQtDICOMPlugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    typedef itk::GDCMImageIO                        ImageIOType;
    typedef itk::MetaDataDictionary                 DictionaryType;
    typedef itk::MetaDataObject< std::string >      MetaDataStringType;
  
    /**
        \fn milxQtDICOMPlugin::milxQtDICOMPlugin(QObject *theParent = 0)
        \brief Default destructor
    */
    milxQtDICOMPlugin(QObject *theParent = 0);
    virtual ~milxQtDICOMPlugin();

    /**
        \fn milxQtDICOMPlugin::name()
        \brief Get the Name of the plugin. [Implement this in your plugin]
    */
    virtual QString name();

    /**
        \fn milxQtDICOMPlugin::hasOpenSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasOpenSupport()
    { return false; }
    /**
        \fn milxQtDICOMPlugin::openFileSupport()
        \brief Get the file support string for opening (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString openFileSupport();
    /**
        \fn milxQtDICOMPlugin::openExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList openExtensions();
    /**
        \fn milxQtDICOMPlugin::hasSaveSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasSaveSupport()
    {   return false;    }
    /**
        \fn milxQtDICOMPlugin::saveFileSupport()
        \brief Get the file support string for saving (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString saveFileSupport();
    /**
        \fn milxQtDICOMPlugin::saveExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList saveExtensions();

    /**
        \fn milxQtDICOMPlugin::hasCollectionSupport()
        \brief Does the plugin support collections (PolyData collection etc.). [Implement this in your plugin]
    */
    inline virtual bool hasCollectionSupport()
    {   return false;    }
    /**
        \fn milxQtDICOMPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
        \brief Pass a collection to internal plugin class. [Implement this in your plugin]
    */
    virtual void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);

    /**
        \fn milxQtDICOMPlugin::open(QString filename)
        \brief Open the file using the plugin. [Implement this in your plugin]
    */
    virtual void open(QString filename);
    /**
        \fn milxQtDICOMPlugin::save(QString filename)
        \brief Save the result as a file using the plugin. [Implement this in your plugin]
    */
    virtual void save(QString filename);

    /**
        \fn milxQtDICOMPlugin::genericResult()
        \brief Get the generic result, which is a milxQtRenderWindow. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtRenderWindow* genericResult();
    /**
        \fn milxQtDICOMPlugin::modelResult()
        \brief Get the model result. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtModel* modelResult();
    /**
        \fn milxQtDICOMPlugin::imageResult()
        \brief Get the image result. The result can then be displayed in milxQtMain etc.[Implement this in your plugin]
    */
    virtual milxQtImage* imageResult();
    /**
        \fn milxQtDICOMPlugin::dockWidget()
        \brief Return the dock widget (if one is provided by plugin). [Implement this in your plugin]
    */
    virtual QDockWidget* dockWidget();

    /**
        \fn milxQtDICOMPlugin::dockDefaultArea()
        \brief Return the default dock widget area (if one is provided by plugin). [Implement this in your plugin]
    */
    inline virtual Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::LeftDockWidgetArea;    }

    /**
        \fn milxQtDICOMPlugin::isPluginWindow(QWidget *window)
        \brief Is the window provided a plugin generated window? In this case a milxQtShapeModel window. [Implement this in your plugin]
    */
    virtual bool isPluginWindow(QWidget *window);

#if (ITK_VERSION_MAJOR > 3) //Review only members
    /*!
      \fn milxQtDICOMPlugin::ExportDICOM_RT(const std::string directoryPath, const std::string rsFileName, const std::string outPath, typename itk::SmartPointer<TImage> &data, std::string &seriesName, const bool reorient = true)
      \brief Opens a DICOM RT from the UID/Series name given. Returns the image volume from the series read by ITK/GDCM.

      You can get the seriesName via the GetDICOMSeriesUIDs() member, which is read from the filenames. The DICOM tags are read and seriesName is replaced with DICOM tag version.

      Returns true if successful. Image is also NOT flipped for VTK.
    */
    template<class TImage>
    static bool ExportDICOM_RT(const std::string directoryPath, const std::string rsFileName, const std::string outPath, typename itk::SmartPointer<TImage> &data, std::string &seriesName, const bool reorient = true);
#endif // (ITK_VERSION_MAJOR > 3)

public slots:
    /**
        \fn milxQtDICOMPlugin::loadExtension()
        \brief Load the extension. [Implement this in your plugin]
    */
    virtual void loadExtension();
    /**
        \fn milxQtDICOMPlugin::update()
        \brief Update the plugin. [Implement this in your plugin]

        This generic call is called after plugin is loaded and is designed to be used to update the plugin
        internals such as manager displays etc.
    */
    virtual void update() {}
    /**
        \fn milxQtDICOMPlugin::preStartTasks()
        \brief Tasks to complete before running or starting the thread. [Implement this]
    */
    virtual void preStartTasks() {}
    /**
        \fn milxQtDICOMPlugin::postStartTasks()
        \brief Tasks to complete after running or starting the thread. [Implement this]
    */
    virtual void postStartTasks() {}

    /**
        \fn milxQtDICOMPlugin::viewTags()
        \brief View DICOM tags of image series in the manager.
    */
    void viewTags();
    /**
    \fn milxQtDICOMPlugin::openSeries()
    \brief Open DICOM image series.
    */
    void openSeries();
    /**
        \fn milxQtDICOMPlugin::openStructureSet()
        \brief Open DICOM RT image series.
    */
    void openStructureSet();

    /**
        \fn milxQtDICOMPlugin::convert()
        \brief Convert a DICOM series (or a path to a number of series) to Nifti (*.nii.gz) images.

        The function assumes either the path provided in the wizard is a path containing a DICOM series or
        that the path provided contains a number of directories each having a DICOM series.
    */
    void convert();
    void anonymize();

    void showRSFileDialog();
    void showInputFileDialog();
    void showRTInputFileDialog();
    void showOutputFileDialog();
    void showRTOutputFileDialog();
    void showInputFileDialogAnonymize();
    void showOutputFileDialogAnonymize();

    /**
     * \brief Function called when the wizard for anonymization is validated
     * It will verify that everything is in order and set a boolean value to true if yes
     */
    void affectValues();

    /**
        \fn milxQtDICOMPlugin::anonymizeDicomImage()
        \brief Function used to anonymize a single dicom image
        \param input: the input image
        \param subject_output_folder: the output folder for the subject
        \param rel_dir: the relative path from the of the subject and the initial dicom
        \param index_subject: the subject id
        \param index_dicom: the id of the current dicom
        \param isFirst: is it the first dicom anonymize for the subject?
        \return true if successfully anonymized
    */
    bool anonymizeDicomImage(const std::string &input, const QString &subject_output_folder, const QString &rel_dir, unsigned int index_subject, unsigned int index_dicom, bool &isFirst);
    
    /**
        \fn milxQtDICOMPlugin::makeFilename()
        \brief Small helper to make filename from dicom tags
        \param path: output directory
        \param gdcmImageIO: the dicom object
        \param index: the dicom index
        \param index: the dicom index
        \param filename: reference to the filename
    */
    void makeFilename(const QString &path, ImageIOType::Pointer gdcmImageIO, unsigned int index, std::string &filename, unsigned int index_subject);

	/**
	\fn milxQtDICOMPlugin::makeFilename()
	\brief Small helper to remove unwanted characters under windows
	\param std: string to strip
	\param charsToRemove: unwanted characters
	*/
	void removeForbiddenChar(std::string &str, char* charsToRemove);

	void writeLog(QString &filename, std::string &output);
    
    /**
        \fn milxQtDICOMPlugin::getTagValue()
        \brief Retrieve tag value fromt he gdcmImageIO object
        \param gdcmImageIO: the dicom object
        \param tag: the dicom tag
        \param tag_value: the tag value
    */
    void getTagValue(ImageIOType::Pointer gdcmImageIO, const std::string &tag, std::string & tag_value);

protected:
    QPointer<milxQtMain> MainWindow;

    QMenu* menuDICOM; //!< DICOM menu
    //----SSM---- (hierarchical deletion)
    QAction* actionOpenSeries; //!< open series action
    QAction* actionTags; //!< open series action
    QAction* actionConvertStructure; //!< convert RT action
    QAction* actionConvert; //!< convert action
    QAction* actionAnonymize; //!< Anonymize action

    //FAI variables for wizard
    QWizard wizard;
    QWizard wizardRT;
    QWizard wizardAnonymize;
    QString inputDirectoryname;
    QLineEdit *txtInputName;
    QString outputDirectoryname;
    QLineEdit *txtOutputName;
    QString inputRTDirectoryname;
    QLineEdit *txtRTInputName;
    QString outputRTDirectoryname;
    QLineEdit *txtRTOutputName;
    QString rsFilename;
    QLineEdit *txtRSName;
    QString inputAnonymizeDirectoryname;
    QLineEdit *txtInputAnonymizeName;
    QString outputAnonymizeDirectoryname;
    QLineEdit *txtOutputAnonymizeName;
    
    QString outputPrefix;
    QLineEdit *txtOutputPrefix;
    
    int outputInitID;
    QLineEdit *txtOutputInitID;
    
    QCheckBox *checkboxPreserveFolderArc;
    QCheckBox *checkboxPatientName;
    QCheckBox *checkboxPatientID;
    QCheckBox *checkboxSeriesDate;
    QCheckBox *checkboxSeriesTime;
    QCheckBox *checkboxStudyID;
    QCheckBox *checkboxStudyDesc;
    QCheckBox *checkboxSeriesNumber;
    QCheckBox *checkboxSequenceName;
    QCheckBox *checkboxProtocolName;
    QCheckBox *checkboxSeriesDescription;
    
    QCheckBox *anonPatientInfo;
    QCheckBox *anonPhysician;
    QCheckBox *anonOperator;
    QCheckBox *anonScanDate;

    //data
    milxQtImage *image; //main window owner
    QPointer<milxQtManager> manager; //!< Manager widget
    QPointer<QDockWidget> dock; //!< Dock widget
    bool valid;

    void createActions();
    void createMenu();
    void createWizard();
    void createRTWizard();
    void createConnections();
    void createWizardAnonymise();

    //remove empty spaces from contour names
    static inline void trim (std::string &str)
    {
        std::string temp;
        for (unsigned int i = 0; i < str.length(); i++)
            if (str[i] != ' ') temp += str[i];
        str = temp;
    }
    //reset image of nD
    template<class TImage>
    static void resetImage(typename itk::SmartPointer<TImage> imageSlice);
    //merge 2D in 3D image
    template<class TImage, class TSliceType>
    static void mergeImages(typename itk::SmartPointer<TSliceType> tempSlice, typename itk::SmartPointer<TImage> finalImage, int iRequiredSlice);

private:

};

class MILXQT_PLUGIN_EXPORT milxQtDICOMPluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQtDICOMPlugin(theParent);  }
};

#if (ITK_VERSION_MAJOR > 3) //Review only members
template<class TImage>
bool milxQtDICOMPlugin::ExportDICOM_RT(const std::string directoryPath, const std::string rsFileName, const std::string outPath, typename itk::SmartPointer<TImage> &data, std::string &seriesName, const bool reorient)
{
  typedef itk::ImageSeriesReader<TImage> ReaderType;
  typedef itk::GDCMImageIO ImageIOType;

  const std::vector<std::string> filenames = milx::File::GetDICOMSeriesFilenames(directoryPath, seriesName);

  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetImageIO( gdcmIO );
    reader->SetFileNames( filenames );
    reader->AddObserver(itk::ProgressEvent(), milx::ProgressUpdates);
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
      std::cerr << "File Exception caught while reading series!" << std::endl;
      std::cerr << excp << std::endl;
      return false;
    }

  ///Print info
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  std::string series_type_id("0008|103e");
  itk::MetaDataDictionary & dic = gdcmIO->GetMetaDataDictionary();
  itk::MetaDataDictionary::ConstIterator series_type_itr = dic.Find(series_type_id);

  std::string caseId = "0010|0020";
  itk::MetaDataDictionary::ConstIterator case_itr = dic.Find( caseId );

  std::string echoNumber = "0018|0086";
  itk::MetaDataDictionary::ConstIterator echoNumber_itr = dic.Find( echoNumber );

  MetaDataStringType::ConstPointer entryValue= dynamic_cast<const MetaDataStringType *>(series_type_itr->second.GetPointer());
  if(entryValue)
  {
    seriesName = entryValue->GetMetaDataObjectValue();
    std::cout << "Series: " << entryValue->GetMetaDataObjectValue() << std::endl;
  }
  MetaDataStringType::ConstPointer entryValue2= dynamic_cast<const MetaDataStringType *>(case_itr->second.GetPointer());
  if(entryValue2)
    std::cout << "Case: " << entryValue2->GetMetaDataObjectValue() << std::endl;
  MetaDataStringType::ConstPointer entryValue3;
  if(dic.HasKey(echoNumber))
  {
    entryValue3 = dynamic_cast<const MetaDataStringType *>(echoNumber_itr->second.GetPointer());
    std::cout << "Echo Number: " << entryValue3->GetMetaDataObjectValue() << std::endl;
  }

  data = reader->GetOutput();
  qApp->processEvents();

  //Begin J. Dowling Code
  typedef itk::Image< typename TImage::PixelType, 2 >   ImageSliceType;
  typename ImageSliceType::Pointer temp2Dimage = ImageSliceType::New();

  gdcm::Reader RTreader;
  RTreader.SetFileName( rsFileName.c_str() );
  if( !RTreader.Read() )
  {
    std::cout << "Problem reading file: " << rsFileName << std::endl;
    return 0;
  }

  resetImage<TImage>(data);

  //we need to create a temporary 2D slice as well...
  typename TImage::RegionType inputRegion = data->GetLargestPossibleRegion();
  typedef itk::ExtractImageFilter< TImage, ImageSliceType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename TImage::SizeType size = inputRegion.GetSize();
  size[2] = 0;
  typename TImage::IndexType start = inputRegion.GetIndex();
  start[2] =0;
  typename TImage::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );
  filter->SetDirectionCollapseToIdentity(); //22.02.2013
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( data );
  filter->AddObserver(itk::ProgressEvent(), milx::ProgressUpdates);
  filter->Update();
  temp2Dimage = filter->GetOutput();

//  const gdcm::FileMetaInformation &h = RTreader.GetFile().GetHeader();
  const gdcm::DataSet& ds = RTreader.GetFile().GetDataSet();
  std::cout << "Parsing: " << rsFileName << std::endl;

  gdcm::MediaStorage ms;
  ms.SetFromFile( RTreader.GetFile() );
  std::cout << "media storage: " << ms << std::endl;

  // (3006,0020) SQ (Sequence with explicit length #=4)      # 370, 1 StructureSetROISequence
  gdcm::Tag tssroisq(0x3006,0x0020);
  if( !ds.FindDataElement( tssroisq ) )
  {
    std::cout << "Problem locating 0x3006,0x0020 - Is this a valid RT Struct file?" << std::endl;
    return 0;
  }
  gdcm::Tag troicsq(0x3006,0x0039);
  if( !ds.FindDataElement( troicsq ) )
  {
    std::cout << "Problem locating 0x3006,0x0039 - Is this a valid RT Struct file?" << std::endl;
    return 0;
  }

  const gdcm::DataElement &roicsq = ds.GetDataElement( troicsq );

  gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = roicsq.GetValueAsSQ();
  if( !sqi || !sqi->GetNumberOfItems() )
  {
    return 0;
  }
  const gdcm::DataElement &ssroisq = ds.GetDataElement( tssroisq );
  gdcm::SmartPointer<gdcm::SequenceOfItems> ssqi = ssroisq.GetValueAsSQ();
  if( !ssqi || !ssqi->GetNumberOfItems() )
  {
    return 0;
  }

  std::cout << "Number of structures found:" << sqi->GetNumberOfItems() << std::endl;

  //loop through structures
  for(unsigned int pd = 0; pd < sqi->GetNumberOfItems(); ++pd)
  {
      const gdcm::Item & item = sqi->GetItem(pd+1); // Item start at #1
      gdcm::Attribute<0x3006,0x0084> roinumber;
      const gdcm::DataSet& nestedds = item.GetNestedDataSet();
      roinumber.SetFromDataElement( nestedds.GetDataElement( roinumber.GetTag() ) );

      qApp->processEvents();

      // find structure_set_roi_sequence corresponding to roi_contour_sequence (by comparing id numbers)
      unsigned int spd = 0;
      gdcm::Item & sitem = ssqi->GetItem(spd+1);
      gdcm::DataSet& snestedds = sitem.GetNestedDataSet();

      gdcm::Attribute<0x3006,0x0022> sroinumber;

      do
      {
        sitem = ssqi->GetItem(spd+1);
        snestedds = sitem.GetNestedDataSet();

        sroinumber.SetFromDataElement( snestedds.GetDataElement( sroinumber.GetTag() ) );

        spd++;

      } while ( sroinumber.GetValue()  != roinumber.GetValue() );

      gdcm::Tag stcsq(0x3006,0x0026);
      if( !snestedds.FindDataElement( stcsq ) )
      {
        std::cout<<"Did not find sttsq data el " << stcsq << "   continuing..." << std::endl;
        continue; //return 0;
      }
      const gdcm::DataElement &sde = snestedds.GetDataElement( stcsq );

      //(3006,002a) IS [255\192\96]                              # 10,3 ROI Display Color
      gdcm::Tag troidc(0x3006,0x002a);
      gdcm::Attribute<0x3006,0x002a> color = {};
      if( nestedds.FindDataElement( troidc) )
      {
        const gdcm::DataElement &decolor = nestedds.GetDataElement( troidc );
        color.SetFromDataElement( decolor );
      }
      //(3006,0040) SQ (Sequence with explicit length #=8)      # 4326, 1 ContourSequence
      gdcm::Tag tcsq(0x3006,0x0040);
      if( !nestedds.FindDataElement( tcsq ) )
      {
        continue;
      }
      const gdcm::DataElement& csq = nestedds.GetDataElement( tcsq );

      gdcm::SmartPointer<gdcm::SequenceOfItems> sqi2 = csq.GetValueAsSQ();
      if( !sqi2 || !sqi2->GetNumberOfItems() )
      {
        std::cout << "csq: " << csq << std::endl;
        std::cout << "sqi2: " << *sqi2 << std::endl;
        std::cout<<"Did not find sqi2 or no. items == 0   " <<  sqi2->GetNumberOfItems() << "   continuing..." << std::endl;
        continue;
      }
      unsigned int nitems = sqi2->GetNumberOfItems();
      std::cout << "Structure " << pd << ". Number of regions: " << nitems << std::endl;
      std::string str_currentOrgan(sde.GetByteValue()->GetPointer(), sde.GetByteValue()->GetLength());

      //trim to remove spaces in organ name which can cause problems in scripts eg. "CBCT01__BULK  BONE .nii" .  Might need to have this as parameter?
      trim (str_currentOrgan);
      std::cout << pd << ". Structure name: " << str_currentOrgan << std::endl;

      //now loop through each item for this structure (eg one prostate region on a single slice is an item)
      typename TImage::PointType point;
      typename TImage::IndexType pixelIndex;
      typedef itk::PolygonSpatialObject<2> PolygonType;
      typedef itk::SpatialObjectPoint<2> PolygonPointType;
      PolygonType::PointListType pointList ;
      PolygonPointType p;
      PolygonType::Pointer polygon = PolygonType::New();;
      typedef itk::GroupSpatialObject<2> GroupType;
      typedef itk::SpatialObjectToImageFilter<GroupType, ImageSliceType> SpatialObjectToImageFilterType;
      GroupType::Pointer group = GroupType::New();
      typename SpatialObjectToImageFilterType::Pointer imageFilter =SpatialObjectToImageFilterType::New();
      int iPointsOutsideBoundary = 0;
      int iCurrentSlice = 0;
      for(unsigned int i = 0; i < nitems; ++i)
      {
        const gdcm::Item & item2 = sqi2->GetItem(i+1); // Item start at #1

        const gdcm::DataSet& nestedds2 = item2.GetNestedDataSet();
        // (3006,0050) DS [43.57636\65.52504\-10.0\46.043102\62.564945\-10.0\49.126537\60.714... # 398,48 ContourData
        gdcm::Tag tcontourdata(0x3006,0x0050);
        const gdcm::DataElement & contourdata = nestedds2.GetDataElement( tcontourdata );

        qApp->processEvents();

        //const gdcm::ByteValue *bv = contourdata.GetByteValue();
        gdcm::Attribute<0x3006,0x0050> at;
        at.SetFromDataElement( contourdata );
        const double* pts = at.GetValues();
        unsigned int npts = at.GetNumberOfValues() / 3;

        for(unsigned int j = 0; j < npts * 3; j+=3)
        {
          point[0] = pts[j+0];
          point[1] = pts[j+1];
          point[2] = pts[j+2];

          //transform points to image co-ordinates
          if (!(data->TransformPhysicalPointToIndex( point, pixelIndex )))
          {
            //Are there points outside the image boundary.  This may occur with automatically segmented objects such as benches or external body outlines?
            iPointsOutsideBoundary++;
          }

          p.SetPosition(pixelIndex[0] ,pixelIndex[1],pixelIndex[2]);

          p.SetRed(1);
          p.SetBlue(1);
          p.SetGreen(1);
          pointList.push_back(p);
        }

        // we have the points for a contour in a single slice.  We need to join these up and insert into the slice as polygon.
        iCurrentSlice = pixelIndex[2];

        //Insert Region
        std::cout << "Inserting region with " << pointList.size() << " points into slice: " << iCurrentSlice << std::endl;
        //reset 2D image
        resetImage<ImageSliceType>(temp2Dimage);

        //need to create a 2D slice here, put the polygon on it, and insert it back into the 3D volume...
        group->AddSpatialObject(polygon); //add a new polygon group

        try
        {
           polygon->SetPoints(pointList);  //so copy them to a polygon object
           imageFilter->SetInput(group);
           imageFilter->SetSize(temp2Dimage->GetLargestPossibleRegion().GetSize());
           imageFilter->AddObserver(itk::ProgressEvent(), milx::ProgressUpdates);
           imageFilter->Update();
           temp2Dimage=imageFilter->GetOutput();
         }
        catch( itk::ExceptionObject & err )
        {
          std::cerr << "Problem setting polygon->SetPoints for this region (non-planar)" << std::endl;
          std::cerr << err << std::endl;
        }

        //merge new polygon from temp image into the contour image
        mergeImages<TImage>(temp2Dimage, data, iCurrentSlice);

         //remove the polygon and clean up pointlist
        group->RemoveSpatialObject(polygon);
        pointList.clear();
      }

      if (iPointsOutsideBoundary > 0)
      {
          std::cout <<  " --" << iPointsOutsideBoundary << " contour points detected outside image boundary. Please check the output volume. " ;
          iPointsOutsideBoundary=0;
      }

      //filename
      std::string strNewVolume = outPath + "/";
      strNewVolume += str_currentOrgan;
      strNewVolume += ".nii.gz";

      //write img
      milx::File::WriteImageUsingITK<TImage>(strNewVolume, data);

      resetImage<TImage>(data);  //reset the temporary volume ready for next structure (if any)
  }
  //End J. Dowling Code

  return true;
}
#endif // (ITK_VERSION_MAJOR > 3)

template<class TImage>
void milxQtDICOMPlugin::resetImage(typename itk::SmartPointer<TImage> imageSlice)
{
  itk::ImageRegionIterator<TImage> imageIt(imageSlice, imageSlice->GetLargestPossibleRegion()) ;
  imageIt.GoToBegin();
  while(!imageIt.IsAtEnd())
  {
    imageIt.Set(0);
    ++imageIt;
  }
}

template<class TImage, class TSliceType>
void milxQtDICOMPlugin::mergeImages(typename itk::SmartPointer<TSliceType> tempSlice, typename itk::SmartPointer<TImage> finalImage, int iRequiredSlice )
{
  typename TImage::PixelType pixelValue =0;
  typename TImage::IndexType pixelIndex;
  typename TSliceType::IndexType sliceIndex;
  int iX = finalImage->GetLargestPossibleRegion().GetSize()[0];
  int iY = finalImage->GetLargestPossibleRegion().GetSize()[1];


  if (iRequiredSlice>0)
  {
  for (int i=0;i<iX;i++)
    for (int j=0;j<iY;j++)
    {
      pixelIndex[0] = i;
      pixelIndex[1] = j;
      sliceIndex[0] =i;
      sliceIndex[1] = j;

      pixelValue = tempSlice->GetPixel(sliceIndex);
      pixelIndex[2] = iRequiredSlice;

       //Disable hole filling (if required please uncomment the next line (and comment the following line)).
       //if (pixelValue != 0)  finalImage->SetPixel(pixelIndex, pixelValue  );
      finalImage->SetPixel(pixelIndex, finalImage->GetPixel(pixelIndex) ^ (pixelValue != 0));
      qApp->processEvents();
    }
  }
}

#endif // MILXQTDICOMPLUGIN_H

