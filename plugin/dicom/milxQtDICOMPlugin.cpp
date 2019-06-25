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
#include "milxQtDICOMPlugin.h"

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

milxQtDICOMPlugin::milxQtDICOMPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = true;
    consoleWindow = false;
    extension = true;
    pluginName = "DICOM";
    //~ dataName = "";
    valid = true;

    manager = new milxQtManager(MainWindow);
        manager->hide();
    dock = new QDockWidget(tr("DICOM Manager"), MainWindow);
        dock->setFeatures(QDockWidget::AllDockWidgetFeatures);
        dock->setWidget(manager);
        dock->setObjectName("DICOM Manager");

    createActions();
    createMenu();
    createWizard();
    createRTWizard();
    createConnections();
    createWizardAnonymise();
}

milxQtDICOMPlugin::~milxQtDICOMPlugin()
{
    if(isRunning() && threaded)
        quit();
    cout << "DICOM Plugin Destroyed." << std::endl;
}

QString milxQtDICOMPlugin::name()
{
    return pluginName;
}

QString milxQtDICOMPlugin::openFileSupport()
{
    QString openPythonExt = "";

    return openPythonExt;
}

QStringList milxQtDICOMPlugin::openExtensions()
{
    QStringList exts;

    return exts;
}

QStringList milxQtDICOMPlugin::saveExtensions()
{
    QStringList exts;

    return exts;
}

QString milxQtDICOMPlugin::saveFileSupport()
{
    QString savePythonExt = "";

    return savePythonExt;
}

void milxQtDICOMPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{

}

void milxQtDICOMPlugin::open(QString filename)
{

}

void milxQtDICOMPlugin::save(QString filename)
{

}

milxQtRenderWindow* milxQtDICOMPlugin::genericResult()
{
    return NULL;
} //No image result

milxQtModel* milxQtDICOMPlugin::modelResult()
{
    //~ denoiseModel = new milxQtDICOMModel;

    return NULL;
    //~ return denoiseModel;
} //No image result

milxQtImage* milxQtDICOMPlugin::imageResult()
{
    return NULL;
} //No image result

QDockWidget* milxQtDICOMPlugin::dockWidget()
{
    return dock;
}

bool milxQtDICOMPlugin::isPluginWindow(QWidget *window)
{
    //~ if(pluginWindow(window) == 0)
        return false;
    //~ else
        //~ return true;
}

//~ milxQtDeNoiseModel* milxQtDeNoisePlugin::pluginWindow(QWidget *window)
//~ {
    //~ if(window)
        //~ return qobject_cast<milxQtDeNoiseModel *>(window);
    //~ return 0;
//~ }

void milxQtDICOMPlugin::loadExtension()
{
    //~ if(!MainWindow->isActiveModel())
        //~ return;

    //~ milxQtModel *currentWin = MainWindow->activeModel();
    //~ milxQtDeNoiseModel *denoiseModel = new milxQtDeNoiseModel(MainWindow); //hierarchical deletion
        //~ denoiseModel->setName(currentWin->getName());
        //~ denoiseModel->SetInput(currentWin->GetOutput());
        //~ denoiseModel->generateModel();

    //~ MainWindow->display(denoiseModel);
}

//~ void milxQtDICOMPlugin::run()
//~ {
    //~ QMutexLocker locker(&mutex); //Lock memory

    //~ ///Execute own thread work here

    //~ //exec();
//~ }

//~ void milxQtDICOMPlugin::createConnections()
//~ {
    //~ //QObject::connect(denoiseAct, SIGNAL(triggered(bool)), denoiseModel, SLOT(denoise()));
//~ }

void milxQtDICOMPlugin::viewTags()
{
    bool success = true;
    QString directoryPath;
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
      return;

    ///Get UIDs
    MainWindow->printInfo("Trying to read DICOMs in " + directoryPath);
    qApp->processEvents();
    std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(directoryPath.toStdString());

    ///Check if input path is a directory of series or just a series
    QStringList directories;
    if(UIDs.empty())
    {
        MainWindow->printWarning("Found no series found in input directory. Checking internal directories");
        QDir directory(directoryPath);
        directory.setFilter(QDir::AllDirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);
        directories = directory.entryList();
    }
    else
    {
        MainWindow->printInfo(QString("Total of ") + QString::number(UIDs.size()) + " UID(s) found");
        directories.push_back("");
    }

    emit working(-1);
    size_t count = 0;
    foreach(const QString &dir, directories)
    {
        QString relPath = directoryPath + "/" + dir;
        MainWindow->printInfo("Reading DICOMs in " + dir);

        qApp->processEvents();
        UIDs = milx::File::GetDICOMSeriesUIDs(relPath.toStdString(), true);
        std::vector<std::string> seriesNames(UIDs.size(), "");
        MainWindow->printInfo(QString("Total of ") + QString::number(UIDs.size()) + " UID(s) found in " + dir);

        ///Update case list in manager
        QStringList headingList;
        headingList << "Tag" << "Value"; //!< Case Browser
        int caseTabIndex = manager->newTab("Tags "+dir, headingList);

        for(size_t j = 0; j < UIDs.size(); j ++)
        {
            qApp->processEvents();
            std::string caseID;
            std::vector< std::pair<std::string, std::string> > tags;
            if (!milx::File::GetDICOMTags<floatImageType>(relPath.toStdString(), tags, UIDs[j], caseID))
              continue;

            QStringList entryName;
            entryName << UIDs[j].c_str() << ""; //!< Case Browser

            //cout << "Loading tags into browser." << std::endl;
            QList<QStringList> entryList;
            for(int k = 0; k < tags.size(); k ++)
            {
                QStringList tagList;
                std::string tag = tags[k].first;
                std::string value = tags[k].second;

                tagList << tag.c_str() << value.c_str();

                //manager->addItem(caseTabIndex, tagList, Qt::NoItemFlags);
                entryList.push_back(tagList);
            }
            manager->addTreeItem(caseTabIndex, entryName, entryList, Qt::NoItemFlags);
        }
        count ++;
    }
    MainWindow->printInfo(QString("Processed ") + QString::number(count) + " DICOM folders in " + directoryPath);

    manager->show();
    MainWindow->printInfo("Done.");
    emit done(-1);
}

void milxQtDICOMPlugin::openSeries()
{
  bool success = true;
  QString directoryPath;
  QPointer<QFileDialog> fileOpener = new QFileDialog;
  QSettings settings("Shekhar Chandra", "milxQt");

  ///If filenames list is empty ask for them
  if (directoryPath.isEmpty())
  {
    QString path = settings.value("recentPath").toString();
    directoryPath = fileOpener->getExistingDirectory(NULL, tr("Open DICOM Directory"),
      path,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  }

  if (directoryPath.isEmpty())
    return;

  ///Get UIDs
  std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(directoryPath.toStdString());

  ///Update case list in manager
  QStringList headingList;
  headingList << "Tag" << "Value"; //!< Case Browser
  int caseTabIndex = manager->newTab("Tag Browser", headingList);

  emit working(-1);
  for (size_t j = 0; j < UIDs.size(); j++)
  {
    qApp->processEvents();
    std::string caseID;
    floatImageType::Pointer floatImg;
    std::vector< std::pair<std::string, std::string> > tags;
    if (!milx::File::OpenDICOMSeriesAndTags<floatImageType>(directoryPath.toStdString(), floatImg, tags, UIDs[j], caseID))
      continue;

    QStringList entryName;
    entryName << UIDs[j].c_str() << ""; //!< Case Browser

    //cout << "Loading tags into browser." << std::endl;
    QList<QStringList> entryList;
    for (int k = 0; k < tags.size(); k++)
    {
      QStringList tagList;
      std::string tag = tags[k].first;
      std::string value = tags[k].second;

      tagList << tag.c_str() << value.c_str();

      //manager->addItem(caseTabIndex, tagList, Qt::NoItemFlags);
      entryList.push_back(tagList);
    }
    manager->addTreeItem(caseTabIndex, entryName, entryList, Qt::NoItemFlags);
    qApp->processEvents();

    //Display Image
    QPointer<milxQtImage> resultImg = new milxQtImage;
    resultImg->setName(UIDs[j].c_str());
    resultImg->setData(floatImg);
    resultImg->generateImage();
    MainWindow->display(resultImg);
  }

  manager->show();
  MainWindow->printInfo("Done.");
  emit done(-1);
}

void milxQtDICOMPlugin::openStructureSet()
{
    cout << "Opening DICOM-RT" << std::endl;
    //use wizaed to ask
    wizardRT.restart();

    int ret = wizardRT.exec();

    if(ret == QDialog::Rejected)
    {
        wizardRT.restart();
        return;
    }

    //Ensure the names set correctly
    rsFilename = txtRSName->text();
    inputRTDirectoryname = txtRTInputName->text();
    outputRTDirectoryname = txtRTOutputName->text();
    qApp->processEvents();

    //Open RS file
    emit working(-1);
//    MainWindow->printInfo("Opening RS file ...");
//    floatImageType::Pointer floatImg = milx::File::ReadImageUsingITK<floatImageType>(rsFilename.toStdString());
//    charImageType::Pointer charImg = milx::File::ReadImageUsingITK<charImageType>(rsFilename.toStdString());
    qApp->processEvents();

    ///Get UIDs and filenames
    const std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(inputRTDirectoryname.toStdString());

    if(UIDs.empty())
    {
        MainWindow->printError("No series found in input directory");
        return;
    }

    ///Export Structures
    std::string seriesName = UIDs[0];
    charImageType::Pointer charImg;

#if (ITK_VERSION_MAJOR > 3) //Review only members
    ExportDICOM_RT<charImageType>(inputRTDirectoryname.toStdString(), rsFilename.toStdString(), outputRTDirectoryname.toStdString(), charImg, seriesName);
#endif // (ITK_VERSION_MAJOR > 3)

    MainWindow->printInfo("Done.");
    emit done(-1);
}

void milxQtDICOMPlugin::convert()
{
    cout << "Converting DICOMs" << std::endl;
    //use wizaed to ask
    wizard.restart();

    int ret = wizard.exec();

    if(ret == QDialog::Rejected)
    {
        wizard.restart();
        return;
    }

    //Ensure the names set correctly
    inputDirectoryname = txtInputName->text();
    outputDirectoryname = txtOutputName->text();

    qApp->processEvents();

    ///Open each series and save
    ///See milxQtFile::openImageSeries() for more relevant code
    ///Get UIDs and filenames
    std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(inputDirectoryname.toStdString(), true);

    ///Check if input path is a directory of series or just a series
    QStringList directories;
    if(UIDs.empty())
    {
        MainWindow->printWarning("Found no series found in input directory. Checking internal directories");
        QDir directory(inputDirectoryname);
        directory.setFilter(QDir::AllDirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);
        directories = directory.entryList();
    }
    else
    {
        MainWindow->printInfo(QString("Total of ") + QString::number(UIDs.size()) + " UID(s) found");
        directories.push_back("");
    }

    ///Open each image
    emit working(-1);
    foreach(const QString &dir, directories)
    {
        qApp->processEvents();
        QString relPath = inputDirectoryname + "/" + dir;
        QString relOutputPath = outputDirectoryname + "/" + dir;
        UIDs = milx::File::GetDICOMSeriesUIDs(relPath.toStdString(), true);
        std::vector<std::string> seriesNames(UIDs.size(), "");
        MainWindow->printInfo(QString("Total of ") + QString::number(UIDs.size()) + " UID(s) found in " + dir);

        QDir outputSubjectDirectory(relOutputPath);
        ///Create the subject's output folder
        bool exist = outputSubjectDirectory.exists();
        if (!exist)
        {
          exist = QDir().mkpath(relOutputPath);
          if (!exist)
            continue;
        }

        for(size_t j = 0; j < UIDs.size(); j ++)
        {
            qApp->processEvents();
            seriesNames[j] = UIDs[j];
            const std::vector<std::string> seriesFilenames = milx::File::GetDICOMSeriesFilenames(relPath.toStdString(), seriesNames[j]);

            //Read Header
            size_t dimensions = 3;
            std::string pixelType, componentType;
            if (!milx::File::ReadImageInformation(seriesFilenames.front(), pixelType, componentType, dimensions))
            {
                milx::PrintError("Failed Reading First Image. Check the image type/file. Skipping.");
                continue;
            }
            milx::PrintInfo("UID: " + seriesNames[j]);
            milx::PrintInfo("Pixel Type: " + pixelType);
            milx::PrintInfo("Component Type: " + componentType);
            milx::PrintInfo("Dimensions: " + milx::NumberToString(dimensions));

            //Open image using relevant type
            std::string caseID;
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
                if (!milx::File::OpenDICOMSeriesAndTags<vectorImageType>(relPath.toStdString(), vectorImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                  milx::PrintError("Failed Reading Vector Images. Skipping.");
                  continue;
                }
                vectorImages = true;
            }
            else if (componentType == "unsigned_char" || componentType == "unsigned char")
            {
                milx::PrintInfo("Detected labelled images.");
                if (!milx::File::OpenDICOMSeriesAndTags<charImageType>(relPath.toStdString(), labelledImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                  milx::PrintError("Failed Reading Labelled Images. Skipping.");
                  continue;
                }
                labelledImages = true;
            }
            else if (componentType == "short" || componentType == "int16")
            {
                milx::PrintInfo("Detected short images.");
                if (!milx::File::OpenDICOMSeriesAndTags<shortImageType>(relPath.toStdString(), shortImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                    milx::PrintError("Failed Reading Short Images. Skipping.");
                    continue;
                }
                shortImages = true;
            }
            else if (componentType == "unsigned_short" || componentType == "unsigned short")
            {
                milx::PrintInfo("Detected unsigned short images.");
                if (!milx::File::OpenDICOMSeriesAndTags<ushortImageType>(relPath.toStdString(), ushortImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                    milx::PrintError("Failed Reading Unsigned Short Images. Skipping.");
                    continue;
                }
                ushortImages = true;
            }
            else if (componentType == "int" || componentType == "signed" || componentType == "int32" || componentType == "int64")
            {
                milx::PrintInfo("Detected integer images.");
                if (!milx::File::OpenDICOMSeriesAndTags<intImageType>(relPath.toStdString(), intImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                  milx::PrintError("Failed Reading Integer Images. Skipping.");
                  continue;
                }
                integerImages = true;
            }
            else if (componentType == "unsigned_int" || componentType == "unsigned int" || componentType == "unsigned")
            {
                milx::PrintInfo("Detected unsigned int images.");
                if (!milx::File::OpenDICOMSeriesAndTags<uintImageType>(relPath.toStdString(), uintImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                  milx::PrintError("Failed Reading Unsigned Integer Images. Skipping.");
                  continue;
                }
                uintegerImages = true;
            }
            else
            {
                milx::PrintInfo("Detected floating point images.");
                if (!milx::File::OpenDICOMSeriesAndTags<floatImageType>(relPath.toStdString(), floatImage, tags, seriesNames[j], caseID)) //Error NOT printed inside
                {
                  milx::PrintError("Failed Reading Images. Skipping.");
                  continue;
                }
            }

            //edit the series to ensure valid name and no spaces
            QString nameSimplified = seriesNames[j].c_str();
            nameSimplified = nameSimplified.simplified();
            nameSimplified.replace(" ", ""); //no spaces
            seriesNames[j] = nameSimplified.toStdString();
            QString caseSimplified = caseID.c_str();
            caseSimplified = caseSimplified.simplified();
            caseSimplified.replace(" ", ""); //no spaces
            caseID = caseSimplified.toStdString();

            MainWindow->printInfo(QString("Opened: ") + seriesNames[j].c_str());
            qApp->processEvents();

            //create filename
            std::string filename = relOutputPath.toStdString() + "/" + dir.toStdString() + "_" + caseID + "_" + seriesNames[j] + ".nii.gz";
            if(caseID == "-1")
              filename = relOutputPath.toStdString() + "/" + dir.toStdString() + seriesNames[j] + ".nii.gz";
            MainWindow->printInfo(QString("Saving DICOM as ") + filename.c_str());

            //Save as relevant type
            if(vectorImages)
                milx::File::SaveImage<vectorImageType>(filename, vectorImage);
            else if(labelledImages)
                milx::File::SaveImage<charImageType>(filename, labelledImage);
            else if(shortImages)
                milx::File::SaveImage<shortImageType>(filename, shortImage);
            else if(ushortImages)
                milx::File::SaveImage<ushortImageType>(filename, ushortImage);
            else if(integerImages)
                milx::File::SaveImage<intImageType>(filename, intImage);
            else if(uintegerImages)
                milx::File::SaveImage<uintImageType>(filename, uintImage);
            else
                milx::File::SaveImage<floatImageType>(filename, floatImage);
        }
    }

    qApp->processEvents();
    emit done(-1);
    MainWindow->printInfo("Done");
}

void milxQtDICOMPlugin::anonymize()
{
  cout << "Anonymizing DICOMs" << std::endl;
  typedef itk::GDCMSeriesFileNames      NamesGeneratorType;
  typedef std::vector< std::string >    FileNamesContainer;
  typedef std::vector< std::string >    SeriesIdContainer;
  
  //use wizaed to ask
  wizardAnonymize.restart();

  int ret = wizardAnonymize.exec();

  if(ret == QDialog::Rejected)
  {
      wizardAnonymize.restart();
      return;
  }
  
  if (!valid)
  {
    return;
  }

  qApp->processEvents();
  
  MainWindow->printInfo(QString("DICOM directory ") + inputAnonymizeDirectoryname);
  MainWindow->printInfo(QString("Output directory ") + outputAnonymizeDirectoryname);
  MainWindow->printInfo(QString("Prefix anonymisation ") + outputPrefix);
  MainWindow->printInfo(QString("Starting at ID ") + QString::number(outputInitID));

  
  emit working(-1);
  MainWindow->printInfo(QString("Working..."));
  ///Open each series and save
  unsigned int initial_id = outputInitID;
  QDir subjectsDir(inputAnonymizeDirectoryname);  
  subjectsDir.setFilter(QDir::Dirs | QDir::NoDotAndDotDot);

  ///Loop over all the subjects
  QStringList entries = subjectsDir.entryList();
  for( QStringList::ConstIterator entry=entries.begin(); entry!=entries.end(); ++entry) 
  {
      bool first = true;
      qApp->processEvents();
      QFileInfo finfo(subjectsDir, *entry);
      
//        std::cout << "Entry -> " << (*entry).toStdString() << std::std::endl;
//        std::cout << "Processing -> " << finfo.absoluteFilePath().toStdString() << std::std::endl;
      QDir subdir(finfo.absoluteFilePath());
      QString subjectPath = subdir.absolutePath();
//        std::cout << "Processing -> " << subjectPath.toStdString() << std::std::endl;
      QDir::setCurrent(subjectPath);
      
      ///Retrieve the dicom series
      try
      {
          NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
          nameGenerator->SetUseSeriesDetails( true );
          nameGenerator->SetRecursive(true);
          nameGenerator->AddSeriesRestriction("0008|0021" ); // -> need check if possible split by subject
          nameGenerator->SetDirectory(subjectPath.toStdString().c_str());

          ///If there are dicom series for the subject
          ///Create the subject's output folder
          const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
          if (seriesUID.size() > 0)
          {
              QString outputSubjectPath = outputAnonymizeDirectoryname+ QDir::separator() + outputPrefix + QString::number(initial_id);
              QDir outputSubjectDirectory(outputSubjectPath);
      
              ///Create the subject's output folder
              bool exist = outputSubjectDirectory.exists();
              if (!exist)
              {
                  exist = QDir().mkpath(outputSubjectPath);
                  if (! exist)
                  {
                      continue;
                  }
              }
              
              SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
              SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
              
              ///Loop over the dicom series
              while( seriesItr != seriesEnd )
              {
                  qApp->processEvents();
                  std::cout << seriesItr->c_str() << std::endl;
                  FileNamesContainer temp = nameGenerator->GetFileNames(seriesItr->c_str());
                  
                  std::vector<std::string>::iterator fileIterator;
                  unsigned int dcm_index = 0;
                  
                  ///Loop over the dicom files and anonymize them
                  for (fileIterator = temp.begin(); fileIterator != temp.end(); fileIterator++)
                  {
                      /// Need to retrieve path relative to subject (in case keep folder architecture)
                      QString rel_path = subdir.relativeFilePath(QString::fromStdString(*fileIterator));
                      QFileInfo img_info(rel_path);
                      QString rel_dir = img_info.path();
                      
                      ///Anonymize
                      bool success = anonymizeDicomImage(*fileIterator, outputSubjectPath, rel_dir, initial_id, dcm_index, first);
                      if (success)
                      {
                          dcm_index++; 
                      }
                      else
                      {
                          std::cout << "Failed to anonymise " << *fileIterator << std::endl;
                      }
                      qApp->processEvents();
                  }
                  seriesItr++;
              }
              initial_id++;
          }
      }
      catch (itk::ExceptionObject &ex)
      {
        std::cerr << ex << std::endl;
        return;
      }
  }

  qApp->processEvents();
  emit done(-1);
  MainWindow->printInfo("Done");
}

void milxQtDICOMPlugin::showRSFileDialog()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString exts = openMedImageExts.c_str();
    QString path = settings.value("recentPath").toString();

    if(!txtRSName->text().isEmpty())
        path = txtRSName->text();

    QFileDialog *fileOpener = new QFileDialog;
    rsFilename = fileOpener->getOpenFileName(&wizardRT,
                            tr("Select RS File to Open"),
                            path,
                            tr(exts.toStdString().c_str()) ); //!< \todo Check and validate extensions support at Open in Main class
    txtRSName->setText(rsFilename);
}

void milxQtDICOMPlugin::showInputFileDialog()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    inputDirectoryname = fileOpener->getExistingDirectory(&wizard,
                                                     tr("Select Open Directory"),
                                                     path,
                                                     QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    txtInputName->setText(inputDirectoryname);
}

void milxQtDICOMPlugin::showRTInputFileDialog()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    inputRTDirectoryname = fileOpener->getExistingDirectory(&wizardRT,
                                                          tr("Select Open Directory"),
                                                          path,
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    txtRTInputName->setText(inputRTDirectoryname);
    txtRSName->setText(inputRTDirectoryname);
}

void milxQtDICOMPlugin::showInputFileDialogAnonymize()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    inputAnonymizeDirectoryname = fileOpener->getExistingDirectory(&wizardAnonymize,
                                                     tr("Select Open Directory"),
                                                     path,
                                                     QFileDialog::DontResolveSymlinks);
    txtInputAnonymizeName->setText(inputAnonymizeDirectoryname);
}

void milxQtDICOMPlugin::showOutputFileDialog()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    outputDirectoryname = fileOpener->getExistingDirectory(&wizard,
                                                     tr("Select Save Directory"),
                                                     path,
                                                     QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    txtOutputName->setText(outputDirectoryname);
}

void milxQtDICOMPlugin::showRTOutputFileDialog()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    outputRTDirectoryname = fileOpener->getExistingDirectory(&wizardRT,
                                                           tr("Select Save Directory"),
                                                           path,
                                                           QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    txtRTOutputName->setText(outputRTDirectoryname);
}

void milxQtDICOMPlugin::showOutputFileDialogAnonymize()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    outputAnonymizeDirectoryname = fileOpener->getExistingDirectory(&wizardAnonymize,
                                                     tr("Select Save Directory"),
                                                     path,
                                                     QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    txtOutputAnonymizeName->setText(outputAnonymizeDirectoryname);
}

void milxQtDICOMPlugin::createActions()
{
    actionOpenSeries = new QAction(MainWindow);
    actionOpenSeries->setIcon(QIcon(":/resources/toolbar/open_series.png"));
    actionOpenSeries->setText(QApplication::translate("MainWindow", "Open Series", 0));
    actionOpenSeries->setShortcut(tr("Ctrl+Alt+o"));
    actionTags = new QAction(MainWindow);
    actionTags->setIcon(QIcon(":/resources/toolbar/search.png"));
    actionTags->setText(QApplication::translate("MainWindow", "View Tags", 0));
    actionTags->setShortcut(tr("Ctrl+Alt+t"));
    actionConvertStructure = new QAction(MainWindow);
    actionConvertStructure->setText(QApplication::translate("DICOMPlugin", "Convert RT/Structure Set ...", 0));
    actionConvertStructure->setShortcut(tr("Ctrl+Alt+s"));
    actionConvertStructure->setDisabled(true);
#if (ITK_VERSION_MAJOR > 3) //Review only members
    actionConvertStructure->setDisabled(false);
#endif // (ITK_VERSION_MAJOR > 3)
    actionConvert = new QAction(MainWindow);
    actionConvert->setText(QApplication::translate("DICOMPlugin", "Convert ...", 0));
    actionConvert->setShortcut(tr("Ctrl+Alt+c"));
    actionAnonymize = new QAction(MainWindow);
    actionAnonymize->setText(QApplication::translate("DICOMPlugin", "Anonymize ...", 0));
    actionAnonymize->setShortcut(tr("Ctrl+Alt+a"));
}

void milxQtDICOMPlugin::createMenu()
{
    menuDICOM = new QMenu(MainWindow);
    menuDICOM->setTitle(QApplication::translate("DICOMPlugin", "DICOM", 0));
    menuDICOM->addAction(actionOpenSeries);
    menuDICOM->addAction(actionTags);
    menuDICOM->addAction(actionConvertStructure);
    menuDICOM->addAction(actionConvert);
    menuDICOM->addAction(actionAnonymize);

    menuToAdd.append(menuDICOM);
}

void milxQtDICOMPlugin::createWizard()
{
  //intro
  QWizardPage *introPage = new QWizardPage;
  introPage->setTitle("Introduction");

  QLabel *label1 = new QLabel("This wizard will help you process DICOM series.");
  label1->setWordWrap(true);
  QLabel *label2 = new QLabel("You will be asked to provide a directory where the DICOM are located.");
  label2->setWordWrap(true);
  QLabel *label3 = new QLabel("You will then be asked to provide an output directory for the process data to be stored.");
  label3->setWordWrap(true);
  QVBoxLayout *introLayout = new QVBoxLayout;
  introLayout->addWidget(label1);
  introLayout->addWidget(label2);
  introLayout->addWidget(label3);
  introPage->setLayout(introLayout);

  //ask for input
  QWizardPage *inputPage = new QWizardPage;
  inputPage->setTitle("Input Directory");

  QLabel *label4 = new QLabel("Please provide the directory to be analysed.");
  label4->setWordWrap(true);
  txtInputName = new QLineEdit;
  QPushButton *btnInputName = new QPushButton;
  connect(btnInputName, SIGNAL(clicked()), this, SLOT(showInputFileDialog()));
  btnInputName->setText("Browse...");
  QHBoxLayout *inputNameLayout = new QHBoxLayout;
  inputNameLayout->addWidget(txtInputName);
  inputNameLayout->addWidget(btnInputName);
  QGroupBox *inputGroupBox = new QGroupBox("Series Directory");
  inputGroupBox->setLayout(inputNameLayout);
  QVBoxLayout *inputLayout = new QVBoxLayout;
  inputLayout->addWidget(label4);
  inputLayout->addWidget(inputGroupBox);
  inputPage->setLayout(inputLayout);

  //ask for output
  QWizardPage *outputPage = new QWizardPage;
  outputPage->setTitle("Output Directory");

  QLabel *label5 = new QLabel("Please provide the output directory for the processing.");
  label5->setWordWrap(true);
  txtOutputName = new QLineEdit;
  QPushButton *btnOutputName = new QPushButton;
  btnOutputName->setText("Browse...");
  connect(btnOutputName, SIGNAL(clicked()), this, SLOT(showOutputFileDialog()));
  QHBoxLayout *ouputNameLayout = new QHBoxLayout;
  ouputNameLayout->addWidget(txtOutputName);
  ouputNameLayout->addWidget(btnOutputName);
  QGroupBox *outputGroupBox = new QGroupBox("Output Directory");
  outputGroupBox->setLayout(ouputNameLayout);
  QVBoxLayout *outputLayout = new QVBoxLayout;
  outputLayout->addWidget(label5);
  outputLayout->addWidget(outputGroupBox);
  outputPage->setLayout(outputLayout);

  //dicom wizard class var
  wizard.addPage(introPage);
  wizard.addPage(inputPage);
  wizard.addPage(outputPage);
  wizard.setWindowTitle("DICOM Wizard");
}

void milxQtDICOMPlugin::createRTWizard()
{
    //intro
    QWizardPage *introPage = new QWizardPage;
    introPage->setTitle("Introduction");

    QLabel *label1 = new QLabel("This wizard will help you process DICOM series.");
    label1->setWordWrap(true);
    QLabel *label2 = new QLabel("You will be asked to provide a directory where the DICOM are located.");
    label2->setWordWrap(true);
    QLabel *label3 = new QLabel("You will then be asked to provide an output directory for the process data to be stored.");
    label3->setWordWrap(true);
    QVBoxLayout *introLayout = new QVBoxLayout;
    introLayout->addWidget(label1);
    introLayout->addWidget(label2);
    introLayout->addWidget(label3);
    introPage->setLayout(introLayout);

    //ask for input
    QWizardPage *inputPage = new QWizardPage;
    inputPage->setTitle("Input Directory");

    QLabel *label4 = new QLabel("Please provide the directory to be analysed.");
    label4->setWordWrap(true);
    txtRTInputName = new QLineEdit;
    QPushButton *btnInputName = new QPushButton;
    connect(btnInputName, SIGNAL(clicked()), this, SLOT(showRTInputFileDialog()));
    btnInputName->setText("Browse...");
    QHBoxLayout *inputNameLayout = new QHBoxLayout;
    inputNameLayout->addWidget(txtRTInputName);
    inputNameLayout->addWidget(btnInputName);
    QGroupBox *inputGroupBox = new QGroupBox("Series Directory");
    inputGroupBox->setLayout(inputNameLayout);
    QVBoxLayout *inputLayout = new QVBoxLayout;
    inputLayout->addWidget(label4);
    inputLayout->addWidget(inputGroupBox);
    inputPage->setLayout(inputLayout);

    //ask for output
    QWizardPage *outputPage = new QWizardPage;
    outputPage->setTitle("Output Directory");

    QLabel *label5 = new QLabel("Please provide the output directory for the processing.");
    label5->setWordWrap(true);
    txtRTOutputName = new QLineEdit;
    QPushButton *btnOutputName = new QPushButton;
    btnOutputName->setText("Browse...");
    connect(btnOutputName, SIGNAL(clicked()), this, SLOT(showRTOutputFileDialog()));
    QHBoxLayout *ouputNameLayout = new QHBoxLayout;
    ouputNameLayout->addWidget(txtRTOutputName);
    ouputNameLayout->addWidget(btnOutputName);
    QGroupBox *outputGroupBox = new QGroupBox("Output Directory");
    outputGroupBox->setLayout(ouputNameLayout);
    QVBoxLayout *outputLayout = new QVBoxLayout;
    outputLayout->addWidget(label5);
    outputLayout->addWidget(outputGroupBox);
    outputPage->setLayout(outputLayout);

    //ask for RS file
    QWizardPage *rsPage = new QWizardPage;
    rsPage->setTitle("RS File");

    QLabel *label6 = new QLabel("Please provide the RS file for processing.");
    label6->setWordWrap(true);
    txtRSName = new QLineEdit;
    QPushButton *btnRSName = new QPushButton;
    btnRSName->setText("Browse...");
    connect(btnRSName, SIGNAL(clicked()), this, SLOT(showRSFileDialog()));
    QHBoxLayout *rsNameLayout = new QHBoxLayout;
    rsNameLayout->addWidget(txtRSName);
    rsNameLayout->addWidget(btnRSName);
    QGroupBox *rsGroupBox = new QGroupBox("RS File");
    rsGroupBox->setLayout(rsNameLayout);
    QVBoxLayout *rsLayout = new QVBoxLayout;
    rsLayout->addWidget(label6);
    rsLayout->addWidget(rsGroupBox);
    rsPage->setLayout(rsLayout);

    //dicom RT wizard class var
    wizardRT.addPage(introPage);
    wizardRT.addPage(inputPage);
    wizardRT.addPage(outputPage);
    wizardRT.addPage(rsPage);
    wizardRT.setWindowTitle("DICOM RT Wizard");
}

void milxQtDICOMPlugin::createWizardAnonymise()
{
    ///Intro
    QWizardPage *introPage = new QWizardPage;
    introPage->setTitle("Introduction");
  
    QLabel *label1 = new QLabel("This wizard will help you anonymize DICOM series.");
    label1->setWordWrap(true);
    QLabel *label2 = new QLabel("You will be asked to provide the 'top'' directory in which the folders of the patients are located.");
    label2->setWordWrap(true);
    QLabel *label3 = new QLabel("You will then be asked to provide an output directory.");
    label3->setWordWrap(true);
    QLabel *label3a = new QLabel("Then, you will then be asked to provide a prefix that will be used to replace the patient name eg. 'ANON'.");
    label3a->setWordWrap(true);
    QLabel *label3b = new QLabel("You will finally be asked to provide an initial ID for the first case to be anonymized (this ID will be incremented per subject).");
    label3b->setWordWrap(true);
    QVBoxLayout *introLayout = new QVBoxLayout;
    introLayout->addWidget(label1);
    introLayout->addWidget(label2);
    introLayout->addWidget(label3);
    introLayout->addWidget(label3a);
    introLayout->addWidget(label3b);
    introPage->setLayout(introLayout);
  
    //ask for input
    ///Page for intput directory (existence will be checked at the end)
    QWizardPage *inputPage = new QWizardPage;
    inputPage->setTitle("Input Directory");
  
    QLabel *label4 = new QLabel("Please provide the top directory in which the subjects' folders are located.");
    label4->setWordWrap(true);
    txtInputAnonymizeName = new QLineEdit;
    txtInputAnonymizeName->setText(QDir::current().absolutePath());
    QPushButton *btnInputName = new QPushButton;
    connect(btnInputName, SIGNAL(clicked()), this, SLOT(showInputFileDialogAnonymize()));
    btnInputName->setText("Browse...");
    QHBoxLayout *inputNameLayout = new QHBoxLayout;
    inputNameLayout->addWidget(txtInputAnonymizeName);
    inputNameLayout->addWidget(btnInputName);
    QGroupBox *inputGroupBox = new QGroupBox("Series Directory");
    inputGroupBox->setLayout(inputNameLayout);
    QVBoxLayout *inputLayout = new QVBoxLayout;
    inputLayout->addWidget(label4);
    inputLayout->addWidget(inputGroupBox);
    inputPage->setLayout(inputLayout);
  
    ///Page for output directory (existence will be checked at the end)
    QWizardPage *outputPage = new QWizardPage;
    outputPage->setTitle("Output Directory");
  
    QLabel *label5 = new QLabel("Please provide the output directory.");
    label5->setWordWrap(true);
    txtOutputAnonymizeName = new QLineEdit;
    txtOutputAnonymizeName->setText(QDir::current().absolutePath());
    QPushButton *btnOutputName = new QPushButton;
    btnOutputName->setText("Browse...");
    connect(btnOutputName, SIGNAL(clicked()), this, SLOT(showOutputFileDialogAnonymize()));
    QHBoxLayout *ouputNameLayout = new QHBoxLayout;
    ouputNameLayout->addWidget(txtOutputAnonymizeName);
    ouputNameLayout->addWidget(btnOutputName);
    QGroupBox *outputGroupBox = new QGroupBox("Output Directory");
    outputGroupBox->setLayout(ouputNameLayout);
    QVBoxLayout *outputLayout = new QVBoxLayout;
    outputLayout->addWidget(label5);
    outputLayout->addWidget(outputGroupBox);
    outputPage->setLayout(outputLayout);
    
    ///Page for choosing prefix (eg. anon)
    QWizardPage *prefixPage = new QWizardPage;
    prefixPage->setTitle("Prefix anonymization");
    
    QLabel *label6 = new QLabel("Please provide the pefix for naming anonymized subject. If left emply 'Anon' will be used");
    label6->setWordWrap(true);
    txtOutputPrefix = new QLineEdit;
    txtOutputPrefix->setText("Anon");
    QHBoxLayout *prfxNameLayout = new QHBoxLayout;
    prfxNameLayout->addWidget(txtOutputPrefix);
    QGroupBox *prefixGroupBox = new QGroupBox("Output prefix");
    prefixGroupBox->setLayout(prfxNameLayout);
    QVBoxLayout *prefixLayout = new QVBoxLayout;
    prefixLayout->addWidget(label6);
    prefixLayout->addWidget(prefixGroupBox);
    prefixPage->setLayout(prefixLayout);
    
    ///Page for choosing initial ID (will then increment per subject)
    QWizardPage *outputIDPage = new QWizardPage;
    outputIDPage->setTitle("Output Directory");
  
    QLabel *label7 = new QLabel("Please provide the initial ID for anonymisation (integer). If left empty, '0' will be used");
    label7->setWordWrap(true);
    txtOutputInitID = new QLineEdit;
    txtOutputInitID->setText("0");
    QHBoxLayout *InitIDNameLayout = new QHBoxLayout;
    InitIDNameLayout->addWidget(txtOutputInitID);
    //InitIDNameLayout->addWidget(btnInitIDName);
    QGroupBox *InitIDGroupBox = new QGroupBox("Initial ID");
    InitIDGroupBox->setLayout(InitIDNameLayout);
    QVBoxLayout *InitIDLayout = new QVBoxLayout;
    InitIDLayout->addWidget(label7);
    InitIDLayout->addWidget(InitIDGroupBox);
    outputIDPage->setLayout(InitIDLayout);
    
    ///Page for choosing options for anonymization / filenames / folder arch
    QWizardPage *optionPage = new QWizardPage;
    optionPage->setTitle("Options");
  
    QHBoxLayout *splitLayout = new QHBoxLayout;
    
    QVBoxLayout *optionLayout = new QVBoxLayout;
    QLabel *label8 = new QLabel("Output options:");
    label8->setWordWrap(true);
    
    checkboxPreserveFolderArc = new QCheckBox("Preserve folder architecture", &wizardAnonymize);
    checkboxPreserveFolderArc->setChecked(true);
    
    QVBoxLayout *anonOptionLayout = new QVBoxLayout;
    QLabel *label10 = new QLabel("Anonymization options:");
    label10->setWordWrap(true);
    
    anonPatientInfo = new QCheckBox("Anonymize Patient info:", &wizardAnonymize);
    anonPatientInfo->setChecked(true);
    anonPhysician = new QCheckBox("Anonymize Physician(s) info:", &wizardAnonymize);
    anonOperator = new QCheckBox("Anonymize Operator info:", &wizardAnonymize);
    optionLayout->addWidget(label10);
    optionLayout->addWidget(checkboxPreserveFolderArc);
    optionLayout->addWidget(anonPatientInfo);
    optionLayout->addWidget(anonPhysician);
    optionLayout->addWidget(anonOperator);
    optionLayout->setAlignment(Qt::AlignTop);
  
    QLabel *label9 = new QLabel("Items for default filename:");
    label8->setWordWrap(true);
   
    checkboxPatientName = new QCheckBox("Patient anonymized name:", &wizardAnonymize);
    checkboxPatientName->setChecked(true);
    checkboxPatientID = new QCheckBox("PatientID:", &wizardAnonymize);
    checkboxPatientID->setChecked(true);
    checkboxSeriesDate = new QCheckBox("SeriesDate:", &wizardAnonymize);
    checkboxSeriesDate->setChecked(true);
    checkboxSeriesTime = new QCheckBox("SeriesTime:", &wizardAnonymize);
    checkboxSeriesTime->setChecked(true);
    checkboxStudyID = new QCheckBox("StudyID:", &wizardAnonymize);
    checkboxStudyDesc = new QCheckBox("StudyDesc:", &wizardAnonymize);
    checkboxStudyDesc->setChecked(true);
    checkboxSeriesNumber = new QCheckBox("SeriesNumber:", &wizardAnonymize);
    checkboxSeriesNumber->setChecked(true);
    checkboxSeriesNumber->setEnabled(false);
    checkboxSequenceName = new QCheckBox("SequenceName:", &wizardAnonymize);
    checkboxSequenceName->setChecked(true);
    checkboxProtocolName = new QCheckBox("ProtocolName:", &wizardAnonymize);
    checkboxProtocolName->setChecked(true);
    checkboxSeriesDescription = new QCheckBox("SeriesDescription:", &wizardAnonymize);
   
    anonOptionLayout->addWidget(label8);
    anonOptionLayout->addWidget(label9);
    anonOptionLayout->addWidget(checkboxPatientName);
    anonOptionLayout->addWidget(checkboxPatientID);
    anonOptionLayout->addWidget(checkboxSeriesDate);
    anonOptionLayout->addWidget(checkboxSeriesTime);
    anonOptionLayout->addWidget(checkboxStudyID);
    anonOptionLayout->addWidget(checkboxStudyDesc);
    anonOptionLayout->addWidget(checkboxSeriesNumber);
    anonOptionLayout->addWidget(checkboxSequenceName);
    anonOptionLayout->addWidget(checkboxProtocolName);
    anonOptionLayout->addWidget(checkboxSeriesDescription);
    anonOptionLayout->setAlignment(Qt::AlignTop);
    
    splitLayout->addLayout(optionLayout);
    splitLayout->addLayout(anonOptionLayout);
    optionPage->setLayout(splitLayout);
  
    //wizard class var
    wizardAnonymize.addPage(introPage);
    wizardAnonymize.addPage(inputPage);
    wizardAnonymize.addPage(outputPage);
    wizardAnonymize.addPage(prefixPage);
    wizardAnonymize.addPage(outputIDPage);
    wizardAnonymize.addPage(optionPage);
    
    wizardAnonymize.setWindowTitle("DICOM Anonymisation Wizard");
    
    QAbstractButton *final = wizardAnonymize.button(QWizard::FinishButton);
    connect(final, SIGNAL(clicked()), this, SLOT(affectValues()));
}

void milxQtDICOMPlugin::createConnections()
{
//    connect(MainWindow, SIGNAL(windowActivated(QWidget*)), this, SLOT(updateManager(QWidget*)));
    connect(actionOpenSeries, SIGNAL(activated()), this, SLOT(openSeries()));
    connect(actionTags, SIGNAL(activated()), this, SLOT(viewTags()));
    connect(actionConvertStructure, SIGNAL(activated()), this, SLOT(openStructureSet()));
    connect(actionConvert, SIGNAL(activated()), this, SLOT(convert()));
    connect(actionAnonymize, SIGNAL(activated()), this, SLOT(anonymize()));

    connect(this, SIGNAL(working(int)), MainWindow, SLOT(working(int)));
    connect(this, SIGNAL(done(int)), MainWindow, SLOT(done(int)));
}

bool milxQtDICOMPlugin::anonymizeDicomImage(const std::string &input, const QString &subject_output_folder, const QString &rel_dir, unsigned int index_subject, unsigned int index_dicom, bool &isFirst)
{
	QString logName= outputAnonymizeDirectoryname + QDir::separator() + "anonymization.log";
	std::string log = "Anonymizing image: " + input;
	writeLog(logName, log);
	MainWindow->printInfo(QString("Anonymizing: ") + input.c_str());

    //TODO: This should be retrieved directly from dicom.
	// And specific image type created then.
    typedef signed short shortPixelType;
    typedef itk::Image<shortPixelType, milx::imgDimension> shortImageType;

    ///Create the anonymization string -  This is the folder of the subject
    std::ostringstream index;
    index << index_subject;
    std::string anonymization_value = outputPrefix.toStdString() + index.str();
    
    ///Create 2 readers --> need two in case of T2 maps
	///TODO: Again this should be retrieved from the DICOM
    typedef itk::ImageFileReader< shortImageType >  ReaderType;
    typedef itk::ImageFileReader< rgbImageType >    rgbReaderType;
    
    ///First the image need to be read with a random type to able to read the dicom tags
	///and checkthe image type for saving later one
    ReaderType::Pointer reader = ReaderType::New();
    ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
    reader->SetFileName(input.c_str());
    reader->SetImageIO( gdcmImageIO );
    try
    {
        reader->Update();  
    }
    catch (itk::ExceptionObject &ex)
    {
        std::cerr << ex << std::endl;
        return false;
    }
    shortImageType::Pointer floatImg = reader->GetOutput();
  
    ///Write the mapping of the subject's real name to his identifier
    if (isFirst)
    {
        QString fileMappingName = outputAnonymizeDirectoryname + QDir::separator() + "subjects_mapping.txt";
        QFile file_mapping(fileMappingName);
        if (file_mapping.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Append))
        {
            isFirst = false;
            ///Here retrieve real name
            std::string name_tag("0010|0010");
            std::string name_value;
            getTagValue(gdcmImageIO, name_tag, name_value);
            
            ///And write
            QTextStream out(&file_mapping);
                out.setCodec("UTF-8");
                out << QString::fromStdString(name_value) << "," << QString::fromStdString(anonymization_value) << "\n";
                file_mapping.close();    
        } 
    }
    
    ///If it is a T2 map need change image type / and re-read (unfortunate)
	///This should be read from the function rather than dicom tags but works
	///as well
    std::string spl_per_px("0028|0002");
    std::string spl_per_pxl_value;
    getTagValue(gdcmImageIO, spl_per_px, spl_per_pxl_value);
    bool isRGB = false;
    if (spl_per_pxl_value == "3")
    {
		isRGB = true;
    }

    ///List of DICOM tags to strip
    std::vector<std::string>  dicomTags;
    if (anonPatientInfo->isChecked())
    {
        dicomTags.push_back("0010|0010");
        dicomTags.push_back("0010|0020");
        dicomTags.push_back("0010|1005");
        dicomTags.push_back("0010|1040");
        dicomTags.push_back("0010|2154");
    }
  
    if (anonPhysician->isChecked())
    {
        dicomTags.push_back("0008|0090");
        dicomTags.push_back("0008|0092");
        dicomTags.push_back("0008|0094");
        dicomTags.push_back("0008|1048");
        dicomTags.push_back("0008|1050");
        dicomTags.push_back("0008|1060");
        dicomTags.push_back("0032|1032");
        dicomTags.push_back("0040|0006");
        dicomTags.push_back("4008|0114");
    }
  
    if (anonOperator->isChecked())
    {
        dicomTags.push_back("0008|1070");
    }
    
    ///Create output directory
    QString outputSequence = subject_output_folder + QDir::separator();
    if (checkboxPreserveFolderArc->isChecked())
    {
        outputSequence = outputSequence + rel_dir;
    }
    else
    {
        std::string tag1("0020|0011");
        std::string append1;
        getTagValue(gdcmImageIO, tag1, append1);
		removeForbiddenChar(append1, "\\/:*?\"<>|");
        outputSequence = outputSequence + QString::fromStdString(append1) + "_";
    
        std::string tag2("0018|0024");
        std::string append2;
        getTagValue(gdcmImageIO, tag2, append2);
		removeForbiddenChar(append2, "\\/:*?\"<>|");
        outputSequence = outputSequence + QString::fromStdString(append2);
    }
  
    QDir outputSequenceDirectory(outputSequence);
    bool exist = outputSequenceDirectory.exists();
    if (!exist)
    {
        exist = QDir().mkpath(outputSequence);
        if (! exist)
        {
			MainWindow->printError("Failed to create directory: " + outputSequence);
            std::cout << "Create fail: " << anonymization_value << std::endl;
            return false;
        }
    }

	///Write some debug on the type of images
	log = "Component: [" + gdcmImageIO->GetComponentTypeAsString(gdcmImageIO->GetInternalComponentType()) + "]";
	log += ", Pixel type: [" + gdcmImageIO->GetPixelTypeAsString(gdcmImageIO->GetPixelType()) + "]\n";
	writeLog(logName, log);
    
    ///Change dicom header and re-write file
    if (!isRGB) // If normal image
	{
		log = "Changing header";
		writeLog(logName, log);
        DictionaryType & dictionary = floatImg->GetMetaDataDictionary();  
        std::vector<std::string>::iterator dicomTagIterator;
        for (dicomTagIterator = dicomTags.begin(); dicomTagIterator != dicomTags.end(); dicomTagIterator++)
        {
            itk::EncapsulateMetaData<std::string>(dictionary, *dicomTagIterator, anonymization_value);  
        }
      
        //make output filename
        std::string filename;
        makeFilename(outputSequence, gdcmImageIO, index_dicom, filename, index_subject);

        
        typedef itk::ImageFileWriter< shortImageType >  Writer1Type;
        gdcmImageIO->KeepOriginalUIDOn();
        Writer1Type::Pointer writer1 = Writer1Type::New();
        writer1->SetInput(floatImg);
        writer1->SetFileName(filename.c_str());
        writer1->SetImageIO( gdcmImageIO);
		log = "Writing image as " + filename;
		writeLog(logName, log);
		MainWindow->printInfo(QString("Writing image as: ") + filename.c_str());
        try
        {
            writer1->Update();
        }
        catch (itk::ExceptionObject &ex)
        {
			MainWindow->printError(ex.GetDescription());
            std::cerr << ex << std::endl;
            return false;
        } 
    }
    else // If RGB image
    {
		rgbReaderType::Pointer readerrgb = rgbReaderType::New();

		gdcmImageIO = NULL;
		gdcmImageIO = ImageIOType::New();
		readerrgb->SetFileName(input.c_str());
		readerrgb->SetImageIO(gdcmImageIO);
		try
		{
			readerrgb->Update();
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cerr << ex << std::endl;
			return false;
		}
		rgbImageType::Pointer rgbImg = readerrgb->GetOutput();

		ImageIOType::Pointer gdcmImageIO2 = ImageIOType::New();

		log = "Changing header";
		writeLog(logName, log);

        DictionaryType & dictionary = rgbImg->GetMetaDataDictionary();  
		std::vector<std::string>::iterator dicomTagIterator;
        for (dicomTagIterator = dicomTags.begin(); dicomTagIterator != dicomTags.end(); dicomTagIterator++)
        {
            itk::EncapsulateMetaData<std::string>(dictionary, *dicomTagIterator, anonymization_value);  
        }
      
        //make output filename
        std::string filename;
        makeFilename(outputSequence, gdcmImageIO, index_dicom, filename, index_subject);
        
        typedef itk::ImageFileWriter< rgbImageType >  Writer1Type;
		gdcmImageIO2->SetMetaDataDictionary(dictionary);
		gdcmImageIO2->SetUseCompression(gdcmImageIO->GetUseCompression());
		gdcmImageIO2->SetIORegion(gdcmImageIO->GetIORegion());
		gdcmImageIO2->SetUIDPrefix(gdcmImageIO->GetUIDPrefix());
		gdcmImageIO2->KeepOriginalUIDOn();

        Writer1Type::Pointer writer1 = Writer1Type::New();
        writer1->SetInput(rgbImg);
		//writer1->SetMetaDataDictionary(dictionary);
        writer1->SetFileName(filename.c_str());
		writer1->SetImageIO(gdcmImageIO2);
		writer1->UseInputMetaDataDictionaryOff();
		log = "Writing image as " + filename;
		writeLog(logName, log);
		MainWindow->printInfo(QString("Writing image as: ") + filename.c_str());
        try
        {
            writer1->Update();
        }
        catch (itk::ExceptionObject &ex)
        {
			MainWindow->printError(ex.GetDescription());
            std::cerr << ex << std::endl;
            return false;
        }
    }
	log = "Success\n";
	writeLog(logName, log);
    return true;
}

void milxQtDICOMPlugin::makeFilename(const QString &path, ImageIOType::Pointer gdcmImageIO, unsigned int index, std::string &filename, unsigned int index_subject)
{
	std::string tag_idx("0020|0013");
	std::string append_idx;
	getTagValue(gdcmImageIO, tag_idx, append_idx);

	std::ostringstream index_dcm;
	index_dcm << std::setfill('0') << std::setw(4) << append_idx;


	QString str_separator(QDir::separator());
	filename = path.toStdString() + str_separator.toStdString();

	if (checkboxPatientName->isChecked())
	{
		std::string append = outputPrefix.toStdString();
		removeForbiddenChar(append, "\\/:*?\"<>|");
		filename = filename + append + "_";
	}

	if (checkboxPatientID->isChecked())
	{
		std::ostringstream SubID;
		SubID << index_subject;
		filename = filename + SubID.str() + "_";
	}

	if (checkboxSeriesDate->isChecked())
	{
		std::string tag("0008|0021");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxSeriesTime->isChecked())
	{
		std::string tag("0008|0031");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxStudyID->isChecked())
	{
		std::string tag("0020|0010");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		removeForbiddenChar(append, "\\/:*?\"<>|");
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxStudyDesc->isChecked())
	{
		std::string tag("0008|1030");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		removeForbiddenChar(append, "\\/:*?\"<>|");
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxSeriesNumber->isChecked())
	{
		std::string tag("0020|0011");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxSequenceName->isChecked())
	{
		std::string tag("0018|0024");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		removeForbiddenChar(append, "\\/:*?\"<>|");
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxProtocolName->isChecked())
	{
		std::string tag("0018|1030");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		removeForbiddenChar(append, "\\/:*?\"<>|");
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	if (checkboxSeriesDescription->isChecked())
	{
		std::string tag("0008|103e");
		std::string append;
		getTagValue(gdcmImageIO, tag, append);
		removeForbiddenChar(append, "\\/:*?\"<>|");
		if (append != "")
		{
			filename = filename + append + "_";
		}
	}

	filename = filename + index_dcm.str() + ".IMA";
}

void milxQtDICOMPlugin::removeForbiddenChar(std::string &str, char* charsToRemove)
{
	for (unsigned int i = 0; i < strlen(charsToRemove); ++i)
	{
		str.erase(remove(str.begin(), str.end(), charsToRemove[i]), str.end());
	}
}

void milxQtDICOMPlugin::getTagValue(ImageIOType::Pointer gdcmImageIO, const std::string &tag, std::string & tag_value)
{
    itk::MetaDataDictionary & dic = gdcmImageIO->GetMetaDataDictionary();
    itk::MetaDataDictionary::ConstIterator tag_iter = dic.Find(tag);
  
    if (tag_iter == dic.End()) //If we don't find the tag
    {
        return;    
    }
    
    MetaDataStringType::ConstPointer entryValue= dynamic_cast<const MetaDataStringType *>(tag_iter->second.GetPointer());
    if(entryValue)
    {
        tag_value = entryValue->GetMetaDataObjectValue();
        QString tag_simplified = QString::fromStdString(tag_value);
        tag_simplified = tag_simplified.simplified();
        tag_value = tag_simplified.toStdString();
        std::replace(tag_value.begin(), tag_value.end(),' ','_');
    }
}

void milxQtDICOMPlugin::writeLog(QString &filename, std::string &output)
{
	QFile file_log(filename);
	if (file_log.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Append))
	{
		QTextStream out_log(&file_log);
		out_log.setCodec("UTF-8");
		out_log << QString::fromStdString(output) << "\n";
		file_log.close();
	}
}

void milxQtDICOMPlugin::affectValues()
{
    outputAnonymizeDirectoryname = txtOutputAnonymizeName->text();
    inputAnonymizeDirectoryname = txtInputAnonymizeName->text();
    QDir output(outputAnonymizeDirectoryname);
    QDir input(inputAnonymizeDirectoryname);
    if (output.exists() && input.exists())
    {
        valid = true;
    }
    else
    {
        valid = false;
    }
    
    if (!txtOutputPrefix->text().isEmpty())
    {
        outputPrefix = txtOutputPrefix->text();    
    }
    else
    {
        outputPrefix = "Anon";
    }
    
    QString value = txtOutputInitID->text();
    if (!value.isEmpty())
        outputInitID = value.toInt();
    else
        outputInitID = 1;
}

//Q_EXPORT_PLUGIN2(DICOMPlugin, milxQtDICOMPluginFactory);

