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

#include <qplugin.h>

milxQtDICOMPlugin::milxQtDICOMPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = false;
    consoleWindow = false;
    extension = true;
    pluginName = "DICOM";
    //~ dataName = "";
    valid = true;

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
    cout << "DICOM Plugin Destroyed." << endl;
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
    return NULL;
} //No Dock result

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

void milxQtDICOMPlugin::openStructureSet()
{
    cout << "Opening DICOM-RT" << endl;
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
    ExportDICOM_RT<charImageType>(inputRTDirectoryname.toStdString(), rsFilename.toStdString(), outputRTDirectoryname.toStdString(), charImg, seriesName);

    MainWindow->printInfo("Done.");
    emit done(-1);
}

void milxQtDICOMPlugin::convert()
{
    cout << "Converting DICOMs" << endl;
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

//    if(MainWindow->getNumberOfTabs() > 1)
//        MainWindow->newTab();
    ///Open each series and save
    ///See milxQtFile::openImageSeries() for more relevant code

    ///Get UIDs and filenames
    const std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(inputDirectoryname.toStdString());
    std::vector<std::string> seriesNames(UIDs.size(), "");
    MainWindow->printInfo(QString("Total of ") + QString::number(UIDs.size()) + " UID(s) found");

    if(UIDs.empty())
    {
        MainWindow->printError("No series found in input directory");
        return;
    }

    ///Open each image
    emit working(-1);
    for(size_t j = 0; j < UIDs.size(); j ++)
    {
        qApp->processEvents();
//        const std::vector<std::string> filenames = milx::File::GetDICOMSeriesFilenames(directoryPath, UIDs[j]);

        //Open
        seriesNames[j] = UIDs[j];
        floatImageType::Pointer floatImg;
        MainWindow->printInfo(QString("Opening: ") + seriesNames[j].c_str());
        milx::File::OpenDICOMSeries<floatImageType>(inputDirectoryname.toStdString(), floatImg, seriesNames[j]);
        qApp->processEvents();

        //create filename
        std::string filename = outputDirectoryname.toStdString() + "/" + seriesNames[j] + ".nii.gz";
        MainWindow->printInfo(QString("Saving DICOM as ") + filename.c_str());

        //Save
        milx::File::SaveImage<floatImageType>(filename, floatImg);
    }

    qApp->processEvents();
    emit done(-1);
    MainWindow->printInfo("Done");
}

void milxQtDICOMPlugin::anonymize()
{
  cout << "Anonymizing DICOMs" << endl;
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
      
//        std::cout << "Entry -> " << (*entry).toStdString() << std::endl;
//        std::cout << "Processing -> " << finfo.absoluteFilePath().toStdString() << std::endl;
      QDir subdir(finfo.absoluteFilePath());
      QString subjectPath = subdir.absolutePath();
//        std::cout << "Processing -> " << subjectPath.toStdString() << std::endl;
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
    actionOpenSeries->setText(QApplication::translate("MainWindow", "Open Series", 0, QApplication::UnicodeUTF8));
    actionOpenSeries->setShortcut(tr("Ctrl+Alt+o"));
    actionConvertStructure = new QAction(MainWindow);
    actionConvertStructure->setText(QApplication::translate("DICOMPlugin", "Convert RT/Structure Set ...", 0, QApplication::UnicodeUTF8));
    actionConvertStructure->setShortcut(tr("Ctrl+Alt+s"));
    actionConvert = new QAction(MainWindow);
    actionConvert->setText(QApplication::translate("DICOMPlugin", "Convert ...", 0, QApplication::UnicodeUTF8));
    actionConvert->setShortcut(tr("Ctrl+Alt+c"));
    actionAnonymize = new QAction(MainWindow);
    actionAnonymize->setText(QApplication::translate("DICOMPlugin", "Anonymize ...", 0, QApplication::UnicodeUTF8));
    actionAnonymize->setShortcut(tr("Ctrl+Alt+a"));
}

void milxQtDICOMPlugin::createMenu()
{
    menuDICOM = new QMenu(MainWindow);
    menuDICOM->setTitle(QApplication::translate("DICOMPlugin", "DICOM", 0, QApplication::UnicodeUTF8));
    menuDICOM->addAction(actionOpenSeries);
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
    connect(actionOpenSeries, SIGNAL(activated()), MainWindow, SLOT(openSeries()));
    connect(actionConvertStructure, SIGNAL(activated()), this, SLOT(openStructureSet()));
    connect(actionConvert, SIGNAL(activated()), this, SLOT(convert()));
    connect(actionAnonymize, SIGNAL(activated()), this, SLOT(anonymize()));

    connect(this, SIGNAL(working(int)), MainWindow, SLOT(working(int)));
    connect(this, SIGNAL(done(int)), MainWindow, SLOT(done(int)));
}

bool milxQtDICOMPlugin::anonymizeDicomImage(const std::string &input, const QString &subject_output_folder, const QString &rel_dir, unsigned int index_subject, unsigned int index_dicom, bool &isFirst)
{
    //May be useful if something decides to fail at some point
    //MainWindow->printInfo(QString("Anonymizing: ") + input.c_str());
    typedef signed short shortPixelType;
    typedef itk::Image<shortPixelType, milx::imgDimension> shortImageType;
    ///Create the anonymization string
    std::ostringstream index;
    index << index_subject;
    std::string anonymization_value = outputPrefix.toStdString() + index.str();
    
    ///Create 2 readers --> need two in case of T2 maps
    typedef itk::ImageFileReader< shortImageType >  ReaderType;
    typedef itk::ImageFileReader< rgbImageType >    rgbReaderType;
    
    ///Read as normal image first
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
    std::string spl_per_px("0028|0002");
    std::string spl_per_pxl_value;
    getTagValue(gdcmImageIO, spl_per_px, spl_per_pxl_value);
    bool isRGB = false;
    rgbImageType::Pointer rgbImg;
    if (spl_per_pxl_value == "3")
    {
        gdcmImageIO = NULL;
        gdcmImageIO = ImageIOType::New();
        rgbReaderType::Pointer readerrgb = rgbReaderType::New();
        readerrgb->SetFileName(input.c_str());
        readerrgb->SetImageIO( gdcmImageIO );
        try
        {
            readerrgb->Update();  
        }
        catch (itk::ExceptionObject &ex)
        {
            std::cerr << ex << std::endl;
            return false;
        }
        isRGB = true;
        rgbImg = readerrgb->GetOutput();
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
        outputSequence = outputSequence + QString::fromStdString(append1) + "_";
    
        std::string tag2("0018|0024");
        std::string append2;
        getTagValue(gdcmImageIO, tag2, append2);
        outputSequence = outputSequence + QString::fromStdString(append2);
    }
  
    QDir outputSequenceDirectory(outputSequence);
    bool exist = outputSequenceDirectory.exists();
    if (!exist)
    {
        exist = QDir().mkpath(outputSequence);
        if (! exist)
        {
            std::cout << "Create fail: " << anonymization_value << std::endl;
            return false;
        }
    }
    
    ///Change dicom header and re-write file
    if (!isRGB) // If normal image
    {
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
        try
        {
            writer1->Update();
        }
        catch (itk::ExceptionObject &ex)
        {
            std::cerr << ex << std::endl;
            return false;
        } 
    }
    else // If RGB image
    {
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
        gdcmImageIO->KeepOriginalUIDOn();
        Writer1Type::Pointer writer1 = Writer1Type::New();
        writer1->SetInput(rgbImg);
        writer1->SetFileName(filename.c_str());
        writer1->SetImageIO( gdcmImageIO);
        try
        {
            writer1->Update();
        }
        catch (itk::ExceptionObject &ex)
        {
            std::cerr << ex << std::endl;
            return false;
        }     
    }
    
    return true;
}

void milxQtDICOMPlugin::makeFilename(const QString &path, ImageIOType::Pointer gdcmImageIO, unsigned int index, std::string &filename, unsigned int index_subject)
{
    std::ostringstream index_dcm;
    index_dcm << std::setfill('0') << std::setw(4) << index;
  
    QString str_separator(QDir::separator());
    filename = path.toStdString() + str_separator.toStdString(); 
  
    if (checkboxPatientName->isChecked())
    {
        filename = filename + outputPrefix.toStdString() + "_";
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
        if (append != "")
        {
            filename = filename + append + "_";
        }
    }
  
    filename = filename + index_dcm.str() + ".IMA";
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

Q_EXPORT_PLUGIN2(DICOMPlugin, milxQtDICOMPluginFactory);
