#include "milxQtDICOMPlugin.h"

#include <qplugin.h>

//#include "milxFile.h"
//#include "milxQtImage.h"

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

    createActions();
    createMenu();
    createWizard();
    createConnections();
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

    qApp->processEvents();
/*
//    if(MainWindow->getNumberOfTabs() > 1)
//        MainWindow->newTab();
    ///Open each series and save
    ///See milxQtFile::openImageSeries() for more relevant code

    ///Get UIDs
    std::vector<std::string> UIDs = milx::File::GetDICOMSeriesUIDs(inputDirectoryname);

    ///Open each image
    for(size_t j = 0; j < UIDs.size(); j ++)
    {
        floatImageType::Pointer floatImg;
        milx::File::OpenDICOMSeries<floatImageType>(inputDirectoryname, floatImg, UIDs[j]);
    }
*/
}

void milxQtDICOMPlugin::anonymize()
{
    cout << "Anonymizing DICOMs" << endl;
    //use wizaed to ask
    wizard.restart();

    int ret = wizard.exec();

    if(ret == QDialog::Rejected)
    {
        wizard.restart();
        return;
    }

    qApp->processEvents();

//    if(MainWindow->getNumberOfTabs() > 1)
//        MainWindow->newTab();
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

void milxQtDICOMPlugin::createActions()
{
    actionOpenSeries = new QAction(MainWindow);
    actionOpenSeries->setIcon(QIcon(":/resources/toolbar/open_series.png"));
    actionOpenSeries->setText(QApplication::translate("MainWindow", "Open Series", 0, QApplication::UnicodeUTF8));
    actionOpenSeries->setShortcut(tr("Ctrl+Alt+o"));
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

  //wizard class var
  wizard.addPage(introPage);
  wizard.addPage(inputPage);
  wizard.addPage(outputPage);
  wizard.setWindowTitle("DICOM Wizard");
}

void milxQtDICOMPlugin::createConnections()
{
//    connect(MainWindow, SIGNAL(windowActivated(QWidget*)), this, SLOT(updateManager(QWidget*)));
    connect(actionOpenSeries, SIGNAL(activated()), MainWindow, SLOT(openSeries()));
    connect(actionConvert, SIGNAL(activated()), this, SLOT(convert()));
    connect(actionAnonymize, SIGNAL(activated()), this, SLOT(anonymize()));

//    connect(btnAtlasName, SIGNAL(clicked()), this, SLOT(showAtlasFileDialog()));
//    connect(btnSurfaceNames, SIGNAL(clicked()), this, SLOT(showSurfacesFileDialog()));
//    connect(btnClearSurfaceNames, SIGNAL(clicked()), comboSurfaceNames, SLOT(clear()));
}

Q_EXPORT_PLUGIN2(DICOMPlugin, milxQtDICOMPluginFactory);
