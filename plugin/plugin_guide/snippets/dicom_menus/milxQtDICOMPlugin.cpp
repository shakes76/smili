#include "milxQtDICOMPlugin.h"

#include <qplugin.h>

milxQtDICOMPlugin::milxQtDICOMPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    //~ threaded = false;
    //~ dockable = false;
    //~ consoleWindow = false;
    //~ extension = true;
    pluginName = "DICOM";
    //~ dataName = "";

    createActions();
    createMenu();
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
}

void milxQtDICOMPlugin::anonymize()
{
    cout << "Anonymizing DICOMs" << endl;
}

void milxQtDICOMPlugin::createActions()
{
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
    menuDICOM->addAction(actionConvert);
    menuDICOM->addAction(actionAnonymize);

    menuToAdd.append(menuDICOM);
}

void milxQtDICOMPlugin::createConnections()
{
//    connect(MainWindow, SIGNAL(windowActivated(QWidget*)), this, SLOT(updateManager(QWidget*)));
    connect(actionConvert, SIGNAL(activated()), this, SLOT(convert()));
    connect(actionAnonymize, SIGNAL(activated()), this, SLOT(anonymize()));

//    connect(btnAtlasName, SIGNAL(clicked()), this, SLOT(showAtlasFileDialog()));
//    connect(btnSurfaceNames, SIGNAL(clicked()), this, SLOT(showSurfacesFileDialog()));
//    connect(btnClearSurfaceNames, SIGNAL(clicked()), comboSurfaceNames, SLOT(clear()));
}

Q_EXPORT_PLUGIN2(DICOMPlugin, milxQtDICOMPluginFactory);
