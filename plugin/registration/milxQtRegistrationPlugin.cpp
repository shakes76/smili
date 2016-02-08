#include "milxQtRegistrationPlugin.h"

#include <qplugin.h>

milxQtRegistrationPlugin::milxQtRegistrationPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = false;
    consoleWindow = false;
    extension = true;
    pluginName = "Registration";

    createActions();
    createMenu();
    createConnections();
}

milxQtRegistrationPlugin::~milxQtRegistrationPlugin()
{
    if(isRunning() && threaded)
        quit();
    cout << "Registration Plugin Destroyed." << endl;
}

QString milxQtRegistrationPlugin::name()
{
    return pluginName;
}

QString milxQtRegistrationPlugin::openFileSupport()
{
    QString openPythonExt = "";

    return openPythonExt;
}

QStringList milxQtRegistrationPlugin::openExtensions()
{
    QStringList exts;

    return exts;
}

QStringList milxQtRegistrationPlugin::saveExtensions()
{
    QStringList exts;

    return exts;
}

QString milxQtRegistrationPlugin::saveFileSupport()
{
    QString savePythonExt = "";

    return savePythonExt;
}

void milxQtRegistrationPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{

}

void milxQtRegistrationPlugin::open(QString filename)
{

}

void milxQtRegistrationPlugin::save(QString filename)
{

}

milxQtRenderWindow* milxQtRegistrationPlugin::genericResult()
{
    return NULL;
} //No image result

milxQtModel* milxQtRegistrationPlugin::modelResult()
{
    //~ denoiseModel = new milxQtRegistrationModel;

    return NULL;
    //~ return denoiseModel;
} //No image result

milxQtImage* milxQtRegistrationPlugin::imageResult()
{
    return NULL;
} //No image result

QDockWidget* milxQtRegistrationPlugin::dockWidget()
{
    return NULL;
} //No Dock result

bool milxQtRegistrationPlugin::isPluginWindow(QWidget *window)
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

void milxQtRegistrationPlugin::loadExtension()
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

//~ void milxQtRegistrationPlugin::run()
//~ {
//~ QMutexLocker locker(&mutex); //Lock memory

//~ ///Execute own thread work here

//~ //exec();
//~ }

void milxQtRegistrationPlugin::createActions()
{
    actionItkAffine = new QAction(MainWindow);
    actionItkAffine->setText(QApplication::translate("RegistrationPlugin", "Affine (Itk)", 0, QApplication::UnicodeUTF8));

    actionItkDemon = new QAction(MainWindow);
    actionItkDemon->setText(QApplication::translate("RegistrationPlugin", "Demon (Itk)"));

#ifdef USE_NIFTI_REG
    actionF3DNifti = new QAction(MainWindow);
    actionF3DNifti->setText(QApplication::translate("RegistrationPlugin", "Free Form Deformation (Nifti)", 0, QApplication::UnicodeUTF8));

    actionAladinNifti = new QAction(MainWindow);
    actionAladinNifti->setText(QApplication::translate("RegistrationPlugin", "Aladin (Nifti)"));
#endif


#ifdef USE_ELASTIX
    actionElastixAffine = new QAction(MainWindow);
    actionElastixAffine->setText(QApplication::translate("RegistrationPlugin", "Affine (Elastix)", 0, QApplication::UnicodeUTF8));

    actionElastixBSpline = new QAction(MainWindow);
    actionElastixBSpline->setText(QApplication::translate("RegistrationPlugin", "BSpline (Elastix)"));
#endif

}

void milxQtRegistrationPlugin::createMenu()
{
    menu = new QMenu(MainWindow);
    menu->setTitle(QApplication::translate("RegistrationPlugin", "Image Registration", 0, QApplication::UnicodeUTF8));

    menu->addAction(actionItkAffine);
    menu->addAction(actionItkDemon);

#ifdef USE_NIFTI_REG
    menu->addAction(actionF3DNifti);
    menu->addAction(actionAladinNifti);
#endif

#ifdef USE_ELASTIX
    menu->addAction(actionElastixAffine);
    menu->addAction(actionElastixBSpline);
#endif

    menuToAdd.append(menu);
}


void milxQtRegistrationPlugin::createConnections()
{
    connect(actionItkAffine, SIGNAL(activated()), this, SLOT(ItkAffineRegistrationSlot()));
    connect(actionItkDemon, SIGNAL(activated()), this, SLOT(ItkDemonRegistrationSlot()));

#ifdef USE_NIFTI_REG
    connect(actionF3DNifti, SIGNAL(activated()), this, SLOT(F3DNiftiRegistrationSlot()));
    connect(actionAladinNifti, SIGNAL(activated()), this, SLOT(AladinNiftiRegistrationSlot()));
#endif

#ifdef USE_ELASTIX
    connect(actionElastixAffine, SIGNAL(activated()), this, SLOT(ElastixAffineRegistrationSlot()));
    connect(actionElastixBSpline, SIGNAL(activated()), this, SLOT(ElastixBSplineRegistrationSlot()));
#endif

}

// Itk Affine Registration slot
void milxQtRegistrationPlugin::ItkAffineRegistrationSlot()
{
    regWindow = new milxQtRegistrationWindow(MainWindow);
    regWindow->setAlgo(AffineItk);
    regWindow->show();
}

// Itk Demon Registration slot
void milxQtRegistrationPlugin::ItkDemonRegistrationSlot()
{
    regWindow = new milxQtRegistrationWindow(MainWindow);
    regWindow->setAlgo(DemonItk);
    regWindow->show();
}

#ifdef USE_NIFTI_REG
// F3DNifti Registration slot
void milxQtRegistrationPlugin::F3DNiftiRegistrationSlot()
{
    regWindow = new milxQtRegistrationWindow(MainWindow);
    regWindow->setAlgo(F3DNifti);
    regWindow->show();
}

// AladinNifti Registration slot
void milxQtRegistrationPlugin::AladinNiftiRegistrationSlot()
{
    regWindow = new milxQtRegistrationWindow(MainWindow);
    regWindow->setAlgo(AladinNifti);
    regWindow->show();
}
#endif

#ifdef USE_ELASTIX
// Elastix Affine Registration slot
void milxQtRegistrationPlugin::ElastixAffineRegistrationSlot()
{
    regWindow->setAlgo(ElastixAffine);
    regWindow->show();
}


// Elastix BSpline Registration slot
void milxQtRegistrationPlugin::ElastixBSplineRegistrationSlot()
{
    regWindow->setAlgo(ElastixBSpline);
    regWindow->show();
}
#endif

Q_EXPORT_PLUGIN2(registrationPlugin, milxQtRegistrationPluginFactory);
