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

	regWindow = new milxQtRegistrationWindow(MainWindow);

    createActions();
    createMenu();
    createConnections();
 }

milxQtRegistrationPlugin::~milxQtRegistrationPlugin()
{
	delete regWindow;
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
    actionFFD = new QAction(MainWindow);
	actionFFD->setText(QApplication::translate("RegistrationPlugin", "Free Form Deformation", 0, QApplication::UnicodeUTF8));

	actionAffine = new QAction(MainWindow);
	actionAffine->setText(QApplication::translate("RegistrationPlugin", "Affine (Aladin)"));
}

void milxQtRegistrationPlugin::createMenu()
{
    menu = new QMenu(MainWindow);
    menu->setTitle(QApplication::translate("RegistrationPlugin", "Image Registration", 0, QApplication::UnicodeUTF8));
	menu->addAction(actionFFD);
	menu->addAction(actionAffine);

    menuToAdd.append(menu);
}


void milxQtRegistrationPlugin::createConnections()
{
	connect(actionFFD, SIGNAL(activated()), this, SLOT(ffdRegistrationSlot()));
	connect(actionAffine, SIGNAL(activated()), this, SLOT(affineRegistrationSlot()));
}


// FFD Registration slot
void milxQtRegistrationPlugin::ffdRegistrationSlot()
{
	regWindow->setup(FFD);
	regWindow->show();
}

// Affine Registration slot
void milxQtRegistrationPlugin::affineRegistrationSlot()
{
	regWindow->setup(Affine);
	regWindow->show();
}




Q_EXPORT_PLUGIN2(registrationPlugin, milxQtRegistrationPluginFactory);
