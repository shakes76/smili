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
#include "milxQtAnimatePlugin.h"

#include <qplugin.h>

milxQtAnimatePlugin::milxQtAnimatePlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = false;
    consoleWindow = false;
    extension = false;
    pluginName = "Animate";
    dataName = "";

    createConnections();
}

milxQtAnimatePlugin::~milxQtAnimatePlugin()
{
    if(isRunning() && threaded)
        quit();
    cerr << "Animate Plugin Destroyed." << endl;
}

QString milxQtAnimatePlugin::name()
{
    return pluginName;
}

QString milxQtAnimatePlugin::openFileSupport()
{
    QString openPythonExt = "";

    return openPythonExt;
}

QStringList milxQtAnimatePlugin::openExtensions()
{
    QStringList exts;

    return exts;
}

QStringList milxQtAnimatePlugin::saveExtensions()
{
    QStringList exts;

    return exts;
}

QString milxQtAnimatePlugin::saveFileSupport()
{
    QString savePythonExt = "";

    return savePythonExt;
}

void milxQtAnimatePlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{
    cout << "Loading Collection via the Animate Plugin" << endl;
    animateModel = new milxQtAnimateModel(MainWindow);
        animateModel->SetInputCollection(collection, filenames);
        animateModel->setConsole(console);
        cout << "Loaded Collection as an animation." << endl;

    ///Update case list in manager
//    QStringList headingList;
//
//    headingList << "Index" << "Case ID"; //!< Case Browser
//    if(!caseManagerCreated)
//    {
//        caseTabIndex = manager->newTab("Case Browser", headingList);
//        caseManagerCreated = true;
//    }
//    else
//        manager->clearTab(caseTabIndex);
//
//    cout << "Loading cases into browser." << endl;
//    QList< int > cases = shapes.last()->getCaseIDs();
//    for(int j = 0; j < cases.size(); j ++)
//    {
//        QStringList caseList;
//
//        caseList << QString::number(j) << "Case " + QString::number(cases[j]);
//
//        manager->addItem(caseTabIndex, caseList, Qt::NoItemFlags);
//    }
//
//    manager->show();
}

void milxQtAnimatePlugin::open(QString filename)
{

}

void milxQtAnimatePlugin::save(QString filename)
{

}

milxQtRenderWindow* milxQtAnimatePlugin::genericResult()
{
    return NULL;
} //No image result

milxQtModel* milxQtAnimatePlugin::modelResult()
{
    return animateModel;
} //No image result

milxQtImage* milxQtAnimatePlugin::imageResult()
{
    return NULL;
} //No image result

QDockWidget* milxQtAnimatePlugin::dockWidget()
{
    return NULL;
} //No Dock result

bool milxQtAnimatePlugin::isPluginWindow(QWidget *window)
{
//    if(pluginWindow(window) == 0)
//        return false;
//    else
        return true;
}

milxQtAnimateModel* milxQtAnimatePlugin::pluginWindow(QWidget *window)
{
    if(window)
        return qobject_cast<milxQtAnimateModel *>(window);
    return 0;
}

void milxQtAnimatePlugin::loadExtension()
{

}

void milxQtAnimatePlugin::run()
{
    QMutexLocker locker(&mutex); //Lock memory

    ///Execute own thread work here

    //exec();
}

void milxQtAnimatePlugin::createConnections()
{
    //QObject::connect(denoiseAct, SIGNAL(triggered(bool)), denoiseModel, SLOT(denoise()));
}

Q_EXPORT_PLUGIN2(animatePlugin, milxQtAnimatePluginFactory);
