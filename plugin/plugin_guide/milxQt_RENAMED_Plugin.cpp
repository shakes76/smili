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
#include "milxQt_RENAMED_Plugin.h"

#include <qplugin.h>

milxQt_RENAMED_Plugin::milxQt_RENAMED_Plugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    //~ MainWindow = qobject_cast<milxQtMain *>(theParent);

    //~ threaded = false;
    //~ dockable = false;
    //~ consoleWindow = false;
    //~ extension = true;
    pluginName = "_RENAMED_";
    //~ dataName = "";

    //~ createConnections();
}

milxQt_RENAMED_Plugin::~milxQt_RENAMED_Plugin()
{
    if(isRunning() && threaded)
        quit();
    cout << "_RENAMED_ Plugin Destroyed." << endl;
}

QString milxQt_RENAMED_Plugin::name()
{
    return pluginName;
}

QString milxQt_RENAMED_Plugin::openFileSupport()
{
    QString openPythonExt = "";

    return openPythonExt;
}

QStringList milxQt_RENAMED_Plugin::openExtensions()
{
    QStringList exts;

    return exts;
}

QStringList milxQt_RENAMED_Plugin::saveExtensions()
{
    QStringList exts;

    return exts;
}

QString milxQt_RENAMED_Plugin::saveFileSupport()
{
    QString savePythonExt = "";

    return savePythonExt;
}

void milxQt_RENAMED_Plugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{

}

void milxQt_RENAMED_Plugin::open(QString filename)
{

}

void milxQt_RENAMED_Plugin::save(QString filename)
{

}

milxQtRenderWindow* milxQt_RENAMED_Plugin::genericResult()
{
    return NULL;
} //No image result

milxQtModel* milxQt_RENAMED_Plugin::modelResult()
{
    //~ denoiseModel = new milxQt_RENAMED_Model;

    return NULL;
    //~ return denoiseModel;
} //No image result

milxQtImage* milxQt_RENAMED_Plugin::imageResult()
{
    return NULL;
} //No image result

QDockWidget* milxQt_RENAMED_Plugin::dockWidget()
{
    return NULL;
} //No Dock result

bool milxQt_RENAMED_Plugin::isPluginWindow(QWidget *window)
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

void milxQt_RENAMED_Plugin::loadExtension()
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

//~ void milxQt_RENAMED_Plugin::run()
//~ {
    //~ QMutexLocker locker(&mutex); //Lock memory

    //~ ///Execute own thread work here

    //~ //exec();
//~ }

//~ void milxQt_RENAMED_Plugin::createConnections()
//~ {
    //~ //QObject::connect(denoiseAct, SIGNAL(triggered(bool)), denoiseModel, SLOT(denoise()));
//~ }

Q_EXPORT_PLUGIN2(_MYPLUGINCLASS_, milxQt_RENAMED_PluginFactory);
