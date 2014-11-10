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
#include "milxQtDiffusionTensorPlugin.h"

#include <qplugin.h>

milxQtDiffusionTensorPlugin::milxQtDiffusionTensorPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = false;
    consoleWindow = false;
    extension = true;
    pluginName = "Diffusion Imaging";
    dataName = "";

    //~ createConnections();
}

milxQtDiffusionTensorPlugin::~milxQtDiffusionTensorPlugin()
{
    if(isRunning() && threaded)
        quit();
    cout << "DiffusionTensor Plugin Destroyed." << endl;
}

QString milxQtDiffusionTensorPlugin::name()
{
    return pluginName;
}

QString milxQtDiffusionTensorPlugin::openFileSupport()
{
    QString openPythonExt = "";

    return openPythonExt;
}

QStringList milxQtDiffusionTensorPlugin::openExtensions()
{
    QStringList exts;

    return exts;
}

QStringList milxQtDiffusionTensorPlugin::saveExtensions()
{
    QStringList exts;

    return exts;
}

QString milxQtDiffusionTensorPlugin::saveFileSupport()
{
    QString savePythonExt = "";

    return savePythonExt;
}

void milxQtDiffusionTensorPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{

}

void milxQtDiffusionTensorPlugin::open(QString filename)
{

}

void milxQtDiffusionTensorPlugin::save(QString filename)
{

}

milxQtRenderWindow* milxQtDiffusionTensorPlugin::genericResult()
{
    return NULL;
} //No render result

milxQtModel* milxQtDiffusionTensorPlugin::modelResult()
{
    diffusionModel = new milxQtDiffusionTensorModel;

    return diffusionModel;
}

milxQtImage* milxQtDiffusionTensorPlugin::imageResult()
{
    diffusionImage = new milxQtDiffusionTensorImage;

    return diffusionImage;
}

QDockWidget* milxQtDiffusionTensorPlugin::dockWidget()
{
    return NULL;
} //No Dock result

bool milxQtDiffusionTensorPlugin::isPluginWindow(QWidget *window)
{
    if(pluginWindow(window) == 0)
        return false;
    else
        return true;
}

milxQtDiffusionTensorModel* milxQtDiffusionTensorPlugin::pluginWindow(QWidget *window)
{
    if(window)
        return qobject_cast<milxQtDiffusionTensorModel *>(window);
    return 0;
}

void milxQtDiffusionTensorPlugin::loadExtension()
{
    if(MainWindow->isActiveModel())
    {
        milxQtModel *currentWin = MainWindow->activeModel();
        diffusionModel = new milxQtDiffusionTensorModel; //hierarchical deletion
            diffusionModel->setName(currentWin->getName());
            diffusionModel->SetInput(currentWin->GetOutput());
            diffusionModel->generateModel();

        MainWindow->display(diffusionModel);
    }
    else if(MainWindow->isActiveImage())
    {
        milxQtImage *currentWin = MainWindow->activeImage();
        diffusionImage = new milxQtDiffusionTensorImage; //hierarchical deletion
            diffusionImage->setName(currentWin->getName());
            diffusionImage->setData(currentWin->GetVectorImage(), false); ///\todo assuming vector image there, checks needed
            diffusionImage->generateImage();
            diffusionImage->disableInterpolateDisplay();

        QObject::connect(diffusionImage, SIGNAL(resultAvailable(milxQtModel*)), MainWindow, SLOT(display(milxQtModel*)));
        QObject::connect(diffusionImage, SIGNAL(resultAvailable(milxQtRenderWindow*)), MainWindow, SLOT(display(milxQtRenderWindow*)));

        MainWindow->display(diffusionImage);
    }
    else
        return;
}

//~ void milxQtDiffusionTensorPlugin::run()
//~ {
    //~ QMutexLocker locker(&mutex); //Lock memory

    //~ ///Execute own thread work here

    //~ //exec();
//~ }

//~ void milxQtDiffusionTensorPlugin::createConnections()
//~ {
    //~ //QObject::connect(denoiseAct, SIGNAL(triggered(bool)), denoiseModel, SLOT(denoise()));
//~ }

Q_EXPORT_PLUGIN2(DTIPlugin, milxQtDiffusionTensorPluginFactory);
