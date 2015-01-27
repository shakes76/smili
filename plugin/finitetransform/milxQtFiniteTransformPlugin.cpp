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
#include "milxQtFiniteTransformPlugin.h"

#include <qplugin.h>
//#include <cstdio>

//FTL
extern "C"
{
    #include <nttw/global.h>
    #include <nttw/image.h>
}

milxQtFiniteTransformPlugin::milxQtFiniteTransformPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    //~ MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = false;
    consoleWindow = false;
    extension = false;
    pluginName = "FiniteTransform";
    //~ dataName = "";

    //~ createConnections();
}

milxQtFiniteTransformPlugin::~milxQtFiniteTransformPlugin()
{
    if(isRunning() && threaded)
        quit();
    images.clear();
    cout << "FiniteTransform Plugin Destroyed." << endl;
}

QString milxQtFiniteTransformPlugin::name()
{
    return pluginName;
}

QString milxQtFiniteTransformPlugin::openFileSupport()
{
    QString openExt = "Portable Graymap Files (*.pgm)";

    return openExt;
}

QStringList milxQtFiniteTransformPlugin::openExtensions()
{
    QStringList exts;

    exts.append(".pgm");

    return exts;
}

QStringList milxQtFiniteTransformPlugin::saveExtensions()
{
    QStringList exts;

    exts.append(".pgm");

    return exts;
}

QString milxQtFiniteTransformPlugin::saveFileSupport()
{
    QString savePythonExt = "Portable Graymap Files (*.pgm)";

    return savePythonExt;
}

void milxQtFiniteTransformPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{

}

void milxQtFiniteTransformPlugin::open(QString filename)
{
    nttw_integer *image;
    int rows, cols;
    size_t n;
    cout << "NTTW Library (nttw_integer) Type Size: " << sizeof(nttw_integer) << ", Alignment: " << ALIGNOF(nttw_integer) << " bytes" << endl;

    //Check if binary
    int binaryInFile = isBinaryPGM(filename.toStdString().c_str());
    //load
    if(!readPGM(&image,&rows,&cols,filename.toStdString().c_str(),binaryInFile))
    {
        cerr << "Error Opening File: " << filename.toStdString() << endl;
        return;
    }
    cout << "Opened " << filename.toStdString() << " as FTL image" << endl;
    cout << "Image of size " << rows << "x" << cols << endl;

    //Convert to VNL
    n = static_cast<unsigned>(rows)*static_cast<unsigned>(cols);
    vnl_vector<nttw_integer> imgData(image, n);
    vnl_vector<float> imgFloatData(n);

    for(size_t j = 0; j < imgData.size(); j ++)
        imgFloatData[j] = static_cast<float>(imgData[j]);

    typedef floatImageType::SizeType SizeType;
    SizeType newSize;
    newSize[0] = static_cast<unsigned>(cols);
    newSize[1] = static_cast<unsigned>(rows);
    newSize[2] = 1;

    cout << "Converting to ITK Image" << endl;
    floatImageType::Pointer tmpImageImage = milx::Image<floatImageType>::ImportVectorToImage<float>(imgFloatData, newSize);
    floatImageType::Pointer Image = milx::Image<floatImageType>::DuplicateImage(tmpImageImage);

    QPointer<milxQtImage> img = new milxQtImage;
        img->setName(filename);
        img->setData(Image);
        img->generateImage();

    images.append(img);
}

void milxQtFiniteTransformPlugin::save(QString filename)
{
    typedef floatImageType::SizeType SizeType;
    SizeType newSize = images.last()->GetFloatImage()->GetLargestPossibleRegion().GetSize();
    const size_t n = static_cast<unsigned>(newSize[1])*static_cast<unsigned>(newSize[0]);
    vnl_vector<nttw_integer> imgData(n);

    for(size_t j = 0; j < imgData.size(); j ++)
        imgData[j] = static_cast<nttw_integer>(images.last()->GetFloatImage()->GetBufferPointer()[j]);

    int binaryFile = FALSE;
    writePGM(imgData.data_block(),newSize[1],newSize[0],255,filename.toStdString().c_str(),binaryFile);
}

milxQtRenderWindow* milxQtFiniteTransformPlugin::genericResult()
{
    return NULL;
} //No image result

milxQtModel* milxQtFiniteTransformPlugin::modelResult()
{
    //~ denoiseModel = new milxQtFiniteTransformModel;

    return NULL;
    //~ return denoiseModel;
} //No image result

milxQtImage* milxQtFiniteTransformPlugin::imageResult()
{
    images.last()->generateImage();

    return images.last();
} //image result

QDockWidget* milxQtFiniteTransformPlugin::dockWidget()
{
    return NULL;
} //No Dock result

bool milxQtFiniteTransformPlugin::isPluginWindow(QWidget *window)
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

void milxQtFiniteTransformPlugin::loadExtension()
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

//void milxQtFiniteTransformPlugin::run()
//{
    //~ QMutexLocker locker(&mutex); //Lock memory

    //~ ///Execute own thread work here

    //~ //exec();
//}

void milxQtFiniteTransformPlugin::preStartTasks()
{

}

void milxQtFiniteTransformPlugin::postStartTasks()
{

}

//~ void milxQtFiniteTransformPlugin::createConnections()
//~ {
    //~ //QObject::connect(denoiseAct, SIGNAL(triggered(bool)), denoiseModel, SLOT(denoise()));
//~ }

Q_EXPORT_PLUGIN2(_MYPLUGINCLASS_, milxQtFiniteTransformPluginFactory);
