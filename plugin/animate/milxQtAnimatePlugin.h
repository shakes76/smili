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
#ifndef MILXQTANIMATEPLUGIN_H
#define MILXQTANIMATEPLUGIN_H

#include <milxQtMain.h>
#include <milxQtPluginInterface.h>

#include "milxQtAnimateModel.h"

/**
    \class milxQtAnimatePlugin
    \brief This is a plugin for milxQt to do animation from collections
*/
class MILXQT_PLUGIN_EXPORT milxQtAnimatePlugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    milxQtAnimatePlugin(QObject *theParent = 0);
    virtual ~milxQtAnimatePlugin();

    QString name();

    inline bool hasOpenSupport()
    {   return false;    }
    QString openFileSupport();
    QStringList openExtensions();
    inline bool hasSaveSupport()
    {   return false;    }
    QString saveFileSupport();
    QStringList saveExtensions();

    inline bool hasCollectionSupport()
    {   return true;    }
    void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);

    void open(QString filename);
    void save(QString filename);

    milxQtRenderWindow* genericResult();
    milxQtModel* modelResult();
    milxQtImage* imageResult();
    QDockWidget* dockWidget();
    inline Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::LeftDockWidgetArea;    }

    bool isPluginWindow(QWidget *window);

    /**
    Casts window to a milxQtDeNoiseConsole class after performing relevant checks
    */
    milxQtAnimateModel* pluginWindow(QWidget *window);

public slots:
    /**
    Load the extension
    */
    void loadExtension();
    void update() {}
    void preStartTasks() {}
    void postStartTasks() {}

protected:
    QPointer<milxQtAnimateModel> animateModel; //!< Model extension
    QPointer<milxQtMain> MainWindow;

    void run();
    void createConnections();

private:

};

class MILXQT_PLUGIN_EXPORT milxQtAnimatePluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "milxQt.Plugins.milxQtPluginFactory/1.0")
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQtAnimatePlugin(theParent);  }
};

#endif // MILXQTANIMATEPLUGIN_H
