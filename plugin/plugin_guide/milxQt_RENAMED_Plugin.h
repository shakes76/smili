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
#ifndef MILXQT_RENAMED_PLUGIN_H
#define MILXQT_RENAMED_PLUGIN_H

#include <QThread>
#include <QMenu>
#include <QDockWidget>
//VTK
#include <vtkPolyDataCollection.h>

#include "milxQtAliases.h"
#include "milxQtPluginInterface.h"
#include "milxQtRenderWindow.h"
#include "milxQtImage.h"
#include "milxQtModel.h"

/**
    \class milxQt_RENAMED_Plugin
    \brief The interface for the _RENAMED_ plugin for milxQt
    \author 

    
*/
class MILXQT_PLUGIN_EXPORT milxQt_RENAMED_Plugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    /**
        \fn milxQt_RENAMED_Plugin::milxQt_RENAMED_Plugin(QObject *theParent = 0)
        \brief Default destructor
    */
    milxQt_RENAMED_Plugin(QObject *theParent = 0);
    virtual ~milxQt_RENAMED_Plugin();

    /**
        \fn milxQt_RENAMED_Plugin::name()
        \brief Get the Name of the plugin. [Implement this in your plugin]
    */
    virtual QString name();

    /**
        \fn milxQt_RENAMED_Plugin::hasOpenSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasOpenSupport()
    { return false; }
    /**
        \fn milxQt_RENAMED_Plugin::openFileSupport()
        \brief Get the file support string for opening (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString openFileSupport();
    /**
        \fn milxQt_RENAMED_Plugin::openExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList openExtensions();
    /**
        \fn milxQt_RENAMED_Plugin::hasSaveSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasSaveSupport()
    {   return false;    }
    /**
        \fn milxQt_RENAMED_Plugin::saveFileSupport()
        \brief Get the file support string for saving (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString saveFileSupport();
    /**
        \fn milxQt_RENAMED_Plugin::saveExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList saveExtensions();

    /**
        \fn milxQt_RENAMED_Plugin::hasCollectionSupport()
        \brief Does the plugin support collections (PolyData collection etc.). [Implement this in your plugin]
    */
    inline virtual bool hasCollectionSupport()
    {   return false;    }
    /**
        \fn milxQt_RENAMED_Plugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
        \brief Pass a collection to internal plugin class. [Implement this in your plugin]
    */
    virtual void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);

    /**
        \fn milxQt_RENAMED_Plugin::open(QString filename)
        \brief Open the file using the plugin. [Implement this in your plugin]
    */
    virtual void open(QString filename);
    /**
        \fn milxQt_RENAMED_Plugin::save(QString filename)
        \brief Save the result as a file using the plugin. [Implement this in your plugin]
    */
    virtual void save(QString filename);

    /**
        \fn milxQt_RENAMED_Plugin::genericResult()
        \brief Get the generic result, which is a milxQtRenderWindow. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtRenderWindow* genericResult();
    /**
        \fn milxQt_RENAMED_Plugin::modelResult()
        \brief Get the model result. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtModel* modelResult();
    /**
        \fn milxQt_RENAMED_Plugin::imageResult()
        \brief Get the image result. The result can then be displayed in milxQtMain etc.[Implement this in your plugin]
    */
    virtual milxQtImage* imageResult();
    /**
        \fn milxQt_RENAMED_Plugin::dockWidget()
        \brief Return the dock widget (if one is provided by plugin). [Implement this in your plugin]
    */
    virtual QDockWidget* dockWidget();

    /**
        \fn milxQt_RENAMED_Plugin::dockDefaultArea()
        \brief Return the default dock widget area (if one is provided by plugin). [Implement this in your plugin]
    */
    inline virtual Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::LeftDockWidgetArea;    }

    /**
        \fn milxQt_RENAMED_Plugin::isPluginWindow(QWidget *window)
        \brief Is the window provided a plugin generated window? In this case a milxQtShapeModel window. [Implement this in your plugin]
    */
    virtual bool isPluginWindow(QWidget *window);

public slots:
    /**
        \fn milxQt_RENAMED_Plugin::loadExtension()
        \brief Load the extension. [Implement this in your plugin]
    */
    virtual void loadExtension();
    /**
        \fn milxQt_RENAMED_Plugin::update()
        \brief Update the plugin. [Implement this in your plugin]

        This generic call is called after plugin is loaded and is designed to be used to update the plugin
        internals such as manager displays etc.
    */
    virtual void update() {}
    /**
        \fn milxQt_RENAMED_Plugin::preStartTasks()
        \brief Tasks to complete before running or starting the thread. [Implement this]
    */
    virtual void preStartTasks() {}
    /**
        \fn milxQt_RENAMED_Plugin::postStartTasks()
        \brief Tasks to complete after running or starting the thread. [Implement this]
    */
    virtual void postStartTasks() {}

protected:
    

private:

};

class MILXQT_PLUGIN_EXPORT milxQt_RENAMED_PluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQt_RENAMED_Plugin(theParent);  }
};

#endif // MILXQT_RENAMED_PLUGIN_H
