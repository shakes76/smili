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
#ifndef MILXQTPLUGININTERFACE_H
#define MILXQTPLUGININTERFACE_H

#include <QThread>
#include <QMenu>
#include <QDockWidget>
#if (QT_VERSION > 4)
    #include <QtPlugin>
#else
    #include <qplugin.h>
#endif

//VTK
#include <vtkPolyDataCollection.h>

#include "milxQtAliases.h"
#include "milxQtRenderWindow.h"
#include "milxQtImage.h"
#include "milxQtModel.h"

/**
    \class milxQtPluginInterface
    \brief The interface for any plugins that can be made for milxQtMain.
    \author Shekhar S. Chandra, 2013

    sMILX has three types of Plugins:
    - A dock window plugin
    - A full plugin for processing and data manipulation, loading and creation etc.
    - An extension
    An extension is a light-weight plugin that provides just a few content menu entries and some extensions to the basic MILX Qt classes like milxQtImage or milxQtModel.
    Examples of a dock plugin, standard plugin and extension is the milxQtPythonPlugin, milxQtSSMPlugin and milxQtAnimatePlugin classes respectively.

    Members marked with [Implement this in your plugin] should be reimplemented in your plugin but maybe empty if not used.
*/
class MILXQT_EXPORT milxQtPluginInterface : public QThread
{
    Q_OBJECT

public:
    /**
        \fn milxQtPluginInterface::milxQtPluginInterface(QObject *theParent = 0)
        \brief Default destructor
    */
    milxQtPluginInterface(QObject *theParent = 0) : QThread(theParent) {};
    virtual ~milxQtPluginInterface() {};

    /**
        \fn milxQtPluginInterface::name()
        \brief Get the Name of the plugin. [Implement this in your plugin]
    */
    virtual QString name() = 0;
    /**
        \fn milxQtPluginInterface::setFileName(const QString filename)
        \brief Set the Name of the data. [Don't Reimplement this]

        Plugins should not reimplement this member
    */
    inline void setFileName(const QString filename)
    {
        dataName = filename;
    }
    /**
        \fn milxQtPluginInterface::setConsole(milxQtConsole *con)
        \brief Set console to be used by plugins for output

        Plugins should not reimplement this member
    */
    inline void setConsole(milxQtConsole *con)
    {
        console = con;
    }

    /**
        \fn milxQtPluginInterface::hasOpenSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    virtual bool hasOpenSupport() = 0;
    /**
        \fn milxQtPluginInterface::openFileSupport()
        \brief Get the file support string for opening (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString openFileSupport() = 0;
    /**
        \fn milxQtPluginInterface::openExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList openExtensions() = 0;
    /**
        \fn milxQtPluginInterface::hasSaveSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    virtual bool hasSaveSupport() = 0;
    /**
        \fn milxQtPluginInterface::saveFileSupport()
        \brief Get the file support string for saving (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString saveFileSupport() = 0;
    /**
        \fn milxQtPluginInterface::saveExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList saveExtensions() = 0;

    /**
        \fn milxQtPluginInterface::hasCollectionSupport()
        \brief Does the plugin support collections (PolyData collection etc.). [Implement this in your plugin]
    */
    virtual bool hasCollectionSupport() = 0;
    /**
        \fn milxQtPluginInterface::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
        \brief Pass a collection to internal plugin class. [Implement this in your plugin]
    */
    virtual void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames) = 0;

    /**
        \fn milxQtPluginInterface::open(QString filename)
        \brief Open the file using the plugin. [Implement this in your plugin]
    */
    virtual void open(QString filename) = 0;
    /**
        \fn milxQtPluginInterface::save(QString filename)
        \brief Save the result as a file using the plugin. [Implement this in your plugin]
    */
    virtual void save(QString filename) = 0;

    /**
        \fn milxQtPluginInterface::genericResult()
        \brief Get the generic result, which is a milxQtRenderWindow. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtRenderWindow* genericResult() = 0;
    /**
        \fn milxQtPluginInterface::modelResult()
        \brief Get the model result. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtModel* modelResult() = 0;
    /**
        \fn milxQtPluginInterface::imageResult()
        \brief Get the image result. The result can then be displayed in milxQtMain etc.[Implement this in your plugin]
    */
    virtual milxQtImage* imageResult() = 0;
    /**
        \fn milxQtPluginInterface::dockWidget()
        \brief Return the dock widget (if one is provided by plugin). [Implement this in your plugin]
    */
    virtual QDockWidget* dockWidget() = 0;

    /**
        \fn milxQtPluginInterface::dockDefaultArea()
        \brief Return the default dock widget area (if one is provided by plugin). [Implement this in your plugin]
    */
    virtual Qt::DockWidgetArea dockDefaultArea() = 0;

    /**
        \fn milxQtPluginInterface::isPluginWindow(QWidget *window)
        \brief Is the window provided a plugin generated window? In this case a milxQtShapeModel window. [Implement this in your plugin]
    */
    virtual bool isPluginWindow(QWidget *window) = 0;

    /**
        \fn milxQtPluginInterface::isThreaded()
        \brief Is the plugin threaded? [Don't Reimplement this]
    */
    inline bool isThreaded()
    {
        return threaded;
    }
    /**
        \fn milxQtPluginInterface::isDockable()
        \brief Is the plugin a dock window? [Don't Reimplement this]
    */
    inline bool isDockable()
    {
        return dockable;
    }
    /**
        \fn milxQtPluginInterface::isConsole()
        \brief Is the plugin a console dock window? If so, it will be added to the console dock window. [Don't Reimplement this]
    */
    inline bool isConsole()
    {
        return consoleWindow;
    }
    /**
        \fn milxQtPluginInterface::isExtension()
        \brief Is the plugin an extension? [Don't Reimplement this]
    */
    inline bool isExtension()
    {
        return extension;
    }

    /**
        \fn milxQtPluginInterface::addToFileMenu()
        \brief Actions to add to the file or context menu. [Don't Reimplement this]
    */
    inline QList<QAction*> addToFileMenu()
    {
        return fileMenuEntries;
    }
    /**
        \fn milxQtPluginInterface::addToMenuBar()
        \brief Menus to add to the menu bar. [Don't Reimplement this]
    */
    inline QList<QMenu*> addToMenuBar()
    {
        return menuToAdd;
    }

public slots:
    /**
        \fn milxQtPluginInterface::loadExtension()
        \brief Load the extension. [Implement this in your plugin]
    */
    virtual void loadExtension() = 0;
    /**
        \fn milxQtPluginInterface::update()
        \brief Update the plugin. [Implement this in your plugin]

        This generic call is called after plugin is loaded and is designed to be used to update the plugin
        internals such as manager displays etc.
    */
    virtual void update() = 0;
    /**
        \fn milxQtPluginInterface::preStartTasks()
        \brief Tasks to complete before running or starting the thread. [Implement this]
    */
    virtual void preStartTasks() = 0;
    /**
        \fn milxQtPluginInterface::postStartTasks()
        \brief Tasks to complete after running or starting the thread. [Implement this]
    */
    virtual void postStartTasks() = 0;

signals:
    /*!
        \fn milxQtPluginInterface::resultAvailable(milxQtRenderWindow*)
        \brief Send signal that Resultant render window is available for showing.
    */
    void resultAvailable(milxQtRenderWindow*);
    /*!
        \fn milxQtPluginInterface::resultAvailable(milxQtModel*)
        \brief Send signal that Resultant model is available for showing.
    */
    void resultAvailable(milxQtModel*);
    /*!
        \fn milxQtPluginInterface::resultAvailable(milxQtImage*)
        \brief Send signal that Resultant image is available for showing.
    */
    void resultAvailable(milxQtImage*);
    /*!
        \fn milxQtPluginInterface::resultAvailable(vtkPolyDataCollection*, QStringList&)
        \brief Send signal that Resultant collection is available for showing.
    */
    void resultAvailable(vtkPolyDataCollection*, QStringList&);

    /*!
        \fn milxQtPluginInterface::working(int value)
        \brief Send signal that computation is in progress. Value carries the progress,
    */
    void working(int value);
    /*!
        \fn milxQtPluginInterface::done(int value)
        \brief Send signal that computation is done. Value carries the progress,
    */
    void done(int value);

protected:
    QString pluginName;
    QString dataName;
    bool threaded; //!< Threaded plugin?
    bool dockable; //!< Dockable plugin?
    bool consoleWindow; //!< console window?
    bool extension; //!< Extension rather than a plugin?

    QMutex mutex;

    QList<QAction*> fileMenuEntries;
    QList<QMenu*> menuToAdd;

    milxQtConsole *console; //!< console for logging

private:

};

class MILXQT_EXPORT milxQtPluginFactory
{
public:
    virtual ~milxQtPluginFactory() {}
    virtual milxQtPluginInterface* newPlugin(QObject *theParent = 0) = 0;
};

Q_DECLARE_INTERFACE(milxQtPluginFactory, "milxQt.Plugins.milxQtPluginFactory/1.0");

#endif // MILXQTPLUGININTERFACE_H
