/*=========================================================================
  Copyright: (c) CSIRO, Australia

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.
=========================================================================*/

#ifndef MILXQTDICOMPLUGIN_H
#define MILXQTDICOMPLUGIN_H

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
#include "milxQtMain.h"

/**
    \class milxQtDICOMPlugin
    \brief The interface for the DICOM plugin for milxQt
    \author 

    
*/
class MILXQT_PLUGIN_EXPORT milxQtDICOMPlugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    /**
        \fn milxQtDICOMPlugin::milxQtDICOMPlugin(QObject *theParent = 0)
        \brief Default destructor
    */
    milxQtDICOMPlugin(QObject *theParent = 0);
    virtual ~milxQtDICOMPlugin();

    /**
        \fn milxQtDICOMPlugin::name()
        \brief Get the Name of the plugin. [Implement this in your plugin]
    */
    virtual QString name();

    /**
        \fn milxQtDICOMPlugin::hasOpenSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasOpenSupport()
    { return false; }
    /**
        \fn milxQtDICOMPlugin::openFileSupport()
        \brief Get the file support string for opening (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString openFileSupport();
    /**
        \fn milxQtDICOMPlugin::openExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList openExtensions();
    /**
        \fn milxQtDICOMPlugin::hasSaveSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasSaveSupport()
    {   return false;    }
    /**
        \fn milxQtDICOMPlugin::saveFileSupport()
        \brief Get the file support string for saving (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString saveFileSupport();
    /**
        \fn milxQtDICOMPlugin::saveExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList saveExtensions();

    /**
        \fn milxQtDICOMPlugin::hasCollectionSupport()
        \brief Does the plugin support collections (PolyData collection etc.). [Implement this in your plugin]
    */
    inline virtual bool hasCollectionSupport()
    {   return false;    }
    /**
        \fn milxQtDICOMPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
        \brief Pass a collection to internal plugin class. [Implement this in your plugin]
    */
    virtual void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);

    /**
        \fn milxQtDICOMPlugin::open(QString filename)
        \brief Open the file using the plugin. [Implement this in your plugin]
    */
    virtual void open(QString filename);
    /**
        \fn milxQtDICOMPlugin::save(QString filename)
        \brief Save the result as a file using the plugin. [Implement this in your plugin]
    */
    virtual void save(QString filename);

    /**
        \fn milxQtDICOMPlugin::genericResult()
        \brief Get the generic result, which is a milxQtRenderWindow. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtRenderWindow* genericResult();
    /**
        \fn milxQtDICOMPlugin::modelResult()
        \brief Get the model result. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtModel* modelResult();
    /**
        \fn milxQtDICOMPlugin::imageResult()
        \brief Get the image result. The result can then be displayed in milxQtMain etc.[Implement this in your plugin]
    */
    virtual milxQtImage* imageResult();
    /**
        \fn milxQtDICOMPlugin::dockWidget()
        \brief Return the dock widget (if one is provided by plugin). [Implement this in your plugin]
    */
    virtual QDockWidget* dockWidget();

    /**
        \fn milxQtDICOMPlugin::dockDefaultArea()
        \brief Return the default dock widget area (if one is provided by plugin). [Implement this in your plugin]
    */
    inline virtual Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::LeftDockWidgetArea;    }

    /**
        \fn milxQtDICOMPlugin::isPluginWindow(QWidget *window)
        \brief Is the window provided a plugin generated window? In this case a milxQtShapeModel window. [Implement this in your plugin]
    */
    virtual bool isPluginWindow(QWidget *window);

public slots:
    /**
        \fn milxQtDICOMPlugin::loadExtension()
        \brief Load the extension. [Implement this in your plugin]
    */
    virtual void loadExtension();
    /**
        \fn milxQtDICOMPlugin::update()
        \brief Update the plugin. [Implement this in your plugin]

        This generic call is called after plugin is loaded and is designed to be used to update the plugin
        internals such as manager displays etc.
    */
    virtual void update() {}
    /**
        \fn milxQtDICOMPlugin::preStartTasks()
        \brief Tasks to complete before running or starting the thread. [Implement this]
    */
    virtual void preStartTasks() {}
    /**
        \fn milxQtDICOMPlugin::postStartTasks()
        \brief Tasks to complete after running or starting the thread. [Implement this]
    */
    virtual void postStartTasks() {}

    void convert();
    void anonymize();

//    void startWizard();
    void showInputFileDialog();
    void showOutputFileDialog();

protected:
    QPointer<milxQtMain> MainWindow;

    QMenu* menuDICOM; //!< DICOM menu
    //----SSM---- (hierarchical deletion)
    QAction* actionOpenSeries; //!< open series action
    QAction* actionConvert; //!< convert action
    QAction* actionAnonymize; //!< Anonymize action

    //FAI variables for wizard
    QWizard wizard;
    QString inputDirectoryname;
    QLineEdit *txtInputName;
    QString outputDirectoryname;
    QLineEdit *txtOutputName;

    void createActions();
    void createMenu();
    void createWizard();
    void createConnections();

private:

};

class MILXQT_PLUGIN_EXPORT milxQtDICOMPluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQtDICOMPlugin(theParent);  }
};

#endif // MILXQTDICOMPLUGIN_H
