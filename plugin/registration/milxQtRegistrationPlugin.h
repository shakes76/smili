#ifndef MILXQTRegistrationPLUGIN_H
#define MILXQTRegistrationPLUGIN_H

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
#include "milxQtRegistrationWindow.h"


/**
    \class milxQtRegistrationPlugin
    \brief The interface for the Registration plugin for milxQt
    \author 

    
*/
class MILXQT_PLUGIN_EXPORT milxQtRegistrationPlugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    /**
        \fn milxQtRegistrationPlugin::milxQtRegistrationPlugin(QObject *theParent = 0)
        \brief Default destructor
    */
    milxQtRegistrationPlugin(QObject *theParent = 0);
    virtual ~milxQtRegistrationPlugin();

    /**
        \fn milxQtRegistrationPlugin::name()
        \brief Get the Name of the plugin. [Implement this in your plugin]
    */
    virtual QString name();

    /**
        \fn milxQtRegistrationPlugin::hasOpenSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasOpenSupport()
    { return false; }
    /**
        \fn milxQtRegistrationPlugin::openFileSupport()
        \brief Get the file support string for opening (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString openFileSupport();
    /**
        \fn milxQtRegistrationPlugin::openExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList openExtensions();
    /**
        \fn milxQtRegistrationPlugin::hasSaveSupport()
        \brief Does the plugin support opening files? [Implement this in your plugin]
    */
    inline virtual bool hasSaveSupport()
    {   return false;    }
    /**
        \fn milxQtRegistrationPlugin::saveFileSupport()
        \brief Get the file support string for saving (extension wildcard list). [Implement this in your plugin]
    */
    virtual QString saveFileSupport();
    /**
        \fn milxQtRegistrationPlugin::saveExtensions()
        \brief Get a list of supported file format extensions. [Implement this in your plugin]
    */
    virtual QStringList saveExtensions();

    /**
        \fn milxQtRegistrationPlugin::hasCollectionSupport()
        \brief Does the plugin support collections (PolyData collection etc.). [Implement this in your plugin]
    */
    inline virtual bool hasCollectionSupport()
    {   return false;    }
    /**
        \fn milxQtRegistrationPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
        \brief Pass a collection to internal plugin class. [Implement this in your plugin]
    */
    virtual void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);

    /**
        \fn milxQtRegistrationPlugin::open(QString filename)
        \brief Open the file using the plugin. [Implement this in your plugin]
    */
    virtual void open(QString filename);
    /**
        \fn milxQtRegistrationPlugin::save(QString filename)
        \brief Save the result as a file using the plugin. [Implement this in your plugin]
    */
    virtual void save(QString filename);

    /**
        \fn milxQtRegistrationPlugin::genericResult()
        \brief Get the generic result, which is a milxQtRenderWindow. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtRenderWindow* genericResult();
    /**
        \fn milxQtRegistrationPlugin::modelResult()
        \brief Get the model result. The result can then be displayed in milxQtMain etc. [Implement this in your plugin]
    */
    virtual milxQtModel* modelResult();
    /**
        \fn milxQtRegistrationPlugin::imageResult()
        \brief Get the image result. The result can then be displayed in milxQtMain etc.[Implement this in your plugin]
    */
    virtual milxQtImage* imageResult();
    /**
        \fn milxQtRegistrationPlugin::dockWidget()
        \brief Return the dock widget (if one is provided by plugin). [Implement this in your plugin]
    */
    virtual QDockWidget* dockWidget();

    /**
        \fn milxQtRegistrationPlugin::dockDefaultArea()
        \brief Return the default dock widget area (if one is provided by plugin). [Implement this in your plugin]
    */
    inline virtual Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::LeftDockWidgetArea;    }

    /**
        \fn milxQtRegistrationPlugin::isPluginWindow(QWidget *window)
        \brief Is the window provided a plugin generated window? In this case a milxQtShapeModel window. [Implement this in your plugin]
    */
    virtual bool isPluginWindow(QWidget *window);

	/**
		\fn milxQtRegistrationPlugin::registration(RegType type)
		\brief Register the opened images with the selected algorithm
	*/
    void registration(RegType type);

public slots:
    /**
        \fn milxQtRegistrationPlugin::loadExtension()
        \brief Load the extension. [Implement this in your plugin]
    */
    virtual void loadExtension();
    /**
        \fn milxQtRegistrationPlugin::update()
        \brief Update the plugin. [Implement this in your plugin]

        This generic call is called after plugin is loaded and is designed to be used to update the plugin
        internals such as manager displays etc.
    */
    virtual void update() {}
    /**
        \fn milxQtRegistrationPlugin::preStartTasks()
        \brief Tasks to complete before running or starting the thread. [Implement this]
    */
    virtual void preStartTasks() {}
    /**
        \fn milxQtRegistrationPlugin::postStartTasks()
        \brief Tasks to complete after running or starting the thread. [Implement this]
    */
    virtual void postStartTasks() {}

	/**
		\fn ffdRegistrationSlot()
		\brief Slot for the ffd registration button
	*/
	void ffdRegistrationSlot();

	/**
		\fn affineRegistrationSlot()
		\brief Slot for the affine registration button
	*/
	void affineRegistrationSlot();

protected:
    
    QPointer<milxQtMain> MainWindow;

    QMenu* menu; //!< Registration menu
    QAction* actionFFD; //!< FFD registration action
	QAction* actionAffine; // Affine registration action
	milxQtRegistrationWindow * regWindow; // registration window

    void createActions();
    void createMenu();
    void createConnections();

};

class MILXQT_PLUGIN_EXPORT milxQtRegistrationPluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQtRegistrationPlugin(theParent);  }
};

#endif // MILXQTRegistrationPLUGIN_H
