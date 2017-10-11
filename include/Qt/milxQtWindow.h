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
#ifndef MILXQTSUBWINDOW_H
#define MILXQTSUBWINDOW_H

#include <QtGui/QtGui>
//VTK Headers
#include <QVTKWidget.h>
#include <vtkSmartPointer.h>
#include <vtkPointPicker.h>
//milxQt Specific
#include "milxQtAliases.h"
#include "milxQtConsole.h"

/*!
    \class milxQtWindow
    \brief This class represents the MILX Qt Window Display object using QVTK.
    \author Shekhar S. Chandra, 2013

    This is a base class for all windows. It encapsulates all the basic properties of the windows. You should
    use the milxQtImage, milxQtModel etc. classes instead of using this class directly.

    The rendering is encapsulated within a QVTK widget.
*/
class MILXQT_EXPORT milxQtWindow : public QVTKWidget
{
    Q_OBJECT

public:
    /*!
        \fn milxQtWindow::milxQtWindow(QWidget *theParent = 0)
        \brief The standard constructor
    */
    milxQtWindow(QWidget *theParent = 0);
    /*!
        \fn milxQtWindow::~milxQtWindow()
        \brief The standard destructor
    */
    virtual ~milxQtWindow();

public slots:
    //Naming
    /*!
        \fn milxQtWindow::setName(const QString filename)
        \brief Set the name of the data.
    */
    void setName(const QString filename);
    /*!
        \fn milxQtWindow::getName()
        \brief Returns the name of the data.
    */
    inline QString getName()
    {
        return name;
    }
    /*!
        \fn milxQtWindow::rename()
        \brief Renames the data.
    */
    void rename();
    /*!
        \fn milxQtWindow::strippedName()
        \brief Returns the stripped (path removed) name of the data.
    */
    inline QString strippedName()
    {
        return QFileInfo(name).fileName();
    }
    /*!
        \fn milxQtWindow::strippedBaseName()
        \brief Returns the stripped (path removed) base (no suffix) name of the data.
    */
    inline QString strippedBaseName()
    {
        return QFileInfo(name).baseName();
    }
    /*!
        \fn milxQtWindow::strippedNamePrefix()
        \brief Returns the stripped (path removed) name of the data with "Generic" prefix.

        This member should be overloaded in derived classes to assign correct prefixes.
    */
    virtual inline QString strippedNamePrefix()
    {
        return prefix + QFileInfo(name).fileName();
    }
    /*!
        \fn milxQtWindow::setNamePrefix()
        \brief Sets the prefix of the name of the data provided.
    */
    inline void setNamePrefix(const QString newPrefix)
    {
        prefix = newPrefix;
    }

    /*!
        \fn milxQtWindow::addToContextMenu(QAction* act)
        \brief Adds (prepends) the action to the context menu. Connections are assumed to be made before hand.
    */
    inline void addToContextMenu(QAction *act)
    {
        actionsToAdd.append(act);
    }
    /*!
        \fn milxQtWindow::appendToContextMenu(QAction* act)
        \brief Adds (appends) the action to the context menu. Connections are assumed to be made before hand.
    */
    inline void appendToContextMenu(QAction *act)
    {
        actionsToAppend.append(act);
    }
    /*!
        \fn milxQtWindow::addToContextMenu(QMenu *newMenu)
        \brief Adds the menu to the context menu. Connections are assumed to be made before hand.
    */
    inline void addToContextMenu(QMenu *newMenu)
    {
        addMenuToContextMenu(newMenu);
    }
    /*!
        \fn milxQtWindow::addExtensionAction(QAction* act)
        \brief Adds (in extension section) the action as an extension action. Connections are assumed to be made before hand.
    */
    inline void addExtensionAction(QAction *act)
    {
        extActionsToAdd.append(act);
    }
    /*!
        \fn milxQtWindow::addMenuToContextMenu(QMenu* menus)
        \brief Adds (prepends) the menu to the context menu. Connections are assumed to be made before hand.
    */
    inline void addMenuToContextMenu(QMenu *newMenu)
    {
        menusToAdd.append(newMenu);
    }
    /*!
        \fn milxQtWindow::appendMenuToContextMenu(QMenu* menus)
        \brief Adds (appends) the menu to the context menu. Connections are assumed to be made before hand.
    */
    inline void appendMenuToContextMenu(QMenu *newMenu)
    {
        menusToAppend.append(newMenu);
    }

    /*!
        \fn milxQtWindow::copyToContextMenu(QMenu *copyMenu)
        \brief Copies the menu, by duplicating the entries, to the context menu. Connections are assumed to be made before hand.
    */
    virtual void copyToContextMenu(QMenu *copyMenu);

    /*!
        \fn milxQtWindow::setVerboseMode(bool verbose)
        \brief Verbose mode for message output.
    */
    inline void setVerboseMode(bool verbose)
    {
        verboseMode = verbose;
    }
    /*!
        \fn milxQtWindow::setDeletableOnClose(bool delOnClose)
        \brief Set if the window deletable on close. Default is true.
    */
    inline void setDeletableOnClose(bool delOnClose)
    {
        deletableOnClose = delOnClose;
        setAttribute(Qt::WA_DeleteOnClose, deletableOnClose);
    }
    /*!
        \fn milxQtWindow::isDeletableOnClose()
        \brief Is the window deletable on close?
    */
    inline bool isDeletableOnClose()
    {
        return deletableOnClose;
    }
    /*!
        \fn milxQtWindow::setConsole(milxQtConsole *con)
        \brief Set the console for log output
    */
    inline void setConsole(milxQtConsole *con)
    {
        consoleAssigned = true;
        console = con;
    }

    //Connectors
    /*!
        \fn milxQtWindow::consumeVTKEvent(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
        \brief Consume the event so that VTK interactor style doesn't get it. Safeguard for if remove observer directives have missed some bindings.
    */
    inline void consumeVTKEvent(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
    {
        command->AbortFlagOn();
    }
    /*!
        \brief Update the Qt events, used to keep UI responsive
    */
    inline void updateQtEvents()
    {
        qApp->processEvents();
    }

signals:
    /*!
        \fn milxQtWindow::closing(QWidget *win)
        \brief Send signal that the window is closing.
    */
    void closing(QWidget *win);
    /*!
        \fn milxQtWindow::nameChanged(const QString newName)
        \brief Send signal that the data has been renamed.
    */
    void nameChanged(const QString newName);
    /*!
        \fn milxQtWindow::working(int value)
        \brief Send signal that computation is in progress. Value carries the progress,
    */
    void working(int value);
    /*!
        \fn milxQtWindow::done(int value)
        \brief Send signal that computation is done. Value carries the progress,
    */
    void done(int value);

protected:
    QString name; //!< Name of the data
    QString prefix; //!< Prefix of the data

    //Flags
    bool verboseMode; //!< Verbose message output mode
    bool deletableOnClose; //!< Delete on close allowed? Allowed by default
    bool consoleAssigned; //!< Console assigned for output?

    //Context Menus
    QMenu *contextMenu; //!< Context Menu
    QList<QAction*> actionsToAdd; //!< Context actions to add.
    QList<QAction*> actionsToAppend; //!< Context actions to append.
    QList<QAction*> extActionsToAdd; //!< Extension actions to add.
    QList<QMenu*> menusToAdd; //!< Context Menu's to add.
    QList<QMenu*> menusToAppend; //!< Context Menu's to append.

    milxQtConsole *console; //!< Console for log outputs

    //Internal Members
    /*!
        \fn milxQtWindow::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
        \fn milxQtWindow::closeEvent(QCloseEvent *clEvent)
        \brief When closing, execute this member
    */
    void closeEvent(QCloseEvent *clEvent);
    /*!
        \fn milxQtWindow::printError(QString msg)
        \brief Error message wrapper for console.
    */
    void printError(QString msg);
    /*!
        \fn milxQtWindow::printWarning(QString msg)
        \brief Warning message wrapper for console.
    */
    void printWarning(QString msg);
    /*!
        \fn milxQtWindow::printDebug(QString msg)
        \brief Debug message wrapper for console.
    */
    void printDebug(QString msg);
    /*!
        \fn milxQtWindow::printInfo(QString msg)
        \brief Info message wrapper for console.
    */
    void printInfo(QString msg);

private:

};

#endif // MILXQTSUBWINDOW_H
