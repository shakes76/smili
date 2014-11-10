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
#ifndef MILXQTCONSOLE
#define MILXQTCONSOLE

#include <QTabWidget>
#include <QTextEdit>
#include <QDockWidget>
#include <QAction>
#include <QPointer>
//milxQt Specific
#include "milxQtAliases.h"

/**
    \class milxQtConsole
    \brief A console (tabbed) widget class for displaying information such as logs, terminals and consoles.
    \author Shekhar S. Chandra, 2013
*/
class MILXQT_EXPORT milxQtConsole : public QTabWidget
{
    Q_OBJECT

public:
    /** Default constructor */
    milxQtConsole(QWidget *theParent = 0);
    /** Default destructor */
    virtual ~milxQtConsole();

    inline void setTimestamps(const bool timestamp)
    {   timestamps = timestamp;   }
    void setTab(QWidget *newWidget);
    //! output message
    void consoleMessage(const QString & message);
    //! output message as HTML
    void consoleHTMLMessage(const QString & message);

    /**
        \fn milxQtConsole::dockWidget()
        \brief Return the dock widget of the current tabs.
    */
    inline QDockWidget* dockWidget()
    {   return dock;    }
    inline Qt::DockWidgetArea dockDefaultArea()
    {   return dockArea;    }

public slots:
    inline void setDockDefaultArea(Qt::DockWidgetArea area)
    {   dockArea = area;    }
    //Internal Print Members
    /*!
        \fn milxQtConsole::printError(QString msg)
        \brief Error message wrapper for console.
    */
    void printError(QString msg);
    /*!
        \fn milxQtConsole::printWarning(QString msg)
        \brief Warning message wrapper for console.
    */
    void printWarning(QString msg);
    /*!
        \fn milxQtConsole::printDebug(QString msg)
        \brief Debug message wrapper for console.
    */
    void printDebug(QString msg);
    /*!
        \fn milxQtConsole::printInfo(QString msg)
        \brief Info message wrapper for console.
    */
    void printInfo(QString msg);
    /*!
        \fn milxQtConsole::consoleWasModified()
        \brief Change the window title to show that the tab has been modified.
    */
    void consoleWasModified();
    /*!
        \fn milxQtConsole::consoleSwitched(int index)
        \brief Update member for when the tab changes.
    */
    void consoleSwitched(int index);

protected:
    bool timestamps; //!< Timestamp messages?

    QAction *copyAct; //!< Copy action for text
//    QAction *cutAct; //!< Cut action for text
//    QAction *pasteAct; //!< Paste action for text
    QAction *clearAct; //!< clear action for text

    QPointer<QTextEdit> logWindow; //!< Log messages window
    QPointer<QDockWidget> dock; //!< Dock widget

    Qt::DockWidgetArea dockArea;

    //! Create all the actions for the console
    void createActions();
    //! Create connections for the console
    void createConnections();

private:
};

#endif // MILXQTCONSOLE
