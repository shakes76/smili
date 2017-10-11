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
#include "milxQtConsole.h"

#include <QDateTime>

milxQtConsole::milxQtConsole(QWidget *theParent) : QTabWidget(theParent)
{
    timestamps = true;

    dockArea = Qt::BottomDockWidgetArea;

    dock = new QDockWidget(tr("Console"));
    dock->setFeatures(QDockWidget::AllDockWidgetFeatures);
    dock->setWidget(this);
    dock->setObjectName("Console");

    logWindow = new QTextEdit;
    logWindow->clear();
    logWindow->setReadOnly(true);
    logWindow->setWindowTitle("Log");
    logWindow->setTextInteractionFlags(Qt::TextBrowserInteraction);
    logWindow->ensureCursorVisible();
//    logWindow->setWindowModified(false);
    setTab(logWindow);
    setCurrentWidget(logWindow); //Hierachy deletion

    createActions();

    createConnections();
}

milxQtConsole::~milxQtConsole()
{
    //dtor
}

void milxQtConsole::setTab(QWidget *newWidget)
{
    QTextEdit *editor = qobject_cast<QTextEdit *>(newWidget);
    if(editor)
    {
        connect(editor, SIGNAL(textChanged()), this, SLOT(consoleWasModified()));
        connect(this, SIGNAL(currentChanged(int)), this, SLOT(consoleSwitched(int)));
    }

    addTab(newWidget, newWidget->windowTitle());
}

void milxQtConsole::consoleMessage(const QString & message)
{
    logWindow->append(QString());
    logWindow->insertPlainText(message);
}

void milxQtConsole::consoleHTMLMessage(const QString & message)
{
    logWindow->append(QString());
    logWindow->insertHtml(message);
}

void milxQtConsole::printError(QString msg)
{
    QDateTime currentTime = QDateTime::currentDateTime();
    QString timeStr = "[" + currentTime.toString("ddd-d hh:mm:ss") + "] ";

    if(!timestamps)
        timeStr = "";

    cerr << timeStr.toStdString() << msg.toStdString() << endl;
    //<font color="red">This is some text!</font>
    msg.prepend("Error: ");
    msg.prepend(timeStr);
    msg.prepend("<font color='red'>");
    msg.append("</font>");
    consoleHTMLMessage(msg);
}

void milxQtConsole::printWarning(QString msg)
{
    QDateTime currentTime = QDateTime::currentDateTime();
    QString timeStr = "[" + currentTime.toString("ddd-d hh:mm:ss") + "] ";

    if(!timestamps)
        timeStr = "";

    cerr << timeStr.toStdString() << msg.toStdString() << endl;
    msg.prepend("Warning: ");
    msg.prepend(timeStr);
    msg.prepend("<font color='blue'>");
    msg.append("</font>");
    consoleHTMLMessage(msg);
}

void milxQtConsole::printDebug(QString msg)
{
    QDateTime currentTime = QDateTime::currentDateTime();
    QString timeStr = "[" + currentTime.toString("ddd-d hh:mm:ss") + "] ";

    if(!timestamps)
        timeStr = "";

    cerr << timeStr.toStdString() << msg.toStdString() << endl;
    msg.prepend("Debug: ");
    msg.prepend(timeStr);
    msg.prepend("<font color='orange'>");
    msg.append("</font>");
    consoleHTMLMessage(msg);
}

void milxQtConsole::printInfo(QString msg)
{
    QDateTime currentTime = QDateTime::currentDateTime();
    QString timeStr = "[" + currentTime.toString("ddd-d hh:mm:ss") + "] ";

    if(!timestamps)
        timeStr = "";

    cout << timeStr.toStdString() << msg.toStdString() << endl;
    msg.prepend(timeStr);
    msg.prepend("<font color='black'>");
    msg.append("</font>");
    consoleHTMLMessage(msg);
}

void milxQtConsole::consoleWasModified()
{
    for(int j = 0; j < count(); j ++)
    {
        QTextEdit *editor = qobject_cast<QTextEdit *>(widget(j));
        if(editor && j != currentIndex())
        {
//            editor->setWindowModified(editor->document()->isModified());
            if(editor->document()->isModified())
                setTabText(j, editor->windowTitle() + "*");
        }
    }
}

void milxQtConsole::consoleSwitched(int index)
{
    QTextEdit *editor = qobject_cast<QTextEdit *>(widget(index));
    if(editor)
    {
//        editor->setWindowModified(false);
        editor->document()->setModified(false);
        setTabText(index, editor->windowTitle());
    }
}

void milxQtConsole::createActions()
{
    copyAct = new QAction(this);
        copyAct->setIcon(QIcon(":/resources/toolbar/copy.png"));
        copyAct->setText(tr("Copy", 0));
        copyAct->setShortcut(tr("Ctrl+c"));

//    cutAct = new QAction(this);
//        cutAct->setText(QApplication::translate("Cut", 0, QApplication::UnicodeUTF8));
//        cutAct->setShortcut(tr("Ctrl+x"));
//
//    pasteAct = new QAction(this);
//        pasteAct->setText(QApplication::translate("Paste", 0, QApplication::UnicodeUTF8));
//        pasteAct->setShortcut(tr("Ctrl+v"));

    clearAct = new QAction(this);
        clearAct->setText(tr("Clear", 0));
        clearAct->setShortcut(tr("Crtl+z"));
}

void milxQtConsole::createConnections()
{
    //Actions
    connect(copyAct, SIGNAL(triggered()), logWindow, SLOT(copy()));
//    connect(cutAct, SIGNAL(triggered()), this, SLOT(cut()));
//    connect(pasteAct, SIGNAL(triggered()), this, SLOT(paste()));
    connect(clearAct, SIGNAL(triggered()), logWindow, SLOT(clear()));
}
