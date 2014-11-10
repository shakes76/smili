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
#ifndef MILXQTPYTHONPLUGIN_H
#define MILXQTPYTHONPLUGIN_H

#include <QSyntaxHighlighter>

#include <milxQtMain.h>
#include <milxQtFile.h>
#include <milxQtUnifiedWindow.h>
#include <milxQtPluginInterface.h>

#include "milxQtPythonConsole.h"

const int consoleMinHeight = 96;

class MILXQT_PLUGIN_EXPORT milxQtPythonPlugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    milxQtPythonPlugin(QObject *theParent = 0);
    ~milxQtPythonPlugin();

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
    {   return false;    }
    void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);

    void open(QString filename);
    void save(QString filename);

    milxQtRenderWindow* genericResult();
    milxQtModel* modelResult();
    milxQtImage* imageResult();
    QDockWidget* dockWidget();
    inline Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::BottomDockWidgetArea;    }

    bool isPluginWindow(QWidget *window);

    /**
    Casts window to a milxQtPythonConsole class after performing relevant checks
    */
    milxQtPythonConsole* pluginWindow(QWidget *window);

public slots:
    void resizeForFloat(bool floating);
    void addObject(milxQtRenderWindow *rnd)
    {   mainContext.addObject("rnd_" + rnd->strippedBaseName(), rnd);    }
    void addObject(milxQtImage *img)
    {   mainContext.addObject("img_" + img->strippedBaseName(), img);    }
    void addObject(milxQtModel *mdl)
    {   mainContext.addObject("mdl_" + mdl->strippedBaseName(), mdl);    }
    virtual void loadExtension() {}
    virtual void update() {}
    void preStartTasks() {}
    void postStartTasks() {}

protected:
    QPointer<milxQtPythonConsole> pyConsole; //!< Python Console widget
    QPointer<QDockWidget> dock; //!< Dock widget

    QPointer<milxQtMain> MainWindow;

    PythonQtObjectPtr mainContext; //!< Python Context to Python Intepreter.

    void run();
    void setupPython();
    void createConnections();

private:

};

class MILXQT_PLUGIN_EXPORT milxQtPythonPluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQtPythonPlugin(theParent);  }
};

//Syntax Highlighter
/**
    \class milxQtPythonSyntaxHighlighter
    \brief Syntax highlighting for Python to be applied to the console (such as milxQtPythonConsole).

    Escape block is defined as '\@' '@\'.

    Work in progress.
*/
class MILXQT_PLUGIN_EXPORT milxQtPythonSyntaxHighlighter : public QSyntaxHighlighter
{
public:
    milxQtPythonSyntaxHighlighter(QTextEdit *textEdit);

    void highlightBlock(const QString &text);

private:
    struct HighlightingRule
    {
        QRegExp pattern;
        QTextCharFormat format;
    };
    QList<HighlightingRule> highlightingRules;
    QStringList keywords;

    QTextCharFormat keywordFormat;
    QTextCharFormat singleLineCommentFormat;
    QTextCharFormat multiLineCommentFormat;
    QTextCharFormat multiLineEscapeFormat;
    QTextCharFormat quotationFormat;
    QTextCharFormat singleQuotationFormat;

    QRegExp commentStartExpression;
    QRegExp commentEndExpression;
    QRegExp escapeStartExpression;
    QRegExp escapeEndExpression;
};

#endif // MILXQTPYTHONPLUGIN_H
