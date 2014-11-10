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
#include "milxQtPythonPlugin.h"

#include <qplugin.h>

milxQtPythonPlugin::milxQtPythonPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = true;
    consoleWindow = true;
    extension = false;
    pluginName = "Python";
    dataName = "";

    setupPython();

    dock = new QDockWidget(tr("Python Console"));
        dock->setFeatures(QDockWidget::AllDockWidgetFeatures);
        dock->setWidget(pyConsole);

    //Syntax Highlight
    milxQtPythonSyntaxHighlighter *highlighter = new milxQtPythonSyntaxHighlighter(pyConsole);
        highlighter->rehighlight();

    createConnections();
}

milxQtPythonPlugin::~milxQtPythonPlugin()
{
    if(isRunning() && threaded)
        quit();
    cerr << "Python Plugin Destroyed." << endl;
}

QString milxQtPythonPlugin::name()
{
    return pluginName;
}

QString milxQtPythonPlugin::openFileSupport()
{
    QString openPythonExt = "";

    return openPythonExt;
}

QStringList milxQtPythonPlugin::openExtensions()
{
    QStringList exts;

    return exts;
}

QStringList milxQtPythonPlugin::saveExtensions()
{
    QStringList exts;

    return exts;
}

QString milxQtPythonPlugin::saveFileSupport()
{
    QString savePythonExt = "";

    return savePythonExt;
}

void milxQtPythonPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{

}

void milxQtPythonPlugin::open(QString filename)
{

}

void milxQtPythonPlugin::save(QString filename)
{

}

milxQtRenderWindow* milxQtPythonPlugin::genericResult()
{
    return NULL;
} //No image result

milxQtModel* milxQtPythonPlugin::modelResult()
{
    return NULL;
} //No image result

milxQtImage* milxQtPythonPlugin::imageResult()
{
    return NULL;
} //No image result

QDockWidget* milxQtPythonPlugin::dockWidget()
{
    pyConsole->show();

    return dock;
}

bool milxQtPythonPlugin::isPluginWindow(QWidget *window)
{
    if(pluginWindow(window) == 0)
        return false;
    else
        return true;
}

milxQtPythonConsole* milxQtPythonPlugin::pluginWindow(QWidget *window)
{
    if(window)
        return qobject_cast<milxQtPythonConsole *>(window);
    return 0;
}

void milxQtPythonPlugin::run()
{
    QMutexLocker locker(&mutex); //Lock memory

    ///Execute own thread work here
    dock->show();

    exec();
}

void milxQtPythonPlugin::setupPython()
{
    PythonQt::init(PythonQt::IgnoreSiteModule | PythonQt::RedirectStdOut);
//    PythonQt_QtAll::init();
    mainContext = PythonQt::self()->getMainModule();

    pyConsole = new milxQtPythonConsole(NULL, mainContext);
//        pyConsole->resize(512,64);
        pyConsole->setWindowTitle("Python");
//        pyConsole->setMaximumHeight(consoleMinHeight);

    // make the object in to the python shell under the name "MainWindow"
    mainContext.addObject("MainWindow", MainWindow);
    milxQtFile *fileIO = new milxQtFile(dock);
    milxQtUnifiedWindow *unified = new milxQtUnifiedWindow(dock); ///\todo remove this when unified window is accessible from Python
    mainContext.addObject("milxQtFile", fileIO); //Add file to python env
    mainContext.addObject("milxQtUnifiedWindow", unified); //Add file to python env
    mainContext.addObject("qApp", qApp); //Add app pointer to python env
}

void milxQtPythonPlugin::createConnections()
{
    QObject::connect(dock, SIGNAL(topLevelChanged(bool)), this, SLOT(resizeForFloat(bool)));
    QObject::connect(MainWindow, SIGNAL(displayed(milxQtRenderWindow*)), this, SLOT(addObject(milxQtRenderWindow*)));
    QObject::connect(MainWindow, SIGNAL(displayed(milxQtImage*)), this, SLOT(addObject(milxQtImage*)));
    QObject::connect(MainWindow, SIGNAL(displayed(milxQtModel*)), this, SLOT(addObject(milxQtModel*)));
}

void milxQtPythonPlugin::resizeForFloat(bool floating)
{
    if(floating)
    {
        pyConsole->resize(512,128);
        pyConsole->setMaximumHeight(QWIDGETSIZE_MAX);
    }
    else
    {
        pyConsole->resize(256,consoleMinHeight);
        pyConsole->setMaximumHeight(consoleMinHeight);
    }
}

Q_EXPORT_PLUGIN2(PythonPlugin, milxQtPythonPluginFactory);

//Syntax class implementation
milxQtPythonSyntaxHighlighter::milxQtPythonSyntaxHighlighter(QTextEdit *textEdit) : QSyntaxHighlighter(textEdit)
{
	// Reserved keywords in Python 2.4
	keywords << "\\band\\b" << "\\bassert\bb" << "\\bbreak\\b" << "\\bclass\\b" << "\\bcontinue\\b" << "\\bdef\\b"
			 << "\\bdel\\b" << "\\belif\\b" << "\\belse\\b" << "\\bexcept\\b" << "\\bexec\\b" << "\\bfinally\\b"
			 << "\\bfor\\b" << "\\bfrom\\b" << "\\bglobal\\b" << "\\bif\\b" << "\\bimport\\b" << "\\bin\\b"
			 << "\\bis\\b" << "\\blambda\\b" << "\\bnot\\b" << "\\bor\\b" << "\\bpass\\b" << "\\bprint\\b" << "\\braise\\b"
			 << "\\breturn\\b" << "\\btry\\b" << "\\bwhile\\b" << "\\byield\\b";

    HighlightingRule rule;

    keywordFormat.setForeground(Qt::darkBlue);
    keywordFormat.setFontWeight(QFont::Bold);

    multiLineCommentFormat.setForeground(Qt::darkRed);
    multiLineEscapeFormat.setForeground(Qt::black);

//    commentStartExpression = QRegExp("'''");
//    commentEndExpression = QRegExp("'''");
    commentStartExpression = QRegExp("/\\*");
    commentEndExpression = QRegExp("\\*/");
    escapeStartExpression = QRegExp("/\\@");
    escapeEndExpression = QRegExp("\\@/");

    ///Add keywords to highlight
    foreach (const QString &pattern, keywords)
    {
         rule.pattern = QRegExp(pattern);
         rule.format = keywordFormat;
         highlightingRules.append(rule);
    }

    ///Quotation highlighting
    quotationFormat.setForeground(Qt::magenta);
    rule.pattern = QRegExp("\".*\"");
    rule.format = quotationFormat;
    highlightingRules.append(rule);

    ///Single Quotation highlighting
    singleQuotationFormat.setForeground(Qt::darkYellow);
    rule.pattern = QRegExp("'.*'");
    rule.format = singleQuotationFormat;
    highlightingRules.append(rule);

    ///Single Line comment highlighting
    singleLineCommentFormat.setForeground(Qt::darkRed);
    rule.pattern = QRegExp("#[^\n]*");
    rule.format = singleLineCommentFormat;
    highlightingRules.append(rule);
}

void milxQtPythonSyntaxHighlighter::highlightBlock(const QString &text)
{
    foreach (const HighlightingRule &rule, highlightingRules)
    {
        QRegExp expression(rule.pattern);
        int index = expression.indexIn(text);
        while (index >= 0)
        {
            int length = expression.matchedLength();
            setFormat(index, length, rule.format);
            index = expression.indexIn(text, index + length);
        }
    }
    setCurrentBlockState(0);

    ///Comment Block
    int startIndex = 0;
    if (previousBlockState() != 1)
        startIndex = commentStartExpression.indexIn(text);

    while (startIndex >= 0)
    {
        int endIndex = commentEndExpression.indexIn(text, startIndex);
        int commentLength;

        if (endIndex == -1)
        {
            setCurrentBlockState(1);
            commentLength = text.length() - startIndex;
        }
        else
            commentLength = endIndex - startIndex + commentEndExpression.matchedLength();

        setFormat(startIndex, commentLength, multiLineCommentFormat);
        startIndex = commentStartExpression.indexIn(text, startIndex + commentLength);
    }

    ///Escape Block
    startIndex = 0;
    if (previousBlockState() != 2)
        startIndex = escapeStartExpression.indexIn(text);

    while (startIndex >= 0)
    {
        int endIndex = escapeEndExpression.indexIn(text, startIndex);
        int commentLength;

        if (endIndex == -1)
        {
            setCurrentBlockState(2);
            commentLength = text.length() - startIndex;
        }
        else
            commentLength = endIndex - startIndex + escapeEndExpression.matchedLength();

        setFormat(startIndex, commentLength, multiLineEscapeFormat);
        startIndex = escapeStartExpression.indexIn(text, startIndex + commentLength);
    }
}
