#ifndef _MILXQTPythonConsole_H
#define _MILXQTPythonConsole_H

/*
 *
 *  Copyright (C) 2010 MeVis Medical Solutions AG All Rights Reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  Further, this software is distributed without any warranty that it is
 *  free of the rightful claim of any third person regarding infringement
 *  or the like.  Any license provided herein, whether implied or
 *  otherwise, applies only to this software file.  Patent licenses, if
 *  any, provided herein do not apply to combinations of this program with
 *  other software, or any other product whatsoever.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Contact information: MeVis Medical Solutions AG, Universitaetsallee 29,
 *  28359 Bremen, Germany or:
 *
 *  http://www.mevis.de
 *
 */

//----------------------------------------------------------------------------------
/*!
// \file    milxQtPythonConsole.h
// \author  Florian Link
// \author  Last changed by $Author: florian $
// \author  Adapted into milxQt by Shekhar Chandra
// \date    2006-10
*/
//----------------------------------------------------------------------------------

#include "PythonQt.h"
#include <QVariant>
#include <QTextEdit>

#include "milxQtAliases.h"

class QCompleter;

//-------------------------------------------------------------------------------
/**
    \class milxQtPythonConsole
    \brief A simple console for python scripting with some syntax highlighting

    Some useful commands:
    \code
    import os

    os.getcwd() # print current working directory
    \endcode

    To run a script (that takes arguments) in the interactive prompt use:
    \code
    >>>import sys
    >>>sys.argv = ['', 'my', 'args', 'here']
    >>>execfile('myscript.py')
    \endcode
*/
class MILXQT_PLUGIN_EXPORT milxQtPythonConsole : public QTextEdit
{
    Q_OBJECT

public:
    milxQtPythonConsole(QWidget* parent, const PythonQtObjectPtr& context, Qt::WindowFlags i = 0);

    ~milxQtPythonConsole();

public slots:
    //! execute current line
    void executeLine(bool storeOnly);

    //! derived key press event
    void keyPressEvent (QKeyEvent * e);

    //! output from console
    void consoleMessage(const QString & message);

    //! get history
    QStringList history()
    {
        return _history;
    }

    //! set history
    void setHistory(const QStringList& h)
    {
        _history = h;
        _historyPosition = 0;
    }

    //! clear the console
    void clear();

    //! overridden to control which characters a user may delete
    virtual void cut();

    //! output redirection
    void stdOut(const QString& s);
    //! output redirection
    void stdErr(const QString& s);

    void insertCompletion(const QString&);

    //! Appends a newline and command prompt at the end of the document.
    void appendCommandPrompt(bool storeOnly = false);

public:
    //! returns true if python cerr had an error
    bool hadError()
    {
        return _hadError;
    }

    //! returns true if python cerr had an error
    void clearError()
    {
        _hadError = false;
    }

protected:
    QAction *copyAct; //!< Copy action for text
    QAction *cutAct; //!< Cut action for text
    QAction *pasteAct; //!< Paste action for text
    QAction *clearAct; //!< clear action for text

    //! handle the pressing of tab
    void handleTabCompletion();

    //! Returns the position of the command prompt
    int commandPromptPosition();

    //! Returns if deletion is allowed at the current cursor
    //! (with and without selected text)
    bool verifySelectionBeforeDeletion();

    //! Sets the current font
    void setCurrentFont(const QColor& color = QColor(0,0,0), bool bold = false);

    //! change the history according to _historyPos
    void changeHistory();

    //! flush output that was not yet printed
    void flushStdOut();

    //! Read settings (mainly history)
    void readSettings();

    //! Write settings (mainly history)
    void writeSettings();

    //! Create all the actions for the console
    void createActions();

    //! Create connections for the console
    void createConnections();

    //! context menu
    void contextMenuEvent(QContextMenuEvent *currentEvent);


private:
    void executeCode(const QString& code);

    PythonQtObjectPtr _context;

    QStringList _history;
    int         _historyPosition;

    QString _clickedAnchor;
    QString _storageKey;
    QString _commandPrompt;

    QString _currentMultiLineCode;

    QString _stdOut;
    QString _stdErr;

    QTextCharFormat _defaultTextCharacterFormat;
    QCompleter* _completer;

    bool _hadError;
};

#endif //_MILXQTPythonConsole_H
