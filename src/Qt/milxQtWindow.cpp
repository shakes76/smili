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
#include "milxQtWindow.h"

#include <QFileDialog>
///VTK Includes
#include <vtkRenderWindow.h>

milxQtWindow::milxQtWindow(QWidget *theParent) : QVTKWidget(theParent)
{
    ///Init Flags
    verboseMode = true;
    deletableOnClose = true;
    consoleAssigned = false;

    ///Set strings
    prefix = "";
    name = "";

	//Turn off the warnings window
	//vtkObject::GlobalWarningDisplayOff();

    createConnections();
}

milxQtWindow::~milxQtWindow()
{
    //dtor
}

void milxQtWindow::setName(const QString filename)
{
    name = filename;
    setWindowTitle(strippedNamePrefix());

    emit nameChanged(strippedNamePrefix());
}

void milxQtWindow::rename()
{
    bool ok;

    QString newName = QInputDialog::getText(this, tr("Rename Data"), tr("New Name: "), QLineEdit::Normal, name, &ok);

    if(ok)
        setName(newName);
}

void milxQtWindow::copyToContextMenu(QMenu *copyMenu)
{
    QList<QAction *> actionsInMenu = copyMenu->actions();
    QMenu *newMenu = new QMenu(this);

    newMenu->setTitle(copyMenu->title());
    foreach(QAction *action, actionsInMenu)
    {
        QAction *newAction = new QAction(action->parent());
        newAction->setText(action->text());
        newAction->setCheckable(action->isCheckable());
        newAction->setChecked(action->isChecked());

        newMenu->addAction(newAction);
    }

    addToContextMenu(newMenu);
}

//Internal Members
void milxQtWindow::createConnections()
{

}

void milxQtWindow::closeEvent(QCloseEvent *clEvent)
{
    clEvent->accept();
    emit closing(this);
}

void milxQtWindow::printError(QString msg)
{
    if(!verboseMode)
        return;

    if(consoleAssigned)
        console->printError(msg);
    else
        cerr << "ERROR: " << msg.toStdString() << endl;
}

void milxQtWindow::printWarning(QString msg)
{
    if(!verboseMode)
        return;

    if(consoleAssigned)
        console->printWarning(msg);
    else
        cerr << "Warning: " << msg.toStdString() << endl;
}

void milxQtWindow::printDebug(QString msg)
{
    if(!verboseMode)
        return;

    if(consoleAssigned)
        console->printDebug(msg);
    else
        cerr << "Debug: " << msg.toStdString() << endl;
}

void milxQtWindow::printInfo(QString msg)
{
    if(!verboseMode)
        return;

    if(consoleAssigned)
        console->printInfo(msg);
    else
        cerr << msg.toStdString() << endl;
}
