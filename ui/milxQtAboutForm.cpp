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
#include "milxQtAboutForm.h"

#include "milxImage.h"
#include "milxModel.h"

milxQtAboutForm::milxQtAboutForm(QWidget *theParent) : QDialog(theParent)
{
	// Setup About dialog box
    ui.setupUi(this);
	//setWindowModality(Qt::ApplicationModal); //block user input
	setWindowTitle(tr("About sMILX"));
	setFixedSize(this->size());
    setupVersion();
    createConnections();
}

milxQtAboutForm::~milxQtAboutForm()
{
    //dtor
}

void milxQtAboutForm::setupVersion()
{
    QString itkStr = QString::number(ITK_VERSION_MAJOR) + "." + QString::number(ITK_VERSION_MINOR) + "." + QString::number(ITK_VERSION_PATCH);
    QString vtkStr = QString::number(VTK_MAJOR_VERSION) + "." + QString::number(VTK_MINOR_VERSION) + "." + QString::number(VTK_BUILD_VERSION);

    ui.aboutEdit->insertPlainText("ITK Version: " + itkStr + "\n");
    ui.aboutEdit->insertPlainText("VTK Version: " + vtkStr + "\n");
    ui.aboutEdit->insertPlainText("SMILI Version: " + QString::number(milx::Version) + "\n");
    ui.aboutEdit->insertPlainText("Qt Version: " + QString(QT_VERSION_STR) + "\n");
    ui.aboutEdit->insertPlainText("milxQt Version: " + QString::number(milxQtVersion));
    ui.aboutEdit->setReadOnly(true);
}

void milxQtAboutForm::createConnections()
{

}
