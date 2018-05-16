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
#include "milxQtControlsForm.h"

milxQtControlsForm::milxQtControlsForm(QWidget *theParent) : QDialog(theParent)
{
	// Setup Controls dialog box
    ui.setupUi(this);
	setWindowTitle(QString("sMILX Controls"));
	setFixedSize(ui.lblControls->size());
}

milxQtControlsForm::~milxQtControlsForm()
{
	//dtor
}

void milxQtControlsForm::closeEvent(QCloseEvent * event)
{
	this->setParent(NULL);
	event->accept();
}