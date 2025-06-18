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
#ifndef MILXQTABOUTFORM_H
#define MILXQTABOUTFORM_H

#include "ui_about.h"
#include "milxQtAliases.h"

/*!
    \class milxQtAboutForm
    \brief This class represents the about form and other info.

    It has SMILX Logo, author and build libraries and other info.
*/
class MILXQT_EXPORT milxQtAboutForm : public QDialog
{
    Q_OBJECT

public:
    milxQtAboutForm(QWidget *theParent = 0);
    virtual ~milxQtAboutForm();

    void setupVersion();

protected:
    Ui::dlgAbout ui;

    void createConnections();
};

#endif // MILXQTABOUTFORM_H
