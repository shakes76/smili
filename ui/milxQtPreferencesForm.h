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
#ifndef MILXQTPreferencesFORM_H
#define MILXQTPreferencesFORM_H

#include "ui_preferences.h"
#include "milxQtMain.h"

/*!
    \class milxQtPreferencesForm
    \brief This class represents the Preferences form and other info.

    It has SMILX Logo, author and build libraries and other info.
*/
class MILXQT_EXPORT milxQtPreferencesForm : public QDialog
{
    Q_OBJECT

public:
    milxQtPreferencesForm(milxQtMain *theParent = 0);
    virtual ~milxQtPreferencesForm();

    void setupPages();

public slots:
    void changePage(QListWidgetItem *current, QListWidgetItem *previous);
    void accept();

protected:
    Ui::dlgPreferences ui;

    //Options
    //General
    QCheckBox *backgroundCheckBox;
    QCheckBox *humanCheckBox;
    QSpinBox *windowSizeEdit;
    QSpinBox *processorsEdit;
    QSpinBox *magnifyEdit;
    QCheckBox *timestampCheckBox;
    //Imaging
    QCheckBox *interpolationCheckBox;
    QCheckBox *orientationCheckBox;
    //Models
    QCheckBox *interpolationModelCheckBox;
    QCheckBox *scalarBarCheckBox;
    //Plugins

    milxQtMain *MainWindow;

    void createConnections();
};

#endif // MILXQTPreferencesFORM_H
