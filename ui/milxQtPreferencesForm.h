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

#include <QPushButton>
#include <QSpinBox>
#include <QLabel>

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
	QString currentTheme();
	void addTheme(QString themeName);
	void removeTheme(QString themeName);
	bool isTheme(QString themeName);

public slots:
    void changePage(QListWidgetItem *current, QListWidgetItem *previous);
	void changeTheme(int themeIndex);
	void editTheme();
	void newTheme();
	void accept();

protected:
    Ui::dlgPreferences ui;

    //Options
    //Application
	QSpinBox *windowSizeEdit;
    QSpinBox *processorsEdit;
    QSpinBox *magnifyEdit;
    QCheckBox *timestampCheckBox;
	QComboBox *themeList;
	QPushButton *editThemeButton;
	QPushButton *newThemeButton;
	//View
	QCheckBox *backgroundCheckBox;
	QCheckBox *humanCheckBox;
    //Imaging
    QCheckBox *interpolationCheckBox;
    QCheckBox *orientationCheckBox;
    //Models
    QCheckBox *interpolationModelCheckBox;
    QCheckBox *scalarBarCheckBox;
	QCheckBox *colourMapCheckBox;
	QPushButton *editColourMapButton;
	QPushButton *newColourMapButton;
    //Plugins
	QLabel *noPluginMsg;
    //Streaming
	QComboBox *streamingLevelList;
	QCheckBox *customStreamingCheckBox;
	QSpinBox *customSplitsEdit;

	// The main window
    milxQtMain *MainWindow;

	/*!
		\fn
		\brief
	*/
	void setupPrefs();

	/*!
		\fn
		\brief
	*/
	void loadCustomThemes();

	/*!
		\fn
		\brief
	*/
	void createConnections();
};

#endif // MILXQTPreferencesFORM_H
