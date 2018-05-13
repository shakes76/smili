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
#ifndef MILXQTThemeEditorFORM_H
#define MILXQTThemeEditorFORM_H

#include <QDialog>
#include "ui_themeEditor.h"
#include "milxQtPreferencesForm.h"

/*!
\class milxQtThemeEditorForm
\brief This class represents the Theme Editor form.
*/
class MILXQT_EXPORT milxQtThemeEditorForm : public QDialog
{
    Q_OBJECT

public:
	milxQtThemeEditorForm(milxQtPreferencesForm *theParent = 0, milxQtMain *mainWindow = 0, QString *themeName = 0, bool isEdit = 0);
	virtual ~milxQtThemeEditorForm();

	void setupEditor(QString *themeName);

public slots:
	void accept();
	void discardTheme(QAbstractButton *button);
	void updateStyles(int index);

protected:
    Ui::dlgThemeEditor ui;
	milxQtMain *MainWindow;
	milxQtPreferencesForm *prefForm;
	bool isEdit;

	// The theme title
	QLabel *themeNameLabel;
	QLineEdit *themeNameInput;
	QString *themeName;
	
	// The app theme elements
	QLabel *widgetLabel;
	QComboBox *widgetCombo;
	QLabel *pbLabel;
	QComboBox *pbCombo;
	QLabel *pbHoverLabel;
	QComboBox *pbHoverCombo;
	QLabel *pbPressedLabel;
	QComboBox *pbPressedCombo;
	QLabel *listViewLabel;
	QComboBox *listViewCombo;
	QLabel *toolTipLabel;
	QComboBox *toolTipCombo;
	QLabel *menuBarLabel;
	QComboBox *menuBarCombo;
	QLabel *tbTabLabel;
	QComboBox *tbTabCombo;
	QLabel *tbTabTextLabel;
	QComboBox *tbTabTextCombo;
	QLabel *tbTabSelectedLabel;
	QComboBox *tbTabSelectedCombo;
	QLabel *tbTabSelectedTextLabel;
	QComboBox *tbTabSelectedTextCombo;
	QLabel *toolBarLabel;
	QComboBox *toolBarCombo;
	
	QString getStyleColour(QTextStream *in, QString section);
	void createColourPalette(QComboBox *box);
	void copyThemeFile(QTextStream *in, QString section, QStringList *data);
	void loadColourSettings(QString *theme);
	void saveColourSettings();
	void createConnections();
};

#endif // MILXQTThemeEditorFORM_H
