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
	milxQtThemeEditorForm(milxQtPreferencesForm *theParent = 0, milxQtMain *mainWindow = 0, QString *themeName = 0, bool isEdit = false);
	virtual ~milxQtThemeEditorForm();

	/*!
		\fn milxQtThemeEditorForm::setupEditor(QString *themeName)
		\brief Initialises the theme editor dialog, using the given theme name.
	*/
	void setupEditor(QString *themeName);

public slots:
	/*!
		\fn milxQtThemeEditorForm::accept()
		\brief Overrides the QDialog::accept() function, to allow for invalid name error correction.
	*/
	void accept();

	/*!
		\fn milxQtThemeEditorForm::discardTheme(QAbstractButton *button)
		\brief Discards the changes to the current theme (if new). If an existing theme is being edited, it is
			   deleted instead.
	*/
	void discardTheme(QAbstractButton *button);

	/*!
		\fn milxQtThemeEditorForm::updateStyles(int index)
		\brief Applies the style changes made to the whole application.
	*/
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
	QLabel *widgetBGLabel;
	QComboBox *widgetBGCombo;
	QLabel *pbBGLabel;
	QComboBox *pbBGCombo;
	QLabel *pbTextLabel;
	QComboBox *pbTextCombo;
	QLabel *pbBorderLabel;
	QComboBox *pbBorderCombo;
	QLabel *pbHoverBGLabel;
	QComboBox *pbHoverBGCombo;
	QLabel *pbPressedBGLabel;
	QComboBox *pbPressedBGCombo;
	QLabel *listViewBGLabel;
	QComboBox *listViewBGCombo;
	QLabel *listViewTextLabel;
	QComboBox *listViewTextCombo;
	QLabel *labelTextLabel;
	QComboBox *labelTextCombo;
	QLabel *toolTipBGLabel;
	QComboBox *toolTipBGCombo;
	QLabel *toolTipTextLabel;
	QComboBox *toolTipTextCombo;
	QLabel *toolTipBorderLabel;
	QComboBox *toolTipBorderCombo;
	QLabel *menuBarBGLabel;
	QComboBox *menuBarBGCombo;
	QLabel *menuBarTextLabel;
	QComboBox *menuBarTextCombo;
	QLabel *tbTabBGLabel;
	QComboBox *tbTabBGCombo;
	QLabel *tbTabTextLabel;
	QComboBox *tbTabTextCombo;
	QLabel *tbTabSelectedBGLabel;
	QComboBox *tbTabSelectedBGCombo;
	QLabel *tbTabSelectedTextLabel;
	QComboBox *tbTabSelectedTextCombo;
	QLabel *toolBarBGLabel;
	QComboBox *toolBarBGCombo;
	QLabel *checkBoxTextLabel;
	QComboBox *checkBoxTextCombo;
	QLabel *groupBoxTextLabel;
	QComboBox *groupBoxTextCombo;
	QLabel *msgBoxLabel;
	QComboBox *msgBoxCombo;

	QString getStyleColour(QTextStream *in, QString section);
	void createColourPalette(QComboBox *box);
	void copyThemeFile(QTextStream *in, QString section, QStringList *data);
	void loadColourSettings(QString *theme);
	void saveColourSettings();
	void createConnections();
};

#endif // MILXQTThemeEditorFORM_H
