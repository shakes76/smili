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
#include "milxQtThemeEditor.h"

#include <cstdlib>
#include <iostream>
#include <QLineEdit>
#include <QScrollArea>

milxQtThemeEditorForm::milxQtThemeEditorForm(milxQtPreferencesForm *theParent, milxQtMain *mainWindow, QString *themeName, bool isEdit) : QDialog(theParent)
{
	ui.setupUi(this);

	setWindowModality(Qt::ApplicationModal); //block user input
	setWindowTitle(tr("Theme Editor"));
	setFixedSize(this->size());

	MainWindow = mainWindow;
	prefForm = theParent;
	this->themeName = themeName;
	if (!themeName) themeName = new QString("New Theme");
	this->isEdit = isEdit;
	setupEditor(themeName);
	createConnections();
}

milxQtThemeEditorForm::~milxQtThemeEditorForm()
{
	//dtor
}

void milxQtThemeEditorForm::setupEditor(QString *themeName)
{
	// The theme name
	themeNameLabel = new QLabel("Theme Name:");
	themeNameInput = new QLineEdit("New Theme");
	if (themeName) {
		themeNameInput->setText(*themeName);
	}
	QHBoxLayout *nameLayout = new QHBoxLayout;
		nameLayout->addWidget(themeNameLabel);
		nameLayout->addWidget(themeNameInput);

	// Default widget background
	widgetBGLabel = new QLabel(tr("Background - Default:"));
	widgetBGCombo = new QComboBox;
	createColourPalette(widgetBGCombo);
	QHBoxLayout *widgetBGLayout = new QHBoxLayout;
		widgetBGLayout->addWidget(widgetBGLabel);
		widgetBGLayout->addWidget(widgetBGCombo);
	
	// Push Button background
	pbBGLabel = new QLabel(tr("Button - Background:"));
	pbBGCombo = new QComboBox;
	createColourPalette(pbBGCombo);
	QHBoxLayout *pbBGLayout = new QHBoxLayout;
		pbBGLayout->addWidget(pbBGLabel);
		pbBGLayout->addWidget(pbBGCombo);
	
	// Push Button text
	pbTextLabel = new QLabel(tr("Button - Text:"));
	pbTextCombo = new QComboBox;
	createColourPalette(pbTextCombo);
	QHBoxLayout *pbTextLayout = new QHBoxLayout;
		pbTextLayout->addWidget(pbTextLabel);
		pbTextLayout->addWidget(pbTextCombo);

	// Push Button border
	pbBorderLabel = new QLabel(tr("Button - Border:"));
	pbBorderCombo = new QComboBox;
	createColourPalette(pbBorderCombo);
	QHBoxLayout *pbBorderLayout = new QHBoxLayout;
		pbBorderLayout->addWidget(pbBorderLabel);
		pbBorderLayout->addWidget(pbBorderCombo);

	// Push Button (hover) background
	pbHoverBGLabel = new QLabel(tr("Button (Hover) - Background:"));
	pbHoverBGCombo = new QComboBox;
	createColourPalette(pbHoverBGCombo);
	QHBoxLayout *pbHoverBGLayout = new QHBoxLayout;
		pbHoverBGLayout->addWidget(pbHoverBGLabel);
		pbHoverBGLayout->addWidget(pbHoverBGCombo);

	// Push Button (pressed) background
	pbPressedBGLabel = new QLabel(tr("Button (Pressed) - Background:"));
	pbPressedBGCombo = new QComboBox;
	createColourPalette(pbPressedBGCombo);
	QHBoxLayout *pbPressedBGLayout = new QHBoxLayout;
		pbPressedBGLayout->addWidget(pbPressedBGLabel);
		pbPressedBGLayout->addWidget(pbPressedBGCombo);

	// List View background
	listViewBGLabel = new QLabel(tr("List View - Background:"));
	listViewBGCombo = new QComboBox;
	createColourPalette(listViewBGCombo);
	QHBoxLayout *listViewBGLayout = new QHBoxLayout;
		listViewBGLayout->addWidget(listViewBGLabel);
		listViewBGLayout->addWidget(listViewBGCombo);

	// List View text
	listViewTextLabel = new QLabel(tr("List View - Text:"));
	listViewTextCombo = new QComboBox;
	createColourPalette(listViewTextCombo);
	QHBoxLayout *listViewTextLayout = new QHBoxLayout;
		listViewTextLayout->addWidget(listViewTextLabel);
		listViewTextLayout->addWidget(listViewTextCombo);

	// Label text
	labelTextLabel = new QLabel(tr("Label - Text:"));
	labelTextCombo = new QComboBox;
	createColourPalette(labelTextCombo);
	QHBoxLayout *labelTextLayout = new QHBoxLayout;
		labelTextLayout->addWidget(labelTextLabel);
		labelTextLayout->addWidget(labelTextCombo);

	// Tool Tip background
	toolTipBGLabel = new QLabel(tr("Tool Tip - Background:"));
	toolTipBGCombo = new QComboBox;
	createColourPalette(toolTipBGCombo);
	QHBoxLayout *toolTipBGLayout = new QHBoxLayout;
		toolTipBGLayout->addWidget(toolTipBGLabel);
		toolTipBGLayout->addWidget(toolTipBGCombo);

	// Tool Tip text
	toolTipTextLabel = new QLabel(tr("Tool Tip - Text:"));
	toolTipTextCombo = new QComboBox;
	createColourPalette(toolTipTextCombo);
	QHBoxLayout *toolTipTextLayout = new QHBoxLayout;
		toolTipTextLayout->addWidget(toolTipTextLabel);
		toolTipTextLayout->addWidget(toolTipTextCombo);

	// Tool Tip border
	toolTipBorderLabel = new QLabel(tr("Tool Tip - Border:"));
	toolTipBorderCombo = new QComboBox;
	createColourPalette(toolTipBorderCombo);
	QHBoxLayout *toolTipBorderLayout = new QHBoxLayout;
		toolTipBorderLayout->addWidget(toolTipBorderLabel);
		toolTipBorderLayout->addWidget(toolTipBorderCombo);

	// Menu Bar background
	menuBarBGLabel = new QLabel(tr("Menu Bar - Background:"));
	menuBarBGCombo = new QComboBox;
	createColourPalette(menuBarBGCombo);
	QHBoxLayout *menuBarBGLayout = new QHBoxLayout;
		menuBarBGLayout->addWidget(menuBarBGLabel);
		menuBarBGLayout->addWidget(menuBarBGCombo);

	// Menu Bar text
	menuBarTextLabel = new QLabel(tr("Menu Bar - Text:"));
	menuBarTextCombo = new QComboBox;
	createColourPalette(menuBarTextCombo);
	QHBoxLayout *menuBarTextLayout = new QHBoxLayout;
		menuBarTextLayout->addWidget(menuBarTextLabel);
		menuBarTextLayout->addWidget(menuBarTextCombo);

	// Tab Bar background
	tbTabBGLabel = new QLabel(tr("Tab Bar - Background:"));
	tbTabBGCombo = new QComboBox;
	createColourPalette(tbTabBGCombo);
	QHBoxLayout *tbTabBGLayout = new QHBoxLayout;
		tbTabBGLayout->addWidget(tbTabBGLabel);
		tbTabBGLayout->addWidget(tbTabBGCombo);

	// Tab Bar Text
	tbTabTextLabel = new QLabel(tr("Tab Bar - Text:"));
	tbTabTextCombo = new QComboBox;
	createColourPalette(tbTabTextCombo);
	QHBoxLayout *tbTabTextLayout = new QHBoxLayout;
		tbTabTextLayout->addWidget(tbTabTextLabel);
		tbTabTextLayout->addWidget(tbTabTextCombo);

	// Tab Bar (selected) background
	tbTabSelectedBGLabel = new QLabel(tr("Tab Bar (Selected) - Background:"));
	tbTabSelectedBGCombo = new QComboBox;
	createColourPalette(tbTabSelectedBGCombo);
	QHBoxLayout *tbTabSelectedBGLayout = new QHBoxLayout;
		tbTabSelectedBGLayout->addWidget(tbTabSelectedBGLabel);
		tbTabSelectedBGLayout->addWidget(tbTabSelectedBGCombo);

	// Tab Bar (selected) text
	tbTabSelectedTextLabel = new QLabel(tr("Tab Bar (Selected) - Text:"));
	tbTabSelectedTextCombo = new QComboBox;
	createColourPalette(tbTabSelectedTextCombo);
	QHBoxLayout *tbTabSelectedTextLayout = new QHBoxLayout;
		tbTabSelectedTextLayout->addWidget(tbTabSelectedTextLabel);
		tbTabSelectedTextLayout->addWidget(tbTabSelectedTextCombo);

	// Tool Bar background
	toolBarBGLabel = new QLabel(tr("Tool Bar - Background:"));
	toolBarBGCombo = new QComboBox;
	createColourPalette(toolBarBGCombo);
	QHBoxLayout *toolBarBGLayout = new QHBoxLayout;
		toolBarBGLayout->addWidget(toolBarBGLabel);
		toolBarBGLayout->addWidget(toolBarBGCombo);

	// Check Box text
	checkBoxTextLabel = new QLabel(tr("Check Box - Text:"));
	checkBoxTextCombo = new QComboBox;
	createColourPalette(checkBoxTextCombo);
	QHBoxLayout *checkBoxTextLayout = new QHBoxLayout;
		checkBoxTextLayout->addWidget(checkBoxTextLabel);
		checkBoxTextLayout->addWidget(checkBoxTextCombo);

	// Group Box text
	groupBoxTextLabel = new QLabel(tr("Group Box - Text:"));
	groupBoxTextCombo = new QComboBox;
	createColourPalette(groupBoxTextCombo);
	QHBoxLayout *groupBoxTextLayout = new QHBoxLayout;
		groupBoxTextLayout->addWidget(groupBoxTextLabel);
		groupBoxTextLayout->addWidget(groupBoxTextCombo);

	// Group Box text
	msgBoxLabel = new QLabel(tr("Message Box - Background:"));
	msgBoxCombo = new QComboBox;
	createColourPalette(msgBoxCombo);
	QHBoxLayout *msgBoxLayout = new QHBoxLayout;
	msgBoxLayout->addWidget(msgBoxLabel);
	msgBoxLayout->addWidget(msgBoxCombo);

	// Add all the elements together
	QVBoxLayout *generalLayout = new QVBoxLayout;
		generalLayout->addLayout(nameLayout);
		generalLayout->addLayout(widgetBGLayout);
		generalLayout->addLayout(pbBGLayout);
		generalLayout->addLayout(pbTextLayout);
		generalLayout->addLayout(pbBorderLayout);
		generalLayout->addLayout(pbHoverBGLayout);
		generalLayout->addLayout(pbPressedBGLayout);
		generalLayout->addLayout(listViewBGLayout);
		generalLayout->addLayout(listViewTextLayout);
		generalLayout->addLayout(labelTextLayout);
		generalLayout->addLayout(toolTipBGLayout);
		generalLayout->addLayout(toolTipTextLayout);
		generalLayout->addLayout(toolTipBorderLayout);
		generalLayout->addLayout(menuBarBGLayout);
		generalLayout->addLayout(menuBarTextLayout);
		generalLayout->addLayout(tbTabBGLayout);
		generalLayout->addLayout(tbTabTextLayout);
		generalLayout->addLayout(tbTabSelectedBGLayout);
		generalLayout->addLayout(tbTabSelectedTextLayout);
		generalLayout->addLayout(toolBarBGLayout);
		generalLayout->addLayout(checkBoxTextLayout);
		generalLayout->addLayout(groupBoxTextLayout);	
		generalLayout->addLayout(msgBoxLayout);

	// Add the layout to the Editor widget
	ui.scrollAreaWidgetContents_2->setLayout(generalLayout);
	
	// Load the current colour style settings
	loadColourSettings(themeName);
}

void milxQtThemeEditorForm::discardTheme(QAbstractButton *button)
{
	// If "OK" button, do nothing
	if (!button->text().compare("OK")) return;

	// Check if "Discard" button was pressed while editing a theme
	if (!button->text().compare("Discard") && isEdit) {
		// Delete the theme
		QDir dir;
		dir.remove(QString(QDir::currentPath() + "/" + *themeName + ".qss"));
		prefForm->removeTheme(*themeName);
		MainWindow->update();
	}

	// Reset theme colours, close the editor
	prefForm->changeTheme(0);
	QDialog::reject();
}

void milxQtThemeEditorForm::accept()
{
	// Check if the theme is valid
	if ((prefForm->isTheme(themeNameInput->text()) && !isEdit) || 
		(prefForm->isTheme(themeNameInput->text()) && isEdit && themeName->compare(themeNameInput->text()))) {
		// Already a theme - highlight the input box
		themeNameInput->setStyleSheet(themeNameInput->styleSheet() + "QLineEdit{background: tomato;}");
	} else {
		// Save the theme
		saveColourSettings();
		QDialog::accept();
	}
}

void milxQtThemeEditorForm::createColourPalette(QComboBox *box)
{
	// Create a combobox item for each available colour
	for (int i = 0; i < QColor::colorNames().size(); i++)
	{
		box->addItem(QColor::colorNames()[i]);
		box->setItemData(i, QColor(QColor::colorNames()[i]), Qt::DecorationRole);
	}

	// Set the style for the combobox, to make it more viewable
	//box->setStyleSheet("background-color: " + QColor::colorNames()[0] + ";");
}

QString milxQtThemeEditorForm::getStyleColour(QTextStream *in, QString section)
{
	QString line; // A line in the qss file
	QString beginning; // For checking style sections
	QString colour; // The colour

	// Go through the file line by line
	while (!in->atEnd()) {
		line = in->readLine();
		beginning = line.section(" ", 0, 0);

		// Check if section reached
		if (!beginning.compare(section)) {
			break;
		}
	}

	// Get the colour used
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);

	// Return the style colour used
	return colour;
}

void milxQtThemeEditorForm::copyThemeFile(QTextStream *in, QString section, QStringList *data)
{
	QString line; // A line in the qss file
	QString beginning; // For checking style sections

	// Go through the file line by line
	while (!in->atEnd()) {
		line = in->readLine();
		beginning = line.section(" ", 0, 0);
		data->append(line);

		// Check if section reached
		if (!beginning.compare(section)) {
			break;
		}
	}
}

void milxQtThemeEditorForm::updateStyles(int index)
{
	QString colour;

	// Default widget background
	colour = widgetBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + ".QWidget{background-color: " + colour + ";}");
	
	// Push Button background
	colour = pbBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QPushButton{background-color: " + colour + ";}");
	
	// Push Button text
	colour = pbTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QPushButton{color: " + colour + ";}");

	// Push Button border
	colour = pbBorderCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QPushButton{border-color: " + colour + ";}");

	// Button (hover) background
	colour = pbHoverBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QPushButton:hover{background-color: " + colour + ";}");
	
	// Button (pressed) background
	colour = pbPressedBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QPushButton:pressed{background-color: " + colour + ";}");

	// List View background
	colour = listViewBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QListView{background-color: " + colour + ";}");

	// List View text
	colour = listViewTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QListView{color: " + colour + ";}");

	// Label Text
	colour = labelTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QLabel{color: " + colour + ";}");

	// Tool Tip background
	colour = toolTipBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QToolTip{background-color: " + colour + ";}");

	// Tool Tip text
	colour = toolTipTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QToolTip{color: " + colour + ";}");

	// Tool Tip border
	colour = toolTipBorderCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QToolTip{border-color: " + colour + ";}");

	// Menu Bar background
	colour = menuBarBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QMenuBar{background-color: " + colour + ";}");

	// Menu Bar text
	colour = menuBarTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QMenuBar{color: " + colour + ";}");

	// Tab Bar background
	colour = tbTabBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QTabBar::tab{background: " + colour + ";}");

	// Tab Bar text
	colour = tbTabTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QTabBar::tab{color: " + colour + ";}");

	// Tab Bar (selected) background
	colour = tbTabSelectedBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QTabBar::tab:selected{background: " + colour + ";}");

	// Tab Bar (selected) text
	colour = tbTabSelectedTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QTabBar::tab:selected{color: " + colour + ";}");

	// Tool Bar background
	colour = toolBarBGCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QToolBar{background: " + colour + ";}");

	// Check Box text
	colour = checkBoxTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QCheckBox{color: " + colour + ";}");

	// Group Box text
	colour = groupBoxTextCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QGroupBox::title{color: " + colour + ";}");

	// Message Box background
	colour = msgBoxCombo->currentText();
	qApp->setStyleSheet(qApp->styleSheet() + "QMessageBox{background: " + colour + ";}");
}

void milxQtThemeEditorForm::loadColourSettings(QString *theme)
{
	QString themeFile; // The theme filename

	// Open the qss file for the current theme
	if (!isEdit) {
		// A new theme - use "Light" as a template
		themeFile = QString(":/resources/styles/Light.qss");
	}
	else {
		// A custom theme
		themeFile = QString(QDir::currentPath() + "/" + *theme + ".qss");
	}
	QFile qss(themeFile);
	qss.open(QFile::ReadOnly);
	QTextStream *in = new QTextStream(&qss); // Text Stream
	QString colour;

	// Default widget background
	colour = getStyleColour(in, QString(".QWidget"));
	widgetBGCombo->setCurrentText(colour);
	
	// Push Button background
	colour = getStyleColour(in, QString("QPushButton"));
	pbBGCombo->setCurrentText(colour);

	// Push Button text
	QString line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	pbTextCombo->setCurrentText(colour);
	
	// Push Button border
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	pbBorderCombo->setCurrentText(colour);

	// Button (hover) background
	colour = getStyleColour(in, QString("QPushButton:hover"));
	pbHoverBGCombo->setCurrentText(colour);
	
	// Button (pressed) background
	colour = getStyleColour(in, QString("QPushButton:pressed"));
	pbPressedBGCombo->setCurrentText(colour);
	
	// List View background
	colour = getStyleColour(in, QString("QListView"));
	listViewBGCombo->setCurrentText(colour);
	
	// List View text
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	listViewTextCombo->setCurrentText(colour);

	// Label text
	colour = getStyleColour(in, QString("QLabel"));
	labelTextCombo->setCurrentText(colour);

	// Tool Tip background
	colour = getStyleColour(in, QString("QToolTip"));
	toolTipBGCombo->setCurrentText(colour);
	
	// Tool Tip text
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	toolTipTextCombo->setCurrentText(colour);

	// Tool Tip border
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	toolTipBorderCombo->setCurrentText(colour);

	// Menu Bar background
	colour = getStyleColour(in, QString("QMenuBar"));
	menuBarBGCombo->setCurrentText(colour);
	
	// Menu Bar text
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	menuBarTextCombo->setCurrentText(colour);

	// Tab Bar background
	colour = getStyleColour(in, QString("QTabBar::tab"));
	tbTabBGCombo->setCurrentText(colour);
		
	// Tab Bar Text
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	tbTabTextCombo->setCurrentText(colour);

	// Tab Bar (selected) background
	colour = getStyleColour(in, QString("QTabBar::tab:selected"));
	tbTabSelectedBGCombo->setCurrentText(colour);
	
	// Tab Bar (selected) text
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	tbTabSelectedTextCombo->setCurrentText(colour);
	
	// Tool Bar background
	colour = getStyleColour(in, QString("QToolBar"));
	toolBarBGCombo->setCurrentText(colour);

	// Check Box text
	colour = getStyleColour(in, QString("QCheckBox"));
	checkBoxTextCombo->setCurrentText(colour);

	// Group Box text
	colour = getStyleColour(in, QString("QGroupBox::title"));
	groupBoxTextCombo->setCurrentText(colour);

	// Message Box text
	colour = getStyleColour(in, QString("QMessageBox"));
	msgBoxCombo->setCurrentText(colour);

	// Close the style file
	qss.close();
}

void milxQtThemeEditorForm::saveColourSettings()
{
	/** Copy the layout for the "Light" theme **/
	QString themeFile = QString(":/resources/styles/Light.qss");
	QFile tmp(themeFile);
	tmp.open(QFile::ReadOnly);

	// Copy the layout, but use the colours in the new theme instead
	QTextStream *in = new QTextStream(&tmp); // Text Stream
	QStringList themeData = QStringList();

	// Default widget background
	copyThemeFile(in, QString(".QWidget"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + widgetBGCombo->currentText() + ";");

	// Push Button background
	copyThemeFile(in, QString("QPushButton"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + pbBGCombo->currentText() + ";");

	// Push Button text
	in->readLine(); // Eat next line
	themeData.append("    color: " + pbTextCombo->currentText() + ";");

	// Push Button border
	in->readLine(); // Eat next line
	themeData.append("    border-color: " + pbBorderCombo->currentText() + ";");

	// Button (hover) background
	copyThemeFile(in, QString("QPushButton:hover"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + pbHoverBGCombo->currentText() + ";");
	
	// Button (pressed) background
	copyThemeFile(in, QString("QPushButton:pressed"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + pbPressedBGCombo->currentText() + ";");

	// List View background
	copyThemeFile(in, QString("QListView"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + listViewBGCombo->currentText() + ";");

	// List View text
	in->readLine(); // Eat next line
	themeData.append("    color: " + listViewTextCombo->currentText() + ";");

	// Label text
	copyThemeFile(in, QString("QLabel"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    color: " + labelTextCombo->currentText() + ";");

	// Tool Tip background
	copyThemeFile(in, QString("QToolTip"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + toolTipBGCombo->currentText() + ";");

	// Tool Tip text
	in->readLine(); // Eat next line
	themeData.append("    color: " + toolTipTextCombo->currentText() + ";");

	// Tool Tip border
	in->readLine(); // Eat next line
	themeData.append("    border-color: " + toolTipBorderCombo->currentText() + ";");

	// Menu Bar background
	copyThemeFile(in, QString("QMenuBar"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + menuBarBGCombo->currentText() + ";");

	// Menu Bar text
	in->readLine(); // Eat next line
	themeData.append("    color: " + menuBarTextCombo->currentText() + ";");

	// Tab Bar background
	copyThemeFile(in, QString("QTabBar::tab"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background: " + tbTabBGCombo->currentText() + ";");

	// Tab Bar Text
	in->readLine();
	themeData.append("    color: " + tbTabTextCombo->currentText() + ";");
	
	// Tab Bar (selected)
	copyThemeFile(in, QString("QTabBar::tab:selected"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background: " + tbTabSelectedBGCombo->currentText() + ";");

	// Tab Bar (selected) text
	in->readLine(); // Eat next line
	themeData.append("    color: " + tbTabSelectedTextCombo->currentText() + ";");
	
	// Tool Bar background
	copyThemeFile(in, QString("QToolBar"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background: " + toolBarBGCombo->currentText() + ";");
	
	// Check Box text
	copyThemeFile(in, QString("QCheckBox"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    color: " + checkBoxTextCombo->currentText() + ";");

	// Group Box text
	copyThemeFile(in, QString("QGroupBox::title"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    color: " + groupBoxTextCombo->currentText() + ";");

	// Group Box text
	copyThemeFile(in, QString("QMessageBox"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background: " + msgBoxCombo->currentText() + ";\n}");

	// Close the read-in file
	tmp.close();

	/** Write the theme data to a qss file **/
	themeFile = QString(QDir::currentPath() + "/" + themeNameInput->text() + ".qss");
	QFile qss(themeFile);

	// Try to write the theme data
	if (qss.open(QFile::WriteOnly)) {
		QTextStream out(&qss); // Text Stream

		// Write the theme data
		for (int i = 0; i < themeData.size(); i++) {
			out << themeData.at(i) << '\n';
		}
	}
	else {
		// Something went wrong
		MainWindow->printDebug("Error opening output file.");
	}

	// Close the file
	qss.close();

	// Check if name update needed
	if (isEdit) {
		if (themeNameInput->text().compare(themeName)) {
			// Name needs updating
			QDir dir;
			dir.remove(QString(QDir::currentPath() + "/" + themeName + ".qss"));
			prefForm->removeTheme(*themeName);
		}
	}

	// Add the theme to the list of themes
	prefForm->addTheme(themeNameInput->text());
}

void milxQtThemeEditorForm::createConnections()
{
	connect(ui.buttonBox, SIGNAL(clicked(QAbstractButton*)), this, SLOT(discardTheme(QAbstractButton*)));
	connect(ui.buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(widgetBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(pbBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(pbTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(pbBorderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(pbHoverBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(pbPressedBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(listViewBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(listViewTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(labelTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(toolTipBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(toolTipTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(toolTipBorderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(menuBarBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(menuBarTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(tbTabBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(tbTabTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(tbTabSelectedBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(tbTabSelectedTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(toolBarBGCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(checkBoxTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(groupBoxTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(msgBoxCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
}
