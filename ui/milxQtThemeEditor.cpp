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
	widgetLabel = new QLabel(tr("Background - Default:"));
	widgetCombo = new QComboBox;
	createColourPalette(widgetCombo);
	QHBoxLayout *widgetLayout = new QHBoxLayout;
		widgetLayout->addWidget(widgetLabel);
		widgetLayout->addWidget(widgetCombo);
	// Button background
	pbLabel = new QLabel(tr("Button - Default"));
	pbCombo = new QComboBox;
	createColourPalette(pbCombo);
	QHBoxLayout *pbLayout = new QHBoxLayout;
		pbLayout->addWidget(pbLabel);
		pbLayout->addWidget(pbCombo);
	// Button (hover) background
	pbHoverLabel = new QLabel(tr("Button - Hover:"));
	pbHoverCombo = new QComboBox;
	createColourPalette(pbHoverCombo);
	QHBoxLayout *pbHoverLayout = new QHBoxLayout;
		pbHoverLayout->addWidget(pbHoverLabel);
		pbHoverLayout->addWidget(pbHoverCombo);
	// Button (pressed) background
	pbPressedLabel = new QLabel(tr("Button - Pressed:"));
	pbPressedCombo = new QComboBox;
	createColourPalette(pbPressedCombo);
	QHBoxLayout *pbPressedLayout = new QHBoxLayout;
		pbPressedLayout->addWidget(pbPressedLabel);
		pbPressedLayout->addWidget(pbPressedCombo);

	// List View
	listViewLabel = new QLabel(tr("List View:"));
	listViewCombo = new QComboBox;
	createColourPalette(listViewCombo);
	QHBoxLayout *listViewLayout = new QHBoxLayout;
		listViewLayout->addWidget(listViewLabel);
		listViewLayout->addWidget(listViewCombo);

	// Tool Tip
	toolTipLabel = new QLabel(tr("Tool Tip Background:"));
	toolTipCombo = new QComboBox;
	createColourPalette(toolTipCombo);
	QHBoxLayout *toolTipLayout = new QHBoxLayout;
		toolTipLayout->addWidget(toolTipLabel);
		toolTipLayout->addWidget(toolTipCombo);

	// Menu Bar
	menuBarLabel = new QLabel(tr("Menu Bar:"));
	menuBarCombo = new QComboBox;
	createColourPalette(menuBarCombo);
	QHBoxLayout *menuBarLayout = new QHBoxLayout;
		menuBarLayout->addWidget(menuBarLabel);
		menuBarLayout->addWidget(menuBarCombo);

	// Tab Bar
	tbTabLabel = new QLabel(tr("Tab Bar - Background:"));
	tbTabCombo = new QComboBox;
	createColourPalette(tbTabCombo);
	QHBoxLayout *tbTabLayout = new QHBoxLayout;
		tbTabLayout->addWidget(tbTabLabel);
		tbTabLayout->addWidget(tbTabCombo);
	// Tab Bar Text
	tbTabTextLabel = new QLabel(tr("Tab Bar - Text:"));
	tbTabTextCombo = new QComboBox;
	createColourPalette(tbTabTextCombo);
	QHBoxLayout *tbTabTextLayout = new QHBoxLayout;
		tbTabTextLayout->addWidget(tbTabTextLabel);
		tbTabTextLayout->addWidget(tbTabTextCombo);
	// Tab Bar (selected)
	tbTabSelectedLabel = new QLabel(tr("Tab Bar - Selected:"));
	tbTabSelectedCombo = new QComboBox;
	createColourPalette(tbTabSelectedCombo);
	QHBoxLayout *tbTabSelectedLayout = new QHBoxLayout;
		tbTabSelectedLayout->addWidget(tbTabSelectedLabel);
		tbTabSelectedLayout->addWidget(tbTabSelectedCombo);
	// Tab Bar (selected) text
	tbTabSelectedTextLabel = new QLabel(tr("Tab Bar - Selected Text:"));
	tbTabSelectedTextCombo = new QComboBox;
	createColourPalette(tbTabSelectedTextCombo);
	QHBoxLayout *tbTabSelectedTextLayout = new QHBoxLayout;
		tbTabSelectedTextLayout->addWidget(tbTabSelectedTextLabel);
		tbTabSelectedTextLayout->addWidget(tbTabSelectedTextCombo);

	// Tool Bar
	toolBarLabel = new QLabel(tr("Tool Bar:"));
	toolBarCombo = new QComboBox;
	createColourPalette(toolBarCombo);
	QHBoxLayout *toolBarLayout = new QHBoxLayout;
		toolBarLayout->addWidget(toolBarLabel);
		toolBarLayout->addWidget(toolBarCombo);

	// Add all the elements together
	QVBoxLayout *generalLayout = new QVBoxLayout;
		generalLayout->addLayout(nameLayout);
		generalLayout->addLayout(widgetLayout);
		generalLayout->addLayout(pbLayout);
		generalLayout->addLayout(pbHoverLayout);
		generalLayout->addLayout(pbPressedLayout);
		generalLayout->addLayout(listViewLayout);
		generalLayout->addLayout(toolTipLayout);
		generalLayout->addLayout(menuBarLayout);
		generalLayout->addLayout(tbTabLayout);
		generalLayout->addLayout(tbTabTextLayout);
		generalLayout->addLayout(tbTabSelectedLayout);
		generalLayout->addLayout(tbTabSelectedTextLayout);
		generalLayout->addLayout(toolBarLayout);

	// Add the layout to the Editor widget
	ui.generalWidget->setLayout(generalLayout);

	// Load the current colour style settings
	loadColourSettings(themeName);
}

void milxQtThemeEditorForm::discardTheme(QAbstractButton *button)
{
	// Check if "Discard" button was pressed
	if (!button->text().compare("Discard")) {
		// Delete the theme (if editing), reset the theme, and close the editor
		if (isEdit) {
			QDir dir;
			dir.remove(QString(QDir::currentPath() + "/" + *themeName + ".qss"));
			prefForm->removeTheme(*themeName);
			MainWindow->update();
		}
		prefForm->changeTheme(0);
		QDialog::reject();
	}
	else if (!button->text().compare("Cancel")) { // "Cancel" button pressed
		// Reset theme colours, close the editor
		prefForm->changeTheme(0);
		QDialog::reject();
	}
}

void milxQtThemeEditorForm::accept()
{
	// Save the style and close the editor, if valid theme name
	if (prefForm->isTheme(themeNameInput->text())) {
		// Already a theme - highlight the input box
		MainWindow->printDebug("Already a Theme");
		themeNameInput->setStyleSheet("QLineEdit{background: tomato};" + themeNameInput->styleSheet());
	}
	else { // Not a theme - save the theme
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
	box->setStyleSheet("background-color: " + QColor::colorNames()[0] + ";");
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
	colour = widgetCombo->currentText();
	qApp->setStyleSheet(".QWidget{background-color: " + colour + ";}" + qApp->styleSheet());

	// Button background
	colour = pbCombo->currentText();
	qApp->setStyleSheet("QPushButton{background-color: " + colour + ";}" + qApp->styleSheet());
	
	// Button (hover) background
	colour = pbHoverCombo->currentText();
	qApp->setStyleSheet("QPushButton:hover{background-color: " + colour + ";}" + qApp->styleSheet());
	
	// Button (pressed) background
	colour = pbPressedCombo->currentText();
	qApp->setStyleSheet("QPushButton:pressed{background-color: " + colour + ";}" + qApp->styleSheet());

	// List View
	colour = listViewCombo->currentText();
	qApp->setStyleSheet("QListView{background-color: " + colour + ";}" + qApp->styleSheet());

	// Tool Tip
	colour = toolTipCombo->currentText();
	qApp->setStyleSheet("QToolTip{background-color: " + colour + ";}" + qApp->styleSheet());

	// Menu Bar
	colour = menuBarCombo->currentText();
	qApp->setStyleSheet("QMenuBar{background-color: " + colour + ";}" + qApp->styleSheet());

	// Tab Bar
	colour = tbTabCombo->currentText();
	qApp->setStyleSheet("QTabBar::tab{background: " + colour + ";}" + qApp->styleSheet());

	// Tab Bar Text
	colour = tbTabTextCombo->currentText();
	qApp->setStyleSheet("QTabBar::tab{color: " + colour + ";}" + qApp->styleSheet());

	// Tab Bar (selected)
	colour = tbTabSelectedCombo->currentText();
	qApp->setStyleSheet("QTabBar::tab:selected{background: " + colour + ";}" + qApp->styleSheet());

	// Tab Bar (selected) text
	colour = tbTabSelectedTextCombo->currentText();
	qApp->setStyleSheet("QTabBar::tab:selected{color: " + colour + ";}" + qApp->styleSheet());

	// Tool Bar
	colour = toolBarCombo->currentText();
	qApp->setStyleSheet("QToolBar{background-color: " + colour + ";}" + qApp->styleSheet());
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
	widgetCombo->setCurrentText(colour);
	
	// Button background
	colour = getStyleColour(in, QString("QPushButton"));
	pbCombo->setCurrentText(colour);
	
	// Button (hover) background
	colour = getStyleColour(in, QString("QPushButton:hover"));
	pbHoverCombo->setCurrentText(colour);
	
	// Button (pressed) background
	colour = getStyleColour(in, QString("QPushButton:pressed"));
	pbPressedCombo->setCurrentText(colour);
	
	// List View
	colour = getStyleColour(in, QString("QListView"));
	listViewCombo->setCurrentText(colour);
	
	// Tool Tip
	colour = getStyleColour(in, QString("QToolTip"));
	toolTipCombo->setCurrentText(colour);
	
	// Menu Bar
	colour = getStyleColour(in, QString("QMenuBar"));
	menuBarCombo->setCurrentText(colour);
	
	// Tab Bar
	colour = getStyleColour(in, QString("QTabBar::tab"));
	tbTabCombo->setCurrentText(colour);
		
	// Tab Bar Text
	QString line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	tbTabTextCombo->setCurrentText(colour);

	// Tab Bar (selected)
	colour = getStyleColour(in, QString("QTabBar::tab:selected"));
	tbTabSelectedCombo->setCurrentText(colour);
	
	// Tab Bar (selected) text
	line = in->readLine();
	colour = line.section(" ", 1, 1, QString::SectionSkipEmpty).section(";", 0, 0);
	tbTabSelectedTextCombo->setCurrentText(colour);
	
	// Tool Bar
	colour = getStyleColour(in, QString("QToolBar"));
	toolBarCombo->setCurrentText(colour);

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
	themeData.append("    background-color: " + widgetCombo->currentText() + ";");

	// Button background
	copyThemeFile(in, QString("QPushButton"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + pbCombo->currentText() + ";");

	// Button (hover) background
	copyThemeFile(in, QString("QPushButton:hover"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + pbHoverCombo->currentText() + ";");
	
	// Button (pressed) background
	copyThemeFile(in, QString("QPushButton:pressed"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + pbPressedCombo->currentText() + ";");

	// List View
	copyThemeFile(in, QString("QListView"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + listViewCombo->currentText() + ";");

	// Tool Tip
	copyThemeFile(in, QString("QToolTip"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + toolTipCombo->currentText() + ";");

	// Menu Bar
	copyThemeFile(in, QString("QMenuBar"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + menuBarCombo->currentText() + ";");

	// Tab Bar
	copyThemeFile(in, QString("QTabBar::tab"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background: " + tbTabCombo->currentText() + ";");

	// Tab Bar Text
	in->readLine();
	themeData.append("    color: " + tbTabTextCombo->currentText() + ";");
	
	// Tab Bar (selected)
	copyThemeFile(in, QString("QTabBar::tab:selected"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background: " + tbTabSelectedCombo->currentText() + ";");

	// Tab Bar (selected) text
	in->readLine();
	themeData.append("    color: " + tbTabSelectedTextCombo->currentText() + ";");
	
	// Tool Bar
	copyThemeFile(in, QString("QToolBar"), &themeData);
	in->readLine(); // Eat next line
	themeData.append("    background-color: " + toolBarCombo->currentText() + ";");
	
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
	connect(this->widgetCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->pbCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->pbHoverCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->pbPressedCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->listViewCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->toolTipCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->menuBarCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->tbTabCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->tbTabTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->tbTabSelectedCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->tbTabSelectedTextCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
	connect(this->toolBarCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateStyles(int)));
}