#include "milxQtRegistrationAdvancedOptions.h"

#include <QCheckBox>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>

milxQtRegistrationAdvancedOptions::milxQtRegistrationAdvancedOptions(milxQtMain *theParent) : QDialog(theParent)
{
	ui.setupUi(this);

	MainWindow = theParent;

	setWindowModality(Qt::ApplicationModal);
	setWindowTitle(tr("Registration Advanced Options"));

	createConnections();
	reset(FFD);
}

milxQtRegistrationAdvancedOptions::~milxQtRegistrationAdvancedOptions()
{

}

// Reset the interface
void milxQtRegistrationAdvancedOptions::reset(RegType algo)
{
	if (algo == FFD)
	{
		ui.titleSplinesOptions->setVisible(true);
		ui.descriptionSplinesOptions->setVisible(true);
		ui.sx->setVisible(true);
		ui.sy->setVisible(true);
		ui.sz->setVisible(true);
		ui.spinBoxSx->setVisible(true);
		ui.spinBoxSy->setVisible(true);
		ui.spinBoxSz->setVisible(true);
		ui.checkBoxPyramidalApproach->setVisible(true);

		ui.lblPctOfBlock->setVisible(false);
		ui.spinBoxPctBlock->setVisible(false);
		ui.checkBoxRigidRegistrationOnly->setVisible(false);
		ui.checkBoxAffineDirect->setVisible(false);

		ui.spinBoxSx->setValue(5);
		ui.spinBoxSy->setValue(5);
		ui.spinBoxSz->setValue(5);

		ui.spinBoxMaxItLevel->setValue(-1);
		ui.spinBoxNbLevel->setValue(3);
		ui.spinBoxFirstLevels->setValue(3);

		ui.checkBoxPyramidalApproach->setCheckState(Qt::Unchecked);
		ui.checkBoxSymmetricApproach->setCheckState(Qt::Unchecked);
	}

	if (algo == Affine)
	{
		ui.titleSplinesOptions->setVisible(false);
		ui.descriptionSplinesOptions->setVisible(false);
		ui.sx->setVisible(false);
		ui.sy->setVisible(false);
		ui.sz->setVisible(false);
		ui.spinBoxSx->setVisible(false);
		ui.spinBoxSy->setVisible(false);
		ui.spinBoxSz->setVisible(false);
		ui.checkBoxPyramidalApproach->setVisible(false);

		ui.lblPctOfBlock->setVisible(true);
		ui.spinBoxPctBlock->setVisible(true);
		ui.checkBoxRigidRegistrationOnly->setVisible(true);
		ui.checkBoxAffineDirect->setVisible(true);

		ui.spinBoxMaxItLevel->setValue(5);
		ui.spinBoxNbLevel->setValue(3);
		ui.spinBoxFirstLevels->setValue(3);
		ui.spinBoxPctBlock->setValue(50.0f);

		ui.checkBoxAffineDirect->setCheckState(Qt::Unchecked);
		ui.checkBoxRigidRegistrationOnly->setCheckState(Qt::Unchecked);
		ui.checkBoxSymmetricApproach->setCheckState(Qt::Unchecked);
	}
}

// Create connections with the user interface
void milxQtRegistrationAdvancedOptions::createConnections()
{

}

// Get input Sx
double milxQtRegistrationAdvancedOptions::getSx()
{
	return ui.spinBoxSx->value();
}

// Get input Sy
double milxQtRegistrationAdvancedOptions::getSy()
{
	return ui.spinBoxSy->value();
}

// Get input Sz
double milxQtRegistrationAdvancedOptions::getSz()
{
	return ui.spinBoxSz->value();
}

// Get input max it level
int milxQtRegistrationAdvancedOptions::getMaxItLevel()
{
	return ui.spinBoxMaxItLevel->value();
}

// Get input nb level
int milxQtRegistrationAdvancedOptions::getNbLevel()
{
	return ui.spinBoxNbLevel->value();
}

// Get input first levels
int milxQtRegistrationAdvancedOptions::getFirstLevels()
{
	return ui.spinBoxFirstLevels->value();
}

// Get input use pyramidal approach
bool milxQtRegistrationAdvancedOptions::getUsePyramidalApproach()
{
	return ui.checkBoxPyramidalApproach->isChecked();
}

// Get input ise symmetric approach
bool milxQtRegistrationAdvancedOptions::getUseSymmetricApproach()
{
	return ui.checkBoxSymmetricApproach->isChecked();
}

// Btn Ok clicked
void milxQtRegistrationAdvancedOptions::accept()
{
	// Close the window
	this->close();
	this->hide();
}

// Btn Cancel clicked
void milxQtRegistrationAdvancedOptions::reject()
{
    // Close the window
	this->close();
	this->hide();
}