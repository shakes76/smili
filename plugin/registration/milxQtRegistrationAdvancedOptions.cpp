#include "milxQtRegistrationAdvancedOptions.h"

#include <QCheckBox>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>

milxQtRegistrationAdvancedOptions::milxQtRegistrationAdvancedOptions(QDialog *theParent) : QDialog(theParent)
{
	ui.setupUi(this);

	setWindowModality(Qt::ApplicationModal);
	setWindowTitle(tr("Registration Advanced Options"));

	createConnections();
	reset(F3D);
}

milxQtRegistrationAdvancedOptions::~milxQtRegistrationAdvancedOptions()
{

}

// Reset the interface
void milxQtRegistrationAdvancedOptions::reset(RegType algo)
{
	QSettings settings("Smili", "Registration Plugin");
	currentAlgo = algo;


	if (algo == F3D)
	{
		ui.titleSplinesOptions->setVisible(true);
		ui.descriptionSplinesOptions->setVisible(true);
		ui.sx->setVisible(true);
		ui.sy->setVisible(true);
		ui.sz->setVisible(true);
		ui.spinBoxSx->setVisible(true);
		ui.spinBoxSy->setVisible(true);
		ui.spinBoxSz->setVisible(true);
		ui.checkBoxNopy->setVisible(true);
        ui.vSpacerSplines->changeSize(20, 20, QSizePolicy::Minimum, QSizePolicy::Minimum);

		ui.lblPctOfBlock->setVisible(false);
		ui.spinBoxPctBlock->setVisible(false);
		ui.checkBoxRigOnly->setVisible(false);
		ui.checkBoxAffDirect->setVisible(false);

		ui.spinBoxSx->setValue(settings.value("F3D/sx", -5).toFloat());
		ui.spinBoxSy->setValue(settings.value("F3D/sy", -5).toFloat());
		ui.spinBoxSz->setValue(settings.value("F3D/sz", -5).toFloat());

		ui.spinBoxMaxIt->setValue(settings.value("F3D/maxit", 300).toInt());
		ui.spinBoxLn->setValue(settings.value("F3D/ln", 3).toInt());
		ui.spinBoxLp->setValue(settings.value("F3D/lp", 3).toInt());

		ui.checkBoxNopy->setChecked(settings.value("F3D/nopy", false).toBool());
		ui.checkBoxSym->setChecked(settings.value("F3D/sym", false).toBool());
	}

	if (algo == Aladin)
	{
		ui.titleSplinesOptions->setVisible(false);
		ui.descriptionSplinesOptions->setVisible(false);
		ui.sx->setVisible(false);
		ui.sy->setVisible(false);
		ui.sz->setVisible(false);
		ui.spinBoxSx->setVisible(false);
		ui.spinBoxSy->setVisible(false);
		ui.spinBoxSz->setVisible(false);
		ui.checkBoxNopy->setVisible(false);
        ui.vSpacerSplines->changeSize(20, 20, QSizePolicy::Ignored, QSizePolicy::Ignored);
		
		ui.lblPctOfBlock->setVisible(true);
		ui.spinBoxPctBlock->setVisible(true);
		ui.checkBoxRigOnly->setVisible(true);
		ui.checkBoxAffDirect->setVisible(true);

		ui.spinBoxMaxIt->setValue(settings.value("aladin/maxit", 5).toInt());
		ui.spinBoxLn->setValue(settings.value("aladin/ln", 3).toInt());
		ui.spinBoxLp->setValue(settings.value("aladin/lp", 3).toInt());
		ui.spinBoxPctBlock->setValue(settings.value("aladin/%v", 50).toFloat());

		ui.checkBoxAffDirect->setChecked(settings.value("aladin/aF3Direct", false).toBool());
		ui.checkBoxRigOnly->setChecked(settings.value("aladin/rigOnly", false).toBool());
		ui.checkBoxSym->setChecked(settings.value("aladin/sym", false).toBool());
	}

	this->resize(10, 10);
}

// Create connections with the user interface
void milxQtRegistrationAdvancedOptions::createConnections()
{

}


ParamsF3D milxQtRegistrationAdvancedOptions::getParamsF3D()
{
	ParamsF3D params;

	if (currentAlgo != F3D)
	{
		this->reset(F3D);
	}

	// Get input Sx
	params.spacing[0] = ui.spinBoxSx->value();

	// Get input Sy
	params.spacing[1] = ui.spinBoxSy->value();

	// Get input Sz
	params.spacing[2] = ui.spinBoxSz->value();

	// Get input max it level
	params.maxit = ui.spinBoxMaxIt->value();

	// Get input nb level
	params.ln = ui.spinBoxLn->value();

	// Get input first levels
	params.lp = ui.spinBoxLp->value();

	// Get input use pyramidal approach
	params.nopy = ui.checkBoxNopy->isChecked();

	// Get input ise symmetric approach
	params.useSym = ui.checkBoxSym->isChecked();

	return params;
}


ParamsAladin milxQtRegistrationAdvancedOptions::getParamsAladin()
{
	ParamsAladin params;

	if (currentAlgo != Aladin)
	{
		this->reset(Aladin);
	}

	// Get input max it level
	params.maxit = ui.spinBoxMaxIt->value();

	// Get input nb level
	params.ln = ui.spinBoxLn->value();

	// Get input first levels
	params.lp = ui.spinBoxLp->value();

	// Get percentage of block
	params.percentBlock = ui.spinBoxPctBlock->value();

	// Aladin direct
	params.aF3Direct = ui.checkBoxAffDirect->isChecked();

	// Rigid only
	params.rigOnly = ui.checkBoxRigOnly->isChecked();

	// Get input ise symmetric approach
	params.useSym = ui.checkBoxSym->isChecked();

	return params;
}


// Btn Ok clicked
void milxQtRegistrationAdvancedOptions::accept()
{
	// Save parameters
	QSettings settings("Smili", "Registration Plugin");

	if (currentAlgo == F3D)
	{
		settings.setValue("F3D/sx", ui.spinBoxSx->value());
		settings.setValue("F3D/sy", ui.spinBoxSy->value());
		settings.setValue("F3D/sz", ui.spinBoxSz->value());

		settings.setValue("F3D/maxit", ui.spinBoxMaxIt->value());
		settings.setValue("F3D/ln", ui.spinBoxLn->value());
		settings.setValue("F3D/lp", ui.spinBoxLp->value());

		settings.setValue("F3D/nopy", ui.checkBoxNopy->isChecked());
		settings.setValue("F3D/sym", ui.checkBoxSym->isChecked());
	}
	else if (currentAlgo == Aladin)
	{
		settings.setValue("aladin/maxit", ui.spinBoxMaxIt->value());
		settings.setValue("aladin/ln", ui.spinBoxLn->value());
		settings.setValue("aladin/lp", ui.spinBoxLp->value());
		settings.setValue("aladin/%v", ui.spinBoxPctBlock->value());

		settings.setValue("aladin/aF3Direct", ui.checkBoxAffDirect->isChecked());
		settings.setValue("aladin/rigOnly", ui.checkBoxRigOnly->isChecked());
		settings.setValue("aladin/sym", ui.checkBoxSym->isChecked());
	}

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
