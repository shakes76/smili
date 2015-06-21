#include "milxQtRegistrationWindow.h"
#include <QCheckBox>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>

milxQtRegistrationWindow::milxQtRegistrationWindow(milxQtMain *theParent) : QDialog(theParent)
{
	ui.setupUi(this);
	MainWindow = theParent;
	advancedOptionsWindow = new milxQtRegistrationAdvancedOptions(MainWindow);
	niftiReg = new milxQtRegistrationNifti(MainWindow);
	setWindowModality(Qt::ApplicationModal);
	setWindowTitle(tr("Registration option"));

	createConnections();
}

milxQtRegistrationWindow::~milxQtRegistrationWindow()
{
	delete advancedOptionsWindow;
	delete niftiReg;
}

void milxQtRegistrationWindow::setup(RegType regType)
{
	// Get the list of opened images
	getListOfHandledImages();
	if (imageList.size() < 2)
	{
        MainWindow->printError("At least two images need to be opened to perform a registration");
		this->reject();
	}

	// Set up the list of algorithms
	QStringList algoList;
    algoList << "FFD" << "Affine";
	this->ui.comboBoxAlgo->addItems(algoList);
	this->ui.comboBoxAlgo->setCurrentIndex(regType);
    algoComboChange(regType);


	// Set up the combo lists
	for (int i = 0; i < imageList.size(); i++)
	{
		QString itemName = QFileInfo(imageList[i]->getName()).fileName();

		// Add to the reference combobox
		this->ui.comboBoxRef->addItem(itemName);
	
		// Add to the image list
		QListWidgetItem* item = new QListWidgetItem(itemName);
		item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
		item->setCheckState(Qt::Unchecked);

		// If it's the first image by default it is disabled
		if (i == 0) { item->setFlags(item->flags() ^ Qt::ItemIsEnabled); }

		// If it's the second image by default it is checked
		if (i == 1) { item->setCheckState(Qt::Checked); }

		this->ui.listWidget->addItem(item);
	}

	// Set default selected items
	this->ui.comboBoxRef->setCurrentIndex(0);

	// By default we show the deformation field
	this->ui.checkBoxDeformationF->setChecked(true);
}

// Create connections with the user interface
void milxQtRegistrationWindow::createConnections()
{
    // Reference combo box changed
	connect(this->ui.comboBoxRef, SIGNAL(currentIndexChanged(int)), this, SLOT(referenceComboChange(int)));
    connect(this->ui.comboBoxAlgo, SIGNAL(currentIndexChanged(int)), this, SLOT(algoComboChange(int)));
	connect(this->niftiReg, SIGNAL(registrationFinished()), this, SLOT(registrationFinished()));
	connect(this->niftiReg, SIGNAL(cpp2defFinished()), this, SLOT(cpp2defFinished()));
	connect(this->ui.btnAdvancedOptions, SIGNAL(clicked()), this, SLOT(advancedOptionsClicked()));
}

// Advanced option button pushed
void milxQtRegistrationWindow::advancedOptionsClicked()
{
	advancedOptionsWindow->reset((RegType)this->ui.comboBoxAlgo->currentIndex());
	advancedOptionsWindow->show();
}


// Algo combo box changed
void milxQtRegistrationWindow::algoComboChange(int newIndex)
{
    if (newIndex == Affine)
    {
        this->ui.checkBoxDeformationF->setChecked(false);
        this->ui.checkBoxDeformationF->setDisabled(true);
    }
    else
    {
        this->ui.checkBoxDeformationF->setDisabled(false);
    }
}

// Reference combo box changed
void milxQtRegistrationWindow::referenceComboChange(int newIndex)
{
	// Enable all items and disable the new current reference item in the image to register list
	for (int row = 0; row < this->ui.listWidget->count(); row++)
	{
		QListWidgetItem *item = this->ui.listWidget->item(row);

		// If it's the previous reference we enable it
		if (!(item->flags() & Qt::ItemIsEnabled))
		{
			item->setFlags(item->flags() | Qt::ItemIsEnabled);
		}

		// If it's the newIndex we disable and uncheck it
		if (row == newIndex)
		{
			item->setCheckState(Qt::Unchecked);
			item->setFlags(item->flags() ^ Qt::ItemIsEnabled);
		}
	}

    // We check if at least one image is selected, if not we select the first one available
    QList<int> selected = getSelectedImages();
    if(selected.count() == 0)
    {
        for(int row = 0; row < this->ui.listWidget->count(); row++)
        {
            QListWidgetItem *item = this->ui.listWidget->item(row);

            if(item->flags() & Qt::ItemIsEnabled)
            {
                item->setCheckState(Qt::Checked);
                break;
            }
        }
    }
}

// Create a new temporary file, and copy its path in the buffer
bool milxQtRegistrationWindow::createTmpFile(char buffer[FILENAME_MAX + 1])
{
	QTemporaryFile tmpFile(QDir::tempPath() + "/smili_reg_XXXXXX.nii");

	tmpFile.setAutoRemove(false);
	if (!tmpFile.open()) {
		MainWindow->printError("Unable to create the temporary files");
		this->reject();
		return false;
	}
	tmpFile.close();
	strncpy(buffer, tmpFile.fileName().toLatin1().data(), FILENAME_MAX);

	return true;
}

// Btn Ok clicked
void milxQtRegistrationWindow::accept()
{
	// Get the algorithm of the registration
	RegType algo = (RegType)this->ui.comboBoxAlgo->currentIndex();

	// Get the reference image number
	int ref = this->ui.comboBoxRef->currentIndex();

	// Get the deformationField Option
	bool deformationField = this->ui.checkBoxDeformationF->isChecked();
	if (this->ui.checkBoxDeformationF->isEnabled() == false) { deformationField = false; }

	// Get the list of selected image to register
	QList<int> images = getSelectedImages();

	// Save the reference image in a temporary file
	milxQtFile writer;
	char refPath[FILENAME_MAX + 1];
	if (!createTmpFile(refPath)) { this->reject(); return; }
	if (!writer.saveImage(QString(refPath), imageList[ref]))
    {
        MainWindow->printError("Unable to save the reference image");
		this->reject();
        return;
    }

    // Add the images to the registration queue
    for(int i=0; i<images.count(); i++)
    {
		// Save the image in a temporary file
		char imgPath[FILENAME_MAX + 1];
		if (!createTmpFile(imgPath)) { this->reject(); return; }
		if (!writer.saveImage(QString(imgPath), imageList[images[i]]))
		{
			MainWindow->printError("Unable to save image");
			this->reject();
			return;
		}

		// Create the output path
		char outPath[FILENAME_MAX + 1];
		if (!createTmpFile(outPath)) { this->reject(); return; }

		// Fill the parameters of the registration
		RegistrationParams params;
		params.algo = algo;
		params.useCpp2Def = deformationField;
		
		if (algo == FFD)
		{
			strncpy(params.F3D.referenceName, refPath, FILENAME_MAX);
			strncpy(params.F3D.floatingName, imgPath, FILENAME_MAX);
			strncpy(params.F3D.outputWarpedName, outPath, FILENAME_MAX);

			// Create the cpp output file
			char cppPath[FILENAME_MAX + 1];
			if (!createTmpFile(cppPath)) { this->reject(); return; }
			strncpy(params.F3D.outputControlPointGridName, cppPath, FILENAME_MAX);

			params.F3D.maxiterationNumber = -1; 
			params.F3D.spacing[0] = 5;
			params.F3D.spacing[1] = 5;
			params.F3D.spacing[2] = 5;
			params.F3D.levelNumber = 3;
			params.F3D.levelToPerform = 3;
			params.F3D.noPyramid = false;
			params.F3D.useSym = false;

			/*
			params.F3D.maxiterationNumber = advancedOptionsWindow->getMaxItLevel();
			params.F3D.spacing[0] = advancedOptionsWindow->getSx();
			params.F3D.spacing[1] = advancedOptionsWindow->getSy();
			params.F3D.spacing[2] = advancedOptionsWindow->getSz();
			params.F3D.levelNumber = advancedOptionsWindow->getNbLevel();
			params.F3D.levelToPerform = advancedOptionsWindow->getFirstLevels();
			params.F3D.noPyramid = advancedOptionsWindow->getUsePyramidalApproach();
			params.F3D.useSym = advancedOptionsWindow->getUseSymmetricApproach();
			*/

			if (deformationField)
			{
				// Create the deformation field path
				char defPath[FILENAME_MAX + 1];
				if (!createTmpFile(defPath)) { this->reject(); return; }
				strncpy(params.cpp2Def.cpp2defOutputName, defPath, FILENAME_MAX);
				strncpy(params.cpp2Def.cpp2defInputName, cppPath, FILENAME_MAX);
				strncpy(params.cpp2Def.referenceImageName, refPath, FILENAME_MAX);
			} else {
				params.F3D.outputControlPointGridName[0] = '\0';
			}
		}
		else if (algo == Affine)
		{
			strncpy(params.Aladin.referenceName, refPath, FILENAME_MAX);
			strncpy(params.Aladin.floatingName, imgPath, FILENAME_MAX);
			strncpy(params.Aladin.outputWarpedName, outPath, FILENAME_MAX);

			params.Aladin.rigOnly = false;
			params.Aladin.affDirect = false;
			params.Aladin.maxiterationNumber = 5;
			params.Aladin.levelNumber = 3;
			params.Aladin.levelToPerform = 3;
			params.Aladin.useSym = false;
			params.Aladin.percentBlock = 50.0f;
		}

		regQueue.push_front(params);
    }

    // Disable all the fields
    this->ui.checkBoxDeformationF->setDisabled(true);
    this->ui.comboBoxAlgo->setDisabled(true);
    this->ui.comboBoxRef->setDisabled(true);
    this->ui.btnOk->setDisabled(true);
    this->ui.listWidget->setDisabled(true);
	this->ui.btnAdvancedOptions->setDisabled(true);

    // Call the registration function
    registration();

	// Hide the windows
	this->hide();
}

QList<int> milxQtRegistrationWindow::getSelectedImages()
{
    QList<int> images;

    for (int row = 0; row < this->ui.listWidget->count(); row++)
    {
        QListWidgetItem *item = this->ui.listWidget->item(row);

        // If it's the previous reference we enable it
        if ((item->flags() & Qt::ItemIsEnabled) && item->checkState())
        {
            images.append(row);
        }
    }

    return images;
}


// Btn Cancel clicked
void milxQtRegistrationWindow::reject()
{
	// Remove TMP files
	QStringList filenames;
	for (int i = 0; i < regQueue.size(); i++)
	{
		if (regQueue[i].algo == FFD)
		{
			filenames.push_back(QString(regQueue[i].F3D.referenceName));
			filenames.push_back(QString(regQueue[i].F3D.floatingName));
			filenames.push_back(QString(regQueue[i].F3D.outputWarpedName));
			filenames.push_back(QString(regQueue[i].F3D.outputControlPointGridName));

			if (regQueue[i].useCpp2Def)
			{
				filenames.push_back(QString(regQueue[i].cpp2Def.cpp2defOutputName));
			}
		}
		else if (regQueue[i].algo == Affine)
		{
			filenames.push_back(QString(regQueue[i].Aladin.referenceName));
			filenames.push_back(QString(regQueue[i].Aladin.floatingName));
			filenames.push_back(QString(regQueue[i].Aladin.outputWarpedName));
		}
	}

	for (int i = 0; i < filenames.size(); i++)
	{
		if (!QFile::remove(filenames[i]))
		{
			MainWindow->printError("Unable to remove temporary file: " + filenames[i]);
		}
	}

	// Clean all variables
	imageList.clear();
	regQueue.clear();
	filenames.clear();

	// Close and hide
	this->close();
	this->hide();
}

// Cpp2Def is done
void milxQtRegistrationWindow::cpp2defFinished()
{
	MainWindow->printInfo("cpp2def completed"); 
	
	// Open the deformation field
	QString defPath = QString(currentReg.cpp2Def.cpp2defOutputName);
	MainWindow->loadFile(defPath);
	QWidgetList windowsList = MainWindow->getListOfWindows();
	milxQtImage * window = qobject_cast<milxQtImage *>(windowsList[windowsList.size() - 1]);
	window->vectorField();
	window->close();

	// Remove the temporary files once loaded
	QString cppPath = QString(currentReg.cpp2Def.cpp2defInputName);
	QString refPath = QString(currentReg.cpp2Def.referenceImageName);
	if (!QFile::remove(refPath) || !QFile::remove(defPath) || !QFile::remove(cppPath))
	{
		MainWindow->printError("Unable to remove temporary deformation field files");
		this->reject();
	}

	// Perform the next registration
	registration();
}

// Registration is done
void milxQtRegistrationWindow::registrationFinished()
{
    MainWindow->printInfo("Registration completed");

	// Retrieve the path of the temporary files
	QString refPath, imgPath, outPath, cppPath;
	if (currentReg.algo == FFD) {
		refPath = QString(currentReg.F3D.referenceName);
		imgPath = QString(currentReg.F3D.floatingName);
		outPath = QString(currentReg.F3D.outputWarpedName);
		cppPath = QString(currentReg.F3D.outputControlPointGridName);
	}
	else {
		refPath = QString(currentReg.Aladin.referenceName);
		imgPath = QString(currentReg.Aladin.floatingName);
		outPath = QString(currentReg.Aladin.outputWarpedName);
		cppPath = "";
	}

	// Open the result of the registration
	MainWindow->loadFile(outPath);

	// Remove the temporary files once loaded (image and output path
	if (!QFile::remove(imgPath) || !QFile::remove(outPath))
	{
		MainWindow->printError("Unable to remove thes temporary files");
		this->reject();
	}

    // Calculate the deformation field is required
	if (currentReg.useCpp2Def)
    {
        MainWindow->printInfo("Calculation of the deformation field");		
		niftiReg->cpp2def_async(currentReg.cpp2Def);
    }
	else
	{
		// Remove the reference image, and the cpp file
		if (!QFile::remove(refPath))
		{
			MainWindow->printError("Unable to remove the reference path");
			this->reject();
		}

		// Remove the cpp file if F3D 
		if (currentReg.algo = FFD)
		{
			if (!QFile::remove(cppPath))
			{
				MainWindow->printError("Unable to remove the reference path");
				this->reject();
			}
		}

		// Perform the next registration
		registration();
	}
}


// Get the list of all the open images
void milxQtRegistrationWindow::getListOfHandledImages()
{
	QWidgetList windows;
	MainWindow->printInfo("Get the list of all handled images for the registration");
	windows = MainWindow->getListOfWindows();

	for (int i = 0; i < windows.size(); i++)
	{
		if (MainWindow->isImage(windows[i]))
		{
			milxQtImage * win = qobject_cast<milxQtImage *>(windows[i]);

			// We only handle float images
			if (win->isFloatingPointImage())
			{
				imageList.append(win);
			}
		}
	}
}

// Main registration function
void milxQtRegistrationWindow::registration()
{
	// If it was the last image
	if (regQueue.count() == 0)
	{
		// Enable all the fields
		this->ui.checkBoxDeformationF->setDisabled(false);
		this->ui.comboBoxAlgo->setDisabled(false);
		this->ui.comboBoxRef->setDisabled(false);
		this->ui.btnOk->setDisabled(false);
		this->ui.listWidget->setDisabled(false);
		this->ui.btnAdvancedOptions->setDisabled(false);

		return;
	}

    // Get information about the current registration
    currentReg = regQueue.last();
    regQueue.pop_back();

    // Display information
    MainWindow->printInfo("REGISTRATION STARTED");

    // We call the right program depending on the algorithm
    if(currentReg.algo == FFD)
	{ 
		niftiReg->f3d_async(currentReg.F3D);
	}
    else if(currentReg.algo == Affine)
	{
		niftiReg->aladin_async(currentReg.Aladin);
	}
    else
	{
		MainWindow->printError("INVALID REGISTRATION ALGORITHM");
	}

    MainWindow->printInfo("Please wait..");
}
