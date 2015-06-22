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
	workInProcess = false;
	createConnections();
}

milxQtRegistrationWindow::~milxQtRegistrationWindow()
{
	delete advancedOptionsWindow;
	delete niftiReg;
}

void milxQtRegistrationWindow::setup(RegType regType)
{
	// If the work is in process, the interface stays disable
	if (workInProcess) { return; }

	// Reset interface
	this->ui.comboBoxAlgo->clear();
	this->ui.comboBoxRef->clear();
	this->ui.listWidget->clear();

	// Get the list of opened images
	getListOfHandledImages();
	if (imageList.size() < 2)
	{
        MainWindow->printError("At least two images need to be opened to perform a registration");
		this->reject();
	}

	// Set up the list of algorithms
	QStringList algoList;
    algoList << "F3D" << "Aladin";
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

	if (regType == F3D)
	{
		// By default we show the deformation field for F3D
		this->ui.checkBoxDeformationF->setChecked(true);
	}
	else
	{
		this->ui.checkBoxDeformationF->setChecked(false);
	}

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
    if (newIndex == Aladin)
    {
        this->ui.checkBoxDeformationF->setChecked(false);
        this->ui.checkBoxDeformationF->setDisabled(true);
    }
    else
    {
		if (!workInProcess)
		{
			this->ui.checkBoxDeformationF->setDisabled(false);
		}
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
		
		if (algo == F3D)
		{
			params.F3D = advancedOptionsWindow->getParamsF3D();

			strncpy(params.F3D.referenceName, refPath, FILENAME_MAX);
			strncpy(params.F3D.floatingName, imgPath, FILENAME_MAX);
			strncpy(params.F3D.outputWarpedName, outPath, FILENAME_MAX);

			// Create the cpp output file
			char cppPath[FILENAME_MAX + 1];
			if (!createTmpFile(cppPath)) { this->reject(); return; }
			strncpy(params.F3D.outputControlPointGridName, cppPath, FILENAME_MAX);

			// Create file if deformation field
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
		else if (algo == Aladin)
		{
			params.Aladin = advancedOptionsWindow->getParamsAladin();

			strncpy(params.Aladin.referenceName, refPath, FILENAME_MAX);
			strncpy(params.Aladin.floatingName, imgPath, FILENAME_MAX);
			strncpy(params.Aladin.outputResultName, outPath, FILENAME_MAX);
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


// Remove files created during a registration
void milxQtRegistrationWindow::removeFiles(RegistrationParams reg)
{
	QStringList filenames;

	// Create the list of file to delete depending on the algorithm
	if (reg.algo == F3D)
	{
		filenames.push_back(QString(reg.F3D.referenceName));
		filenames.push_back(QString(reg.F3D.floatingName));
		filenames.push_back(QString(reg.F3D.outputWarpedName));
		filenames.push_back(QString(reg.F3D.outputControlPointGridName));

		if (reg.useCpp2Def)
		{
			filenames.push_back(QString(reg.cpp2Def.cpp2defOutputName));
		}
	}
	else if (reg.algo == Aladin)
	{
		filenames.push_back(QString(reg.Aladin.referenceName));
		filenames.push_back(QString(reg.Aladin.floatingName));
		filenames.push_back(QString(reg.Aladin.outputResultName));
	}


	// Remove all files in the list
	for (int i = 0; i < filenames.size(); i++)
	{
		if (!QFile::remove(filenames[i]))
		{
			MainWindow->printError("Unable to remove temporary file: " + filenames[i]);
		}
	}
}

// Btn Cancel clicked
void milxQtRegistrationWindow::reject()
{
	// Remove TMP files
	QStringList filenames;
	for (int i = 0; i < regQueue.size(); i++)
	{
		removeFiles(regQueue[i]);
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

	// Remove the temporary files
	removeFiles(currentReg);

	// Perform the next registration
	registration();
}

// Registration is done
void milxQtRegistrationWindow::registrationFinished()
{
    MainWindow->printInfo("Registration completed");

	// Retrieve the path of the temporary files
	QString outPath;
	if (currentReg.algo == F3D) {
		outPath = QString(currentReg.F3D.outputWarpedName);
	}
	else {
		outPath = QString(currentReg.Aladin.outputResultName);
	}

	// Open the result of the registration
	MainWindow->loadFile(outPath);

    // Calculate the deformation field is required
	if (currentReg.useCpp2Def)
    {
        MainWindow->printInfo("Calculation of the deformation field");		
		niftiReg->cpp2def_async(currentReg.cpp2Def);
    }
	else
	{
		// Remove the temporary files
		removeFiles(currentReg);

		// Perform the next registration
		registration();
	}
}


// Get the list of all the open images
void milxQtRegistrationWindow::getListOfHandledImages()
{
	QWidgetList windows;
	MainWindow->printInfo("Get the list of all handled images for the registration");
	imageList.clear();
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

		// Work in process
		workInProcess = false;

		return;
	}

    // Get information about the current registration
    currentReg = regQueue.last();
    regQueue.pop_back();

    // Display information
    MainWindow->printInfo("REGISTRATION STARTED");
	workInProcess = true;

    // We call the right program depending on the algorithm
    if(currentReg.algo == F3D)
	{ 
		niftiReg->f3d_async(currentReg.F3D);
	}
    else if(currentReg.algo == Aladin)
	{
		niftiReg->aladin_async(currentReg.Aladin);
	}
    else
	{
		MainWindow->printError("INVALID REGISTRATION ALGORITHM");
	}

    MainWindow->printInfo("Please wait..");
}
