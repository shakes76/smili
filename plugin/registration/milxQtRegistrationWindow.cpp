#include "milxQtRegistrationWindow.h"
#include <QCheckBox>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>


/*
		Constructor
*/
milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget * theParent) : QDialog(theParent)
{
    MainWindow = qobject_cast<milxQtMain *>(theParent);
    advancedOptionsWindow = new milxQtRegistrationAdvancedOptions(this);
    niftiReg = new milxQtRegistrationNifti(this);
    workInProgress = false;
    computeAverage = false;
    openResults = false;
    initUI();
    createConnections();
}

/*
	Destructor
*/
milxQtRegistrationWindow::~milxQtRegistrationWindow()
{
    // Free the list of images
    qDeleteAll(images);
    images.clear();

    // Free the advanced options window
    delete(advancedOptionsWindow);

    // Destroy niftiReg
    delete(niftiReg);
}

// Initialise the user interface
void milxQtRegistrationWindow::initUI()
{
    // Setup window
    ui.setupUi(this);
    setWindowModality(Qt::ApplicationModal);
    setWindowTitle(tr("Registration option"));


    // Set up the list of algorithms
    QStringList algoList;
    algoList << "F3D" << "Aladin";
    this->ui.comboBoxAlgo->addItems(algoList);
}

// Change the algorithm selected (and change the interface accordingly)
void milxQtRegistrationWindow::setAlgo(RegType regType)
{
    // Reset interface
    this->ui.checkBoxSaveToDir->setChecked(false);
    this->ui.inputDirectoryBrowser->clear();

    // Set the current algo selected
    this->ui.comboBoxAlgo->setCurrentIndex(regType);

    // By default we show the deformation field for F3D
    if (regType == F3D)
    {
        this->ui.checkBoxDeformationF->setEnabled(true);
        this->ui.checkBoxDeformationF->setChecked(true);
    }
    // Otherwise we disable this field
    else
    {
        this->ui.checkBoxDeformationF->setEnabled(false);
        this->ui.checkBoxDeformationF->setChecked(false);
    }
    updateOpenImages();
    updateImageListCombo();
}


// Update the list of images in the combolist
void milxQtRegistrationWindow::updateImageListCombo()
{
    // If we have work in progress we don't update the image list until completed
    if (workInProgress) {
        return;
    }

    // We disconnect the event to avoid calling it while updateing the reference combo-box
    disconnect(this->ui.comboBoxRef, SIGNAL(currentIndexChanged(int)), this, SLOT(referenceComboChange(int)));

    // Update the parameters of images before displaying them
    this->updateParameters();

    // If we have two images
    if (images.size() >= 2)
    {
        // At least one image should be checked
        // We search if one image is already checked
        bool oneImageChecked = false;
        for (int i = 0; i < images.size(); i++)
        {
            if (images[i]->isChecked()) {
                oneImageChecked = true;
                break;
            }
        }

        // If no image checked, we check the first available
        if (oneImageChecked == false)
        {
            for (int i = 0; i < images.size(); i++)
            {
                if (!images[i]->isRef()) {
                    images[i]->setChecked(true);
                    break;
                }
            }
        }
    }

    // Clear the combos
    this->ui.comboBoxRef->clear();
    this->ui.listWidget->clear();

    // Add each image to the list
    for (int i = 0; i < images.size(); i++)
    {
        QString itemName = images[i]->getPath();

        // Add to the reference combobox
        this->ui.comboBoxRef->addItem(itemName);

        // Add to the image list
        QListWidgetItem* item = new QListWidgetItem(itemName);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);

        // Is the item checked ?
        if (images[i]->isChecked()) {
            item->setCheckState(Qt::Checked);
        } else {
            item->setCheckState(Qt::Unchecked);
        }

        // Is the item the reference (then disable it)
        if (images[i]->isRef()) {
            item->setFlags(item->flags() ^ Qt::ItemIsEnabled);
            this->ui.comboBoxRef->setCurrentIndex(i);
        }


        this->ui.listWidget->addItem(item);
    }


    connect(this->ui.comboBoxRef, SIGNAL(currentIndexChanged(int)), this, SLOT(referenceComboChange(int)));
}

// Update the list of open images that we can use for the registration
void milxQtRegistrationWindow::updateOpenImages()
{
    QWidgetList windows;
    windows = MainWindow->getListOfWindows();

    for (int i = 0; i < windows.size(); i++)
    {
        if (MainWindow->isImage(windows[i]))
        {
            milxQtImage * win = qobject_cast<milxQtImage *>(windows[i]);

            // We only handle float images
            if (win->isFloatingPointImage())
            {
                // If not in the list we add the image
                if (!isImageInList(win->getName()))
                {
                    addImage(new milxQtRegistration(this, win, MainWindow));
                }
            }
        }
    }
}

// Check if the image is in the list
bool milxQtRegistrationWindow::isImageInList(QString path)
{
    // We look is the image is already in the list
    bool inlist = false;

    for (int j = 0; j < images.size(); j++)
    {
        if (images[j]->getPath() == path)
        {
            inlist = true;
            break;
        }
    }
    return inlist;
}

// Disable the user interface
void milxQtRegistrationWindow::disableUI()
{
    // Disable all the fields
    this->ui.checkBoxDeformationF->setDisabled(true);
    this->ui.comboBoxAlgo->setDisabled(true);
    this->ui.comboBoxRef->setDisabled(true);
    this->ui.btnOk->setDisabled(true);
    this->ui.listWidget->setDisabled(true);
    this->ui.btnAdvancedOptions->setDisabled(true);

    // Disable checkboxes
    for (int row = 0; row < this->ui.listWidget->count(); row++)
    {
        QListWidgetItem *item = this->ui.listWidget->item(row);
        item->setFlags(item->flags() ^ Qt::ItemIsEnabled);
    }
}

// Enable the user interface
void milxQtRegistrationWindow::enableUI()
{
    // Disable all the fields
    this->ui.checkBoxDeformationF->setDisabled(false);
    this->ui.comboBoxAlgo->setDisabled(false);
    this->ui.comboBoxRef->setDisabled(false);
    this->ui.btnOk->setDisabled(false);
    this->ui.listWidget->setDisabled(false);
    this->ui.btnAdvancedOptions->setDisabled(false);

    // Disable checkboxes
    for (int row = 0; row < this->ui.listWidget->count(); row++)
    {
        QListWidgetItem *item = this->ui.listWidget->item(row);
        item->setFlags(item->flags() ^ Qt::ItemIsEnabled);
    }

    setAlgo((RegType) this->ui.comboBoxAlgo->currentIndex());
}

// get the parameters of F3D registration
ParamsF3D milxQtRegistrationWindow::getParamsF3D()
{
    ParamsF3D params = advancedOptionsWindow->getParamsF3D();
    params.cpp2Def = this->ui.checkBoxDeformationF->isChecked();
    return params;
}

// get the parameters of aladin registration
ParamsAladin milxQtRegistrationWindow::getParamsAladin()
{
    ParamsAladin params = advancedOptionsWindow->getParamsAladin();
    return params;
}

// update images parameters
void milxQtRegistrationWindow::updateParameters()
{
    // If no image to update, return
    if (images.size() == 0) return;

    // Get the registration type
    RegType type = (RegType)this->ui.comboBoxAlgo->currentIndex();

    // Get the reference image
    int refIndex = this->ui.comboBoxRef->currentIndex();
    milxQtRegistration * ref = NULL;
    if (refIndex >= 0 && refIndex < images.size())
    {
        ref = images[refIndex];
        ref->setIsRef(true);
    }

    // Get the output folder
    QString outFolder = this->ui.inputDirectoryBrowser->text();
    if (!this->ui.checkBoxSaveToDir->isChecked() || !QDir(outFolder).exists())
    {
        outFolder = "";
    }

    // Do we need to open the results
    openResults = this->ui.checkBoxOpenResults->isChecked();

    // Update the parameters of the registration for each image
    for (int i = 0; i < this->ui.listWidget->count(); i++)
    {
        // We look if the image is checked
        QListWidgetItem *item = this->ui.listWidget->item(i);
        if ((item->flags() & Qt::ItemIsEnabled) && item->checkState())
        {
            images[i]->setChecked(true);
        }
        else
        {
            images[i]->setChecked(false);
        }

        // We set the reference image
        images[i]->setReference(ref);

        // We set the output folder
        images[i]->setOutputFolder(outFolder);

        // We look the type of the registration
        images[i]->setRegType(type);

        // We set if we need to open the results of the registration
        images[i]->setOpenResults(openResults);

        // We update the parameters of the registration
        if (type == Aladin)
        {
            images[i]->setParams(getParamsAladin());
        }
        else if (type == F3D)
        {
            images[i]->setParams(getParamsF3D());
        }
    }
}


// Create connections with the user interface
void milxQtRegistrationWindow::createConnections()
{
    connect(this->ui.comboBoxRef, SIGNAL(currentIndexChanged(int)), this, SLOT(referenceComboChange(int)));
    connect(this->ui.comboBoxAlgo, SIGNAL(currentIndexChanged(int)), this, SLOT(algoComboChange(int)));
    connect(this->ui.btnAdvancedOptions, SIGNAL(clicked()), this, SLOT(advancedOptionsClicked()));
    connect(this->ui.btnAddImage, SIGNAL(clicked()), this, SLOT(addImageClicked()));
    connect(this->ui.btnSelectAll, SIGNAL(clicked()), this, SLOT(selectAllClicked()));
    connect(this->ui.btnUnselectAll, SIGNAL(clicked()), this, SLOT(unselectAllClicked()));
    connect(this->ui.btnBrowse, SIGNAL(clicked()), this, SLOT(browseBtnClicked()));
    connect(this->ui.clearList, SIGNAL(clicked()), this, SLOT(clearList()));
    connect(this->niftiReg, SIGNAL(averageCompleted()), this, SLOT(averageComputed()));
}

// Add a new image to the list
void milxQtRegistrationWindow::addImage(milxQtRegistration * image)
{
    if (images.size() == 0)
    {
        image->setIsRef(true);
    }

    images.append(image);
    connect(image, SIGNAL(done()), this, SLOT(regComplete()));
}


/*
					SLOTS
*/

// The average (Atlas) has been computed
void milxQtRegistrationWindow::averageComputed()
{
    // Do we need to open the results
    if (openResults)
    {
        MainWindow->loadFile(this->atlasPath);
    }

    // All the work has been completed
    workCompleted();
}


// Button Clear List clicked
void milxQtRegistrationWindow::clearList()
{
    // Free the list of images
    qDeleteAll(images);
    images.clear();

    // Update list
    this->updateImageListCombo();
}


// Button add image clicked
void milxQtRegistrationWindow::addImageClicked()
{
    QFileDialog fileOpener;
    fileOpener.setFileMode(QFileDialog::ExistingFiles);
    QStringList filenames = fileOpener.getOpenFileNames(this, tr("Select Files"));

    for (int i = 0; i < filenames.size(); i++)
    {
        if (!isImageInList(filenames[i]))
        {
            addImage(new milxQtRegistration(this, filenames[i], MainWindow));
        }
    }

    updateImageListCombo();

}

// Button select all clicked
void milxQtRegistrationWindow::selectAllClicked()
{
    for (int row = 0; row < this->ui.listWidget->count(); row++)
    {
        QListWidgetItem *item = this->ui.listWidget->item(row);

        if ((item->flags() & Qt::ItemIsEnabled))
        {
            item->setCheckState(Qt::Checked);
        }
    }
}

// Button unselect all clicked
void milxQtRegistrationWindow::unselectAllClicked()
{
    for (int row = 0; row < this->ui.listWidget->count(); row++)
    {
        QListWidgetItem *item = this->ui.listWidget->item(row);

        if ((item->flags() & Qt::ItemIsEnabled))
        {
            item->setCheckState(Qt::Unchecked);
        }
    }
}

// Button browse clicked
void milxQtRegistrationWindow::browseBtnClicked()
{
    QFileDialog dialog;
    dialog.setFileMode(QFileDialog::Directory);
    dialog.setOption(QFileDialog::ShowDirsOnly);
    QString folder = dialog.getExistingDirectory(this, tr("Select an output folder"));
    if (!folder.isEmpty())
    {
        this->ui.inputDirectoryBrowser->setText(folder);
    }
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
        setAlgo(Aladin);
    }
    else
    {
        setAlgo(F3D);
    }
}

// Reference combo box changed
void milxQtRegistrationWindow::referenceComboChange(int newIndex)
{
    // Set the current image as reference
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->setIsRef(false);

        if (i == newIndex)
        {
            images[i]->setIsRef(true);
        }
    }

    // We update the image list
    updateImageListCombo();
}

// Btn Ok clicked
void milxQtRegistrationWindow::accept()
{
    // Check if we have at least two images
    if (images.size() < 2) {
        // Message box: Error we need at least two images
        QMessageBox msgBox(QMessageBox::Critical, "Registration", "Error you need at least two images (a reference and and image to register).", QMessageBox::NoButton);
        msgBox.exec();
        return;
    }

    // Check if we have at least one image checked
    bool oneImageChecked = false;
    for (int i = 0; i < images.size(); i++)
    {
        if (images[i]->isChecked() == true) {
            oneImageChecked = true;
            break;
        }
    }

    if (!oneImageChecked)
    {
        // Message box: Error one image must be checked
        QMessageBox msgBox(QMessageBox::Critical, "Registration", "Error, at least one image should be checked", QMessageBox::NoButton);
        msgBox.exec();
        return;
    }

    // Update the parameters
    this->updateParameters();

    // Disable the form
    this->disableUI();

    // Hide the windows
    this->hide();

    // Do we have to compute the average of the registered image
    computeAverage = this->ui.checkBoxCreateAtlas->isChecked();

    // Perform the registration
    this->performRegistrations();
}

// Btn Cancel clicked
void milxQtRegistrationWindow::reject()
{
    // Close and hide
    this->close();
    this->hide();
}

// A registration has been completed
void milxQtRegistrationWindow::regComplete()
{
    performRegistrations();
}

// Perform the next registration
void milxQtRegistrationWindow::performRegistrations()
{
    int i;

    workInProgress = true;
    int test = images.size();
    for (i = 0; i < images.size(); i++)
    {
        if (images[i]->isChecked())
        {
            images[i]->startRegistration();
            break;
        }
    }


    // If all the registration are done
    if (i == images.size())
    {
        // If we have to compute the average
        if (computeAverage)
        {
            computeAtlas();
        }
        else
        {
            // The work is completed
            workCompleted();
        }
    }
}

// compute the average of all the registrations
void milxQtRegistrationWindow::computeAtlas()
{
    // Get the list of images and find the reference image
    QStringList filenames;
    milxQtRegistration * ref;

    for (int i = 0; i < images.size(); i++)
    {
        if (images[i]->isWorkDone())
        {
            filenames.append(images[i]->getOutputPath());
        }

        if (images[i]->isRef())
        {
            ref = images[i];
        }
    }

    // Create the outputfile for the atlas (the output filename will contain the name of the reference)
    atlasPath = ref->createAtlasFile();

    // Start the computation
    niftiReg->average_async(atlasPath, filenames);
}

// everything has been completed
void milxQtRegistrationWindow::workCompleted()
{
    workInProgress = false;
    this->enableUI();

    // Message box: registration completed
    QMessageBox msgBox(QMessageBox::Information, "Registration", "Registration completed !", QMessageBox::NoButton);
    msgBox.exec();
}
