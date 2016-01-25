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
    regAlgos = new milxQtRegistrationAlgos(this);
    workInProgress = false;
    computeAverage = false;
    openResults = false;
    initUI();
    createConnections();

#ifdef USE_ELASTIX
    // This needs to be setup once for elastix to work properly
    int ret = elx::xoutSetup("", false, false);
    if (ret)
    {
        QMessageBox msgBox(QMessageBox::Critical, "Elastix Initialisation Error", "Error while initializing Elastix logging system. Elastix registration might fail.", QMessageBox::NoButton);
        msgBox.exec();
    }
#endif
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

    // Destroy regAlgos
    delete(regAlgos);
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
    algoList << "Affine (ITK)" << "Demon (ITK)";

#ifdef USE_NIFTI_REG
    algoList << "F3D (Nifti)" << "Aladin (Nifti)";
#endif

#ifdef USE_ELASTIX
    algoList << "Affine (Elastix)" << "BSpline (Elastix)";
#endif

    this->ui.comboBoxAlgo->addItems(algoList);

    // If we don't have nifti we can't compute the atlas and we can't compute the deformation field
#ifndef USE_NIFTI_REG
    this->ui.checkBoxCreateAtlas->setVisible(false);
    this->ui.checkBoxDeformationF->setVisible(false);
    this->ui.checkBoxSimilarities->setVisible(false);
#endif
}

// Return the current algorithm selected in the combo box
RegType milxQtRegistrationWindow::getCurrentAlgo()
{
    if (ui.comboBoxAlgo->currentText().compare("Affine (ITK)") == 0) return AffineItk;
    else if (ui.comboBoxAlgo->currentText().compare("Demon (ITK)") == 0) return DemonItk;
    else if (ui.comboBoxAlgo->currentText().compare("F3D (Nifti)") == 0) return F3DNifti;
    else if (ui.comboBoxAlgo->currentText().compare("Aladin (Nifti)") == 0) return AladinNifti;
    else if (ui.comboBoxAlgo->currentText().compare("Affine (Elastix)") == 0) return ElastixAffine;
    else if (ui.comboBoxAlgo->currentText().compare("BSpline (Elastix)") == 0) return ElastixBSpline;

    return None;
}

// Change the algorithm selected (and change the interface accordingly)
void milxQtRegistrationWindow::setAlgo(RegType regType)
{
    // Reset interface
    // Set the output directory
    if (this->ui.inputDirectoryBrowser->text() == "")
    {
        this->ui.inputDirectoryBrowser->setText(getDefaultOutputFolder());
    }

    // Set the current algo selected
    this->ui.comboBoxAlgo->setCurrentIndex(regType);

    // With ITK we don't have advanced options (hide the button)
    if (regType == AffineItk || regType == DemonItk)
    {
        this->ui.btnAdvancedOptions->setVisible(false);
    }
    else
    {
        this->ui.btnAdvancedOptions->setVisible(true);
    }

    // By default we don't show the deformation field for Nifti F3D
    if (regType == F3DNifti)
    {
        this->ui.checkBoxDeformationF->setVisible(true);
        this->ui.checkBoxDeformationF->setEnabled(true);
        this->ui.checkBoxDeformationF->setChecked(false);
    }
    // Otherwise we disable this field
    else
    {
        this->ui.checkBoxDeformationF->setVisible(false);
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

    MainWindow->initialiseWindowTraversal();
    for (int i = 0; i < MainWindow->getNumberOfWindows(); i ++)
    {
        milxQtImage *img = MainWindow->nextImage();
        // We only handle float images
        if (img->isFloatingPointImage())
        {
            // If the image is not in the list we add the image
            if (!isImageInList(img->getName()))
                addImage(new milxQtRegistration(this, img, MainWindow));
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
    this->ui.checkBoxCreateAtlas->setDisabled(true);
    this->ui.checkBoxOpenResults->setDisabled(true);
    this->ui.btnAddImage->setDisabled(true);
    this->ui.btnBrowse->setDisabled(true);
    this->ui.btnSelectAll->setDisabled(true);
    this->ui.btnUnselectAll->setDisabled(true);
    this->ui.clearList->setDisabled(true);
    this->ui.checkBoxSimilarities->setDisabled(true);

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
    this->ui.checkBoxCreateAtlas->setDisabled(false);
    this->ui.checkBoxOpenResults->setDisabled(false);
    this->ui.btnAddImage->setDisabled(false);
    this->ui.btnBrowse->setDisabled(false);
    this->ui.btnSelectAll->setDisabled(false);
    this->ui.btnUnselectAll->setDisabled(false);
    this->ui.clearList->setDisabled(false);
    this->ui.checkBoxSimilarities->setDisabled(false);

    // Disable checkboxes
    for (int row = 0; row < this->ui.listWidget->count(); row++)
    {
        QListWidgetItem *item = this->ui.listWidget->item(row);
        item->setFlags(item->flags() ^ Qt::ItemIsEnabled);
    }

    setAlgo(this->getCurrentAlgo());
}

// get the parameters of Itk Affine registration
milxQtRegistrationParams milxQtRegistrationWindow::getParamsAffineItk()
{
    milxQtRegistrationParams params;
    return params;
}


// get the parameters of Itk Demon registration
milxQtRegistrationParams milxQtRegistrationWindow::getParamsDemonItk()
{
    milxQtRegistrationParams params;
    return params;
}

// get the parameters of Nifti F3D registration
milxQtRegistrationParams milxQtRegistrationWindow::getParamsF3DNifti()
{
    milxQtRegistrationParams params = advancedOptionsWindow->getParamsF3DNifti();
    params.cpp2Def = this->ui.checkBoxDeformationF->isChecked();
    return params;
}

// get the parameters of a Nifti Aladin registration
milxQtRegistrationParams milxQtRegistrationWindow::getParamsAladinNifti()
{
    milxQtRegistrationParams params = advancedOptionsWindow->getParamsAladinNifti();
    return params;
}

// get the parameters of Elastix Affine registration
milxQtRegistrationParams milxQtRegistrationWindow::getParamsElastixAffine()
{
    milxQtRegistrationParams params = advancedOptionsWindow->getParamsElastixAffine();
    return params;
}

// get the parameters of Elastix BSpline registration
milxQtRegistrationParams milxQtRegistrationWindow::getParamsElastixBSpline()
{
    milxQtRegistrationParams params = advancedOptionsWindow->getParamsElastixBSpline();
    return params;
}

QString milxQtRegistrationWindow::getDefaultOutputFolder()
{
    QString outputFoldPath;

    // we try to find the reference image
    int indexRef;
    for (indexRef = 0; indexRef < images.size(); indexRef++)
    {
        if (images[indexRef]->isRef() == true)
        {
            break;
        }
    }

    // If we have images opened
    // the default path is the reference image path + regoutput
    if (images.size() != 0)
    {
        outputFoldPath = QDir(QFileInfo(images[indexRef]->getPath()).absolutePath() + "/regoutput").absolutePath();
    }
    else
    {
        // Otherwise If no image are opened
        // By default it will create a folder on the desktop
        outputFoldPath = QDir(QDesktopServices::storageLocation(QDesktopServices::DesktopLocation) + "/regoutput").absolutePath();
    }

    return outputFoldPath;
}

// update images parameters
void milxQtRegistrationWindow::updateParameters()
{
    // If no image to update, return
    if (images.size() == 0) return;

    // Get the registration type
    RegType type = this->getCurrentAlgo();

    // Get the reference image
    int refIndex = this->ui.comboBoxRef->currentIndex();
    milxQtRegistration * ref = NULL;
    if (refIndex >= 0 && refIndex < images.size())
    {
        ref = images[refIndex];
        ref->setIsRef(true);
    }

    // If we have to compute similarities
    bool computeSimilarities = false;
    if (this->ui.checkBoxSimilarities->isChecked())
    {
        computeSimilarities = true;
    }

    // Get the output folder
    QString outFolder = this->ui.inputDirectoryBrowser->text();
    if (!QDir(outFolder).exists())
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
        images[i]->setOutputFolder(QDir(outFolder).absolutePath());

        // We look the type of the registration
        images[i]->setRegType(type);

        // We set if we need to open the results of the registration
        images[i]->setOpenResults(openResults);

        // We update the parameters of the registration
        if (type == AffineItk)
        {
            images[i]->setParams(getParamsAffineItk());
        }
        else if (type == DemonItk)
        {
            images[i]->setParams(getParamsDemonItk());
        }
        else if (type == AladinNifti)
        {
            images[i]->setParams(getParamsAladinNifti());
        }
        else if (type == F3DNifti)
        {
            images[i]->setParams(getParamsF3DNifti());
        }
        else if (type == ElastixAffine)
        {
            images[i]->setParams(getParamsElastixAffine());
        }
        else if (type == ElastixBSpline)
        {
            images[i]->setParams(getParamsElastixBSpline());
        }

        // We set if we have to compute the similarities
        images[i]->setComputeSimilarities(computeSimilarities);
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

#ifdef USE_NIFTI_REG
    connect(this->regAlgos, SIGNAL(averageCompleted()), this, SLOT(averageComputed()));
#endif

    connect(this->regAlgos, SIGNAL(error(QString, QString)), this, SLOT(regError(QString, QString)));
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
    connect(image, SIGNAL(error(QString, QString)), this, SLOT(regError(QString, QString)));
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
    QStringList filenames = fileOpener.getOpenFileNames(this, "Add File(s)", QString(), "Nifti (*.nii *.nii.gz)");

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
    advancedOptionsWindow->reset(this->getCurrentAlgo());
    advancedOptionsWindow->show();
}


// Algo combo box changed
void milxQtRegistrationWindow::algoComboChange(int newIndex)
{
    if (newIndex == AffineItk)
    {
        setAlgo(AffineItk);
    }
    else if (newIndex == DemonItk)
    {
        setAlgo(DemonItk);
    }
    else if (newIndex == AladinNifti)
    {
        setAlgo(AladinNifti);
    }
    else if (newIndex == F3DNifti)
    {
        setAlgo(F3DNifti);
    }
    else if (newIndex == ElastixAffine)
    {
        setAlgo(ElastixAffine);
    }
    else if (newIndex == ElastixBSpline)
    {
        setAlgo(ElastixBSpline);
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

    // Check if the directory is valid and exist, otherwise we ask to create it
    QString outputDir = this->ui.inputDirectoryBrowser->text();
    if (!QFile::exists(outputDir))
    {
        // If the directory doesn't exist we ask if we should create it
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "Directory Creation", "The directory \"" + outputDir + "\" doesn't exist, do you want to create it ?", QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes)
        {
            QDir().mkdir(outputDir);
        }
    }

    // If the directory is invalid and still doesn't exist, error message
    if (!QDir(outputDir).exists())
    {
        QMessageBox msgBox(QMessageBox::Critical, "Registration", "Error the output directory path doesn't exist or is invalid.", QMessageBox::NoButton);
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
    // Perform the next registration
    performRegistrations();
}

// An error happenned during the registration
void milxQtRegistrationWindow::regError(QString functionName, QString errorMsg)
{
    // Show an error message
    QMessageBox msgBox(QMessageBox::Critical, "Registration error in function " + functionName, errorMsg, QMessageBox::NoButton);
    msgBox.exec();

    // If it's the average function, all the work is done
    if (functionName == "average()" || functionName == "computeAtlas()") {
        // Remove the atlas file it failed
        if (QFile::exists(atlasPath))
        {
            QFile::remove(atlasPath);
        }

        // Work is done
        workCompleted();
    }

    // Perform the next registration
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

#ifdef USE_NIFTI_REG
        // If we have to write the similarities
        if (this->ui.checkBoxSimilarities->isChecked())
        {

            writeSimilarities();
        }
#endif


        // If we have to compute the average
        if (computeAverage)
        {
#ifdef USE_NIFTI_REG
            computeAtlas();
#endif
        }
        else
        {
            // The work is completed
            workCompleted();
        }
    }

}


#ifdef USE_NIFTI_REG

// Write the similarities file
void milxQtRegistrationWindow::writeSimilarities()
{
    // Find the reference images
    int indexRef = 0;
    for (indexRef = 0; indexRef < images.size(); indexRef++)
    {
        if (images[indexRef]->isRef())
        {
            break;
        }
    }

    // Write the similarity file before the registration (i == 0) and after the registration (i == 1)
    for (int i = 0; i < 2; i++)
    {
        // Create the file

        QString similarityFile;

        if (i == 0)
            similarityFile = images[indexRef]->createSimilarityFileBefore();
        else if (i == 1)
            similarityFile = images[indexRef]->createSimilarityFileAfter();

        // Write the CSV
        QFile file(similarityFile);

        if (file.open(QIODevice::ReadWrite))
        {
            QTextStream stream(&file);
            // Create the header
            stream << "Reference;Registered Image;Output Image;Computation Time;NMI;SSD;NCC;LNCC\n";

            // Loop through the images
            for (int j = 0; j < images.size(); j++)
            {
                if (images[j]->isWorkDone() || images[j]->isChecked())
                {
                    milxQtSimilarities similarities;

                    if (i == 0)
                        similarities = images[j]->similarities_before;
                    else if (i == 1)
                        similarities = images[j]->similarities_after;

                    // Write results in textfile
                    stream << images[indexRef]->getPath() << ";";
                    stream << images[j]->getPath() << ";";
                    stream << images[j]->getOutputPath() << ";";

                    if (i == 0)
                        stream << "0;";
                    else if (i == 1)
                        stream << images[j]->getDuration() << ";";


                    stream << similarities.nmi << ";";
                    stream << similarities.ssd << ";";
                    stream << similarities.ncc << ";";
                    stream << similarities.lncc;
                    stream << "\n";
                }
            }

            file.close();
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
    if (atlasPath == "")
    {
        regError("computeAtlas()", "Unable to create the atlas file");
        return;
    }

    // Start the computation
    regAlgos->average_async(atlasPath, filenames);
}
#endif

// everything has been completed
void milxQtRegistrationWindow::workCompleted()
{
    workInProgress = false;
    this->enableUI();

    // Message box: registration completed
    QMessageBox msgBox(QMessageBox::Information, "Registration", "Registration completed !", QMessageBox::NoButton);
    msgBox.exec();

    // Reset images
    for (int i = 0; i < images.size(); i++)
    {
        images[i]->reset();
    }

    // If the files were not opened in SMILI we open the output folder
    if (!this->ui.checkBoxOpenResults->isChecked())
    {
        QString folder = QDir(this->ui.inputDirectoryBrowser->text()).absolutePath();
        QDesktopServices::openUrl(QUrl::fromLocalFile(folder));
    }
}
