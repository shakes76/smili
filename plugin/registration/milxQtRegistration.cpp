#include "milxQtRegistration.h"


void milxQtRegistration::startRegistration()
{
    // Starts an image registration
    MainWindow->printInfo("---------- New Image Registration ----------\n");

    // Create the temporary files for the registration
    createFiles();

    // Add a listener (registration completed)
    connect(this->niftiReg, SIGNAL(registrationCompleted()), this, SLOT(registrationCompleted()));

    // Start the registration
    workDone = false;

    if (type == Aladin)
    {
        MainWindow->printInfo("Start Aladin Registration\n");
        this->niftiReg->aladin_async(paramsAladin);
    }
    else if (type == F3D)
    {
        MainWindow->printInfo("Start F3D Registration\n");
        this->niftiReg->f3d_async(paramsF3D);
    }
}

// Create the files required for the registration
void milxQtRegistration::createFiles()
{
    // Image to register filename
    QString filename = QFileInfo(path).baseName();

    // Create the tmp files
    tmp_img = createFile(QDir::tempPath() + "/smili_reg_img_XXXXXX.nii");
    tmp_ref = createFile(QDir::tempPath() + "/smili_reg_ref_XXXXXX.nii");
    tmp_cpp = createFile(QDir::tempPath() + "/smili_reg_cpp_XXXXXX.nii");

    // Create the output files
    if (type == F3D && paramsF3D.cpp2Def)
    {
        output_def = createFile(this->outputFolder + "/" + filename + "_deffield_XXXXXX.nii");
    }
    output_path = createFile(this->outputFolder + "/" + filename + "_registered_XXXXXX.nii");

    // Copy the filepath in the registration parameters
    if (type == F3D)
    {
        copyPath(paramsF3D.floatingName, tmp_img);
        copyPath(paramsF3D.referenceName, tmp_ref);
        copyPath(paramsF3D.outputControlPointGridName, tmp_cpp);
        copyPath(paramsF3D.outputWarpedName, output_path);
        if (paramsF3D.cpp2Def)
        {
            copyPath(paramsF3D.defOutputName, output_def);
        }
    }
    else if (type == Aladin)
    {
        copyPath(paramsAladin.floatingName, tmp_img);
        copyPath(paramsAladin.referenceName, tmp_ref);
        copyPath(paramsAladin.outputResultName, output_path);
    }

    // Copy the reference image and the image to register
    copyAndReplace(this->path, tmp_img);
    copyAndReplace(this->reference->getPath(), tmp_ref);
}

// Create the file for the Atlas and return the path
QString milxQtRegistration::createAtlasFile()
{
    QString filename = QFileInfo(path).baseName();
    return createFile(this->outputFolder + "/" + filename + "_atlas_XXXXXX.nii");
}

// Copy a path from a QString to a char *
void milxQtRegistration::copyPath(char dest[FILENAME_MAX], QString src)
{
    QByteArray byteArray = src.toAscii();
    strncpy(dest, byteArray.data(), FILENAME_MAX);
}

// Set the image as the reference
void milxQtRegistration::setIsRef(bool istheref)
{
    if (istheref == true) {
        this->checked = false;
    }
    this->isRefImg = istheref;
}

// Is the image the reference image
bool milxQtRegistration::isRef()
{
    return this->isRefImg;
}

// Copy and replace
void milxQtRegistration::copyAndReplace(QString src, QString dst)
{
    if (QFile::exists(dst))
    {
        QFile::remove(dst);
    }

    QFile::copy(src, dst);
}

// Create a temporary file and store the path in a char array
QString milxQtRegistration::createFile(QString pathTemplate)
{
    QTemporaryFile tmpFile(pathTemplate);
    tmpFile.setAutoRemove(false);
    if (!tmpFile.open()) {
        MainWindow->printError("Unable to create temporary file: "+ pathTemplate+"\n");
    }
    tmpFile.close();

    MainWindow->printInfo("File \"" + tmpFile.fileName() + "\" created");

    return tmpFile.fileName();
}

// Delete the temporary files
void milxQtRegistration::deleteTmpFiles()
{
    if (tmp_img != "") {
        QFile::remove(tmp_img);
    }
    if (tmp_ref != "") {
        QFile::remove(tmp_ref);
    }
    if (tmp_cpp != "") {
        QFile::remove(tmp_cpp);
    }

    MainWindow->printInfo("File \"" + tmp_img + "\" removed");
    MainWindow->printInfo("File \"" + tmp_ref + "\" removed");
    MainWindow->printInfo("File \"" + tmp_cpp + "\" removed");

    tmp_img = "";
    tmp_ref = "";
    tmp_cpp = "";
}


void milxQtRegistration::registrationCompleted()
{
    MainWindow->printInfo("---------- Registration Completed ----------\n");

    // We open the result
    if (openResults)
    {
        MainWindow->loadFile(this->output_path);
    }

    // If we have to calculate the cpp2def, we start it
    if (type == F3D && paramsF3D.cpp2Def)
    {
        MainWindow->printInfo("---------- Deformation Field Computation ----------\n");

        // Add a listener (cpp2def completed)
        connect(this->niftiReg, SIGNAL(cpp2defCompleted()), this, SLOT(cpp2defCompleted()));

        // Start cpp2def
        this->niftiReg->cpp2def_async(paramsF3D);

        return;
    }

    // Once done we mark the image as done and we uncheck it
    this->workDone = true;
    this->checked = false;

    // Otherwise we delete the temporary files
    deleteTmpFiles();

    // We emit the done signaml
    emit done();
}

void milxQtRegistration::cpp2defCompleted()
{
    MainWindow->printInfo("Deformation field computed\n");

    // We open the result
    if (openResults)
    {
        MainWindow->loadFile(this->output_def);
        QWidgetList windowsList = MainWindow->getListOfWindows();
        milxQtImage * window = qobject_cast<milxQtImage *>(windowsList[windowsList.size() - 1]);
        window->vectorField();
        window->close();
    }

    // Once done we mark it as done and we unchecked the image
    this->workDone = true;
    this->checked = false;

    deleteTmpFiles();
    emit done();
}

milxQtRegistration::milxQtRegistration(QObject * parent, milxQtImage * imageWindow, milxQtMain * mainW) : QObject(parent)
{
    init(imageWindow->getName());
    MainWindow = mainW;
    window = imageWindow;
    openedImage = true;
}


milxQtRegistration::milxQtRegistration(QObject * parent, QString filepath, milxQtMain * mainW) : QObject(parent)
{
    init(filepath);
    MainWindow = mainW;
}

milxQtRegistration::~milxQtRegistration()
{
    // Delete niftiReg object
    delete niftiReg;
}

void milxQtRegistration::setOpenResults(bool open)
{
    openResults = open;
}

void milxQtRegistration::init(QString filepath)
{
    path = filepath;
    openedImage = false;
    window = NULL;
    checked = true;
    workDone = false;
    openResults = false;
    isRefImg = false;
    type = None;
    outputFolder = QFileInfo(filepath).absolutePath();
    niftiReg = new milxQtRegistrationNifti(this);
    tmp_img = "";
    tmp_ref = "";
    tmp_cpp = "";
    output_def = "";
    output_path = "";
}

void milxQtRegistration::setOutputFolder(QString folder)
{
    if (folder == "") {
        outputFolder = QFileInfo(path).absolutePath();
    }
    else {
        outputFolder = folder;
    }
}

QString milxQtRegistration::getOutputFolder()
{
    return outputFolder;
}


QString milxQtRegistration::getPath()
{
    return path;
}

QString milxQtRegistration::getOutputPath()
{
    return output_path;
}

void milxQtRegistration::setChecked(bool state)
{
    checked = state;
    if (isRef()) {
        checked = false;
    }
}

bool milxQtRegistration::isChecked()
{
    return checked;
}

bool milxQtRegistration::isWorkDone()
{
    return workDone;
}

void milxQtRegistration::setParams(ParamsAladin paladin)
{
    paramsAladin = paladin;
}

void milxQtRegistration::setParams(ParamsF3D pf3d)
{
    paramsF3D = pf3d;
}

void milxQtRegistration::setRegType(RegType regtype)
{
    type = regtype;
}

void milxQtRegistration::setReference(milxQtRegistration * ref)
{
    reference = ref;
}

bool milxQtRegistration::isOpened()
{
    return openedImage;
}
