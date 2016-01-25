#include "milxQtRegistration.h"
#include "milxQtFile.h"


int milxQtRegistration::startRegistration()
{
    // Create the temporary files for the registration
    if (createFiles() == EXIT_FAILURE)
    {
        return EXIT_FAILURE;
    }

    // Start the registration
    workDone = false;

    // Keep the time of start
    this->startTime = QDateTime::currentMSecsSinceEpoch();

    if (this->type == AffineItk)
    {
        MainWindow->printInfo("New Itk Affine registration started.\n");
        this->regAlgos->affine_async(params);
    }
    else if (this->type == DemonItk)
    {
        MainWindow->printInfo("New Itk Demon registration started.\n");
        this->regAlgos->demon_async(params);
    }
#ifdef USE_NIFTI_REG
    else if (this->type == AladinNifti)
    {
        MainWindow->printInfo("New Nifti Aladin registration started.\n");
        this->regAlgos->aladin_async(params);
    }
    else if (this->type == F3DNifti)
    {
        MainWindow->printInfo("New Nifti F3D registration started.\n");
        this->regAlgos->f3d_async(params);
    }
#endif
#ifdef USE_ELASTIX
    else if (this->type == ElastixAffine)
    {
        MainWindow->printInfo("New Elastix Affine registration started.\n");
        this->regAlgos->elastix_async(params);
    }
    else if (this->type == ElastixBSpline)
    {
        MainWindow->printInfo("New Elastix Bspline registration started.\n");
        this->regAlgos->elastix_async(params);
    }
#endif
    return EXIT_SUCCESS;
}

// Create the files required for the registration
int milxQtRegistration::createFiles()
{
    // Image to register filename
    QString filename = QFileInfo(path).baseName();

    // Create the tmp files
    params.floatingName = createFile(QDir::tempPath() + QDir::separator() + "smili_reg_img_XXXXXX.nii.gz");
    params.referenceName = createFile(QDir::tempPath() + QDir::separator() + "smili_reg_ref_XXXXXX.nii.gz");

    // Output path
    params.outputName = createFile(this->outputFolder + QDir::separator() + filename + "_" + this->getAlgoName() + "_registered_XXXXXX.nii.gz");

    // If the file were not created properly
    if (params.referenceName == "" || params.floatingName == "" || params.outputName == "")
    {
        return EXIT_FAILURE;
    }

    // Create special files for nifti
    if (this->type == F3DNifti)
    {
        // tmp file for nifti
        params.cppOutputName = createFile(QDir::tempPath() + QDir::separator() + "smili_reg_cpp_XXXXXX.nii.gz");
        if (params.cppOutputName == "")
        {
            return EXIT_FAILURE;
        }

        // cpp output file for nifti
        if (params.cpp2Def)
        {
            params.defOutputName = createFile(this->outputFolder + QDir::separator() + filename + "_" + this->getAlgoName() + "_deffield_XXXXXX.nii.gz");
            if (params.defOutputName == "")
            {
                return EXIT_FAILURE;
            }
        }
    }

    // Create the special file for Elastix
    if ((this->type == ElastixAffine || this->type == ElastixBSpline) && params.customParameterFile == false)
    {
        params.parameterFile = createFile(QDir::tempPath() + QDir::separator() + "elastix_reg_params_XXXXXX.txt");
        if (params.parameterFile == "")
        {
            return EXIT_FAILURE;
        }
        params.outputFolder = this->outputFolder;
    }

    if ((this->type == ElastixAffine || this->type == ElastixBSpline) && params.customParameterFile == false)
    {
        // write the parameter file
        QFile outputFile(params.parameterFile);
        outputFile.open(QIODevice::WriteOnly);

        // Check it opened OK
        if (!outputFile.isOpen()) {
            // If not open !! error !
            algoError("createFiles()", "Unable to write elastix parameter file: \"" + params.parameterFile + "\"");
            return EXIT_FAILURE;
        }
        QTextStream outStream(&outputFile);
        outStream << params.parametersTxt;
        outputFile.close();
    }

    // Copy the reference image and the image to register
    if(this->window->isLoaded())
        copyAndReplace(this->window, params.floatingName);
    else
        copyAndReplace(this->path, params.floatingName);
    if(this->reference->getImage()->isLoaded())
        copyAndReplace(this->reference->getImage(), params.referenceName);
    else
        copyAndReplace(this->reference->getPath(), params.referenceName);

    return EXIT_SUCCESS;
}

// Return a string with the name of the algo
QString milxQtRegistration::getAlgoName()
{
    QString strTypeAlgo = "unknownalgo";
    if (this->type == AffineItk) strTypeAlgo = "affineitk";
    if (this->type == DemonItk) strTypeAlgo = "demonitk";
    if (this->type == F3DNifti) strTypeAlgo = "f3dnifti";
    if (this->type == AladinNifti) strTypeAlgo = "aladinnifti";
    if (this->type == ElastixAffine) strTypeAlgo = "elastixaffine";
    if (this->type == ElastixBSpline) strTypeAlgo = "elastixbspline";

    return strTypeAlgo;
}

// Create the file for the Atlas and return the path
QString milxQtRegistration::createAtlasFile()
{
    QString filename = QFileInfo(path).baseName();
    return createFile(this->outputFolder + QDir::separator() + filename + "_" + this->getAlgoName() + "_atlas_XXXXXX.nii");
}

// Create the similarity file containing the similarity after the registration was performed and return the filepath
QString milxQtRegistration::createSimilarityFileAfter()
{
    QString filename = QFileInfo(path).baseName();
    return createFile(this->outputFolder + QDir::separator() + filename + "_" + this->getAlgoName() + "_similarities_XXXXXX.csv");
}

// Create the similarity file containing the similarity before the registration was performed and return the filepath
QString milxQtRegistration::createSimilarityFileBefore()
{
    QString filename = QFileInfo(path).baseName();
    return createFile(this->outputFolder + QDir::separator() + filename + "_" + this->getAlgoName() + "_similarities_before_XXXXXX.csv");
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


void milxQtRegistration::setComputeSimilarities(bool value)
{
    this->computeSimilarities = value;
}

// Return the duration of the registration
qint64 milxQtRegistration::getDuration()
{
    return (this->stopTime - this->startTime);
}


// Copy and replace
void milxQtRegistration::copyAndReplace(QString src, QString dst)
{
    milxQtFile inFile, outFile;
    milxQtImage *srcImage;

    inFile.openImage(src, srcImage);

    outFile.saveImage(dst, srcImage);
}

void milxQtRegistration::copyAndReplace(milxQtImage *src, QString dst)
{
    milxQtFile outFile;

    outFile.saveImage(dst, src);
}

// Create a temporary file and store the path in a char array
QString milxQtRegistration::createFile(QString pathTemplate)
{
    QTemporaryFile tmpFile(pathTemplate);
    tmpFile.setAutoRemove(false);
    if (!tmpFile.open()) {
        algoError("createFile()", "Unable to create temporary file with template: \"" + pathTemplate + "\".\n");
        return "";
    }
    tmpFile.close();

    MainWindow->printDebug("File \"" + tmpFile.fileName() + "\" created");

    return tmpFile.fileName();
}

// Delete the temporary files
void milxQtRegistration::deleteTmpFiles()
{
    if (params.floatingName != "") {
        if (QFile::exists(params.floatingName))
        {
            QFile::remove(params.floatingName);
            MainWindow->printInfo("File \"" + params.floatingName + "\" removed");
        }
    }
    if (params.referenceName != "") {
        if (QFile::exists(params.referenceName))
        {
            QFile::remove(params.referenceName);
            MainWindow->printInfo("File \"" + params.referenceName + "\" removed");
        }
    }
    if (params.cppOutputName != "") {
        if (QFile::exists(params.cppOutputName))
        {
            QFile::remove(params.cppOutputName);
            MainWindow->printInfo("File \"" + params.cppOutputName + "\" removed");
        }
    }

    if (params.customParameterFile == false && params.parameterFile != "")
    {
        if (QFile::exists(params.parameterFile))
        {
            QFile::remove(params.parameterFile);
            MainWindow->printInfo("File \"" + params.parameterFile + "\" removed");
        }
    }

    params.floatingName = "";
    params.referenceName = "";
    params.cppOutputName = "";
    params.parameterFile = "";
}


void milxQtRegistration::registrationCompleted()
{
    MainWindow->printInfo("---------- Registration Completed ----------\n");

    // Keep the time of completition
    this->stopTime = QDateTime::currentMSecsSinceEpoch();

    // We open the result
    if (openResults)
    {
        MainWindow->loadFile(params.outputName);
    }

#ifdef USE_NIFTI_REG
    // If we have to calculate the cpp2def, we start it
    if (this->type == F3DNifti && params.cpp2Def)
    {
        MainWindow->printInfo("Deformation field computation");

        // Start cpp2def
        this->regAlgos->cpp2def_async(params);

        return;
    }

    if (computeSimilarities)
    {
        this->regAlgos->similarities_async(this);
        return;
    }
#endif

    // Once done we mark the image as done and we uncheck it
    this->workDone = true;
    this->checked = false;

    // Otherwise we delete the temporary files
    deleteTmpFiles();

    // We emit the done signaml
    emit done();
}

void milxQtRegistration::algoError(QString functionName, QString errorMsg)
{
    // This image failed, but we marked it as done
    this->workDone = true;
    this->checked = false;

    // We delete the temporary files
    deleteTmpFiles();

    // Emit the error
    emit error(functionName, errorMsg);
}

#ifdef USE_NIFTI_REG

void milxQtRegistration::cpp2defCompleted()
{
    MainWindow->printInfo("Deformation field computed\n");

    if (computeSimilarities)
    {
        this->regAlgos->similarities_async(this);
        return;
    }

    // Once done we mark it as done and we unchecked the image
    this->workDone = true;
    this->checked = false;

    deleteTmpFiles();
    emit done();
}

void milxQtRegistration::similaritiesComputed()
{
    MainWindow->printInfo("Similarities computed\n");

    // Once done we mark it as done and we unchecked the image
    this->workDone = true;
    this->checked = false;

    deleteTmpFiles();
    emit done();
}
#endif

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
    // Delete regAlgos object
    delete regAlgos;
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
    regAlgos = new milxQtRegistrationAlgos(this);
    startTime = 0;
    stopTime = 0;
    similarities_after.lncc = 0;
    similarities_after.ncc = 0;
    similarities_after.nmi = 0;
    similarities_after.ssd = 0;
    similarities_before.lncc = 0;
    similarities_before.ncc = 0;
    similarities_before.nmi = 0;
    similarities_before.ssd = 0;
    computeSimilarities = false;

    // Add a listener (registration completed)
    connect(this->regAlgos, SIGNAL(registrationCompleted()), this, SLOT(registrationCompleted()));
    connect(this->regAlgos, SIGNAL(error(QString, QString)), this, SLOT(algoError(QString, QString)));

#ifdef USE_NIFTI_REG
    // Add a listener (cpp2def completed) if we have nifti library
    connect(this->regAlgos, SIGNAL(cpp2defCompleted()), this, SLOT(cpp2defCompleted()));

    // Add a listener for the similarities
    connect(this->regAlgos, SIGNAL(similaritiesComputed()), this, SLOT(similaritiesComputed()));
#endif
}

void milxQtRegistration::reset()
{
    params.reset();
    checked = false;
    workDone = false;
    startTime = 0;
    stopTime = 0;
    computeSimilarities = false;
    similarities_after.lncc = 0;
    similarities_after.ncc = 0;
    similarities_after.nmi = 0;
    similarities_after.ssd = 0;
    similarities_before.lncc = 0;
    similarities_before.ncc = 0;
    similarities_before.nmi = 0;
    similarities_before.ssd = 0;
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

// Get the path of the image to register
QString milxQtRegistration::getPath()
{
    return path;
}

// get the ouput path of the image registeres
QString milxQtRegistration::getOutputPath()
{
    return params.outputName;
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

void milxQtRegistration::setParams(milxQtRegistrationParams regparams)
{
    params = regparams;
}

void milxQtRegistration::setRegType(RegType reg)
{
    this->type = reg;
}

void milxQtRegistration::setReference(milxQtRegistration * ref)
{
    reference = ref;
}

bool milxQtRegistration::isOpened()
{
    return openedImage;
}
