#include "milxQtRegistrationAlgos.h"
#include "milxQtRegistration.h"

milxQtRegistrationAlgos::milxQtRegistrationAlgos(QObject * parent) : QObject(parent)
{
    return;
}

void milxQtRegistrationAlgos::affine_async(milxQtRegistrationParams params)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::affine, params);
}

void milxQtRegistrationAlgos::demon_async(milxQtRegistrationParams params)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::demon, params);
}

int milxQtRegistrationAlgos::affine(milxQtRegistrationParams params)
{
    milx::Registration reg;

    char referenceName[FILENAME_MAX + 1],
         floatingName[FILENAME_MAX + 1],
         outputName[FILENAME_MAX + 1];

    strncpy(referenceName, params.referenceName.toLatin1().constData(), FILENAME_MAX);
    strncpy(floatingName, params.floatingName.toLatin1().constData(), FILENAME_MAX);
    strncpy(outputName, params.outputName.toLatin1().constData(), FILENAME_MAX);


    bool result = reg.RegisterAffineITK(std::string(referenceName),
                                        std::string(floatingName),
                                        std::string(outputName),
                                        false, std::string(""),
                                        std::string(""));

    if (result == true)
        emit registrationCompleted();
    else
        emit error("affine()", "Itk affine registration error");

    return 0;
}

int milxQtRegistrationAlgos::demon(milxQtRegistrationParams params)
{
    milx::Registration reg;

    char referenceName[FILENAME_MAX + 1],
         floatingName[FILENAME_MAX + 1],
         outputName[FILENAME_MAX + 1];

    strncpy(referenceName, params.referenceName.toLatin1().constData(), FILENAME_MAX);
    strncpy(floatingName, params.floatingName.toLatin1().constData(), FILENAME_MAX);
    strncpy(outputName, params.outputName.toLatin1().constData(), FILENAME_MAX);

    bool result = reg.RegisterMilxNonRigid(std::string(referenceName),
                                           std::string(floatingName),
                                           std::string(outputName),
                                           false);


    if (result == true)
        emit registrationCompleted();
    else
        emit error("demon()", "Itk demon registration error");

    return 0;
}

#ifdef USE_NIFTI
void milxQtRegistrationAlgos::cpp2def_async(milxQtRegistrationParams params)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::cpp2def, params);
}

void milxQtRegistrationAlgos::f3d_async(milxQtRegistrationParams params)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::f3d, params);
}

void milxQtRegistrationAlgos::aladin_async(milxQtRegistrationParams params)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::aladin, params);
}

void milxQtRegistrationAlgos::average_async(QString outputName, QStringList filenames)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::average, outputName, filenames);
}


void milxQtRegistrationAlgos::similarities_async(milxQtRegistration * image)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::similarities, image);
}

int milxQtRegistrationAlgos::average(QString outputName, QStringList filenames)
{
    if (filenames.size() < 2) {
        emit error("average()", "Error while computing average, invalid number of file: minimum 2");
        return EXIT_FAILURE;
    }

    // Transformt the QtStrings in a char * array
    char * * arguments;
    int nbargs = filenames.count() + 2;

    // Allocate the array and copy the strings
    arguments = (char * *)malloc(nbargs * sizeof(char *));
    for (int i = 0; i < nbargs; i++)
    {
        arguments[i] = (char *)malloc(sizeof(char) * (FILENAME_MAX + 1));

        QString path;
        if (i == 0) {
            path = "";
        }
        else if (i == 1) {
            path = outputName;
        }
        else {
            path = filenames[i - 2];
        }

        QByteArray byteArray = path.toAscii();
        strncpy(arguments[i], byteArray.data(), FILENAME_MAX);
    }

    //Check the name of the first file to verify if they are analyse or nifti image
    std::string n(arguments[2]);
    if (n.find(".nii.gz") != std::string::npos ||
            n.find(".nii") != std::string::npos ||
            n.find(".hdr") != std::string::npos ||
            n.find(".img") != std::string::npos ||
            n.find(".img.gz") != std::string::npos)
    {
        // Input arguments are image filename
        // Read the first image to average
        nifti_image *tempImage = reg_io_ReadImageHeader(arguments[2]);
        if (tempImage == NULL) {
            // fprintf(stderr, "The following image can not be read: %s\n", arguments[2]);
            for (int i = 0; i < nbargs; i++) {
                free(arguments[i]);
            }
            free(arguments);
            emit error("average()", "Error while computing average, the following image can not be read: \"" + QString(arguments[2]) + "\".");

            return EXIT_FAILURE;
        }
        reg_checkAndCorrectDimension(tempImage);

        // Create the average image
        nifti_image *average_image = nifti_copy_nim_info(tempImage);
        nifti_image_free(tempImage);
        tempImage = NULL;
        average_image->datatype = NIFTI_TYPE_FLOAT32;
        if (sizeof(PrecisionTYPE) == sizeof(double))
            average_image->datatype = NIFTI_TYPE_FLOAT64;
        average_image->nbyper = sizeof(PrecisionTYPE);
        average_image->data = (void *)malloc(average_image->nvox*average_image->nbyper);
        reg_tools_addSubMulDivValue(average_image, average_image, 0.f, 2);

        int imageTotalNumber = 0;
        for (int i = 2; i<nbargs; ++i) {
            nifti_image *tempImage = reg_io_ReadImageFile(arguments[i]);
            if (tempImage == NULL) {
                //fprintf(stderr, "[!] The following image can not be read: %s\n", arguments[i]);
                for (int i = 0; i < nbargs; i++) {
                    free(arguments[i]);
                }
                free(arguments);
                emit error("average()", "Error while computing average, the following image can not be read: \"" + QString(arguments[i]) + "\".");
                return EXIT_FAILURE;
            }
            reg_checkAndCorrectDimension(tempImage);
            if (average_image->nvox != tempImage->nvox) {
                //fprintf(stderr, "[!] All images must have the same size. Error when processing: %s\n", arguments[i]);
                for (int i = 0; i < nbargs; i++) {
                    free(arguments[i]);
                }
                free(arguments);
                emit error("average()", "Error while computing average, all images must have the same size. Error when processing: \"" + QString(arguments[i]) + "\".");
                return EXIT_FAILURE;
            }
            reg_tools_addSubMulDivImages(average_image, tempImage, average_image, 0);
            imageTotalNumber++;
            nifti_image_free(tempImage);
            tempImage = NULL;
        }
        reg_tools_addSubMulDivValue(average_image, average_image, (float)imageTotalNumber, 3);
        reg_io_WriteImageFile(average_image, arguments[1]);
        nifti_image_free(average_image);
    }
    else {
        // input arguments are assumed to be text file name
        // Create an mat44 array to store all input matrices
        const size_t matrixNumber = nbargs - 2;
        mat44 *inputMatrices = (mat44 *)malloc(matrixNumber * sizeof(mat44));
        // Read all the input matrices
        for (size_t m = 0; m<matrixNumber; ++m) {
            if (FILE *aff = fopen(arguments[m + 2], "r")) {
                fclose(aff);
            }
            else {
                //fprintf(stderr, "The specified input affine file (%s) can not be read\n", arguments[m + 2]);
                for (int i = 0; i < nbargs; i++) {
                    free(arguments[i]);
                }
                free(arguments);
                emit error("average()", "Error while computing average, The specified input affine file (\"" + QString(arguments[m + 2]) + "\") can not be read.");
                return EXIT_FAILURE;
            }
            // Read the current matrix file
            std::ifstream affineFile;
            affineFile.open(arguments[m + 2]);
            if (affineFile.is_open()) {
                // Transfer the values into the mat44 array
                int i = 0;
                float value1, value2, value3, value4;
                while (!affineFile.eof()) {
                    affineFile >> value1 >> value2 >> value3 >> value4;
                    inputMatrices[m].m[i][0] = value1;
                    inputMatrices[m].m[i][1] = value2;
                    inputMatrices[m].m[i][2] = value3;
                    inputMatrices[m].m[i][3] = value4;
                    i++;
                    if (i>3) break;
                }
            }
            affineFile.close();
        }
        // All the input matrices are log-ed
        for (size_t m = 0; m<matrixNumber; ++m) {
            inputMatrices[m] = reg_mat44_logm(&inputMatrices[m]);
        }
        // All the exponentiated matrices are summed up into one matrix
        //temporary double are used to avoid error accumulation
        double tempValue[16] = { 0, 0, 0, 0,
                                 0, 0, 0, 0,
                                 0, 0, 0, 0,
                                 0, 0, 0, 0
                               };
        for (size_t m = 0; m<matrixNumber; ++m) {
            tempValue[0] += (double)inputMatrices[m].m[0][0];
            tempValue[1] += (double)inputMatrices[m].m[0][1];
            tempValue[2] += (double)inputMatrices[m].m[0][2];
            tempValue[3] += (double)inputMatrices[m].m[0][3];
            tempValue[4] += (double)inputMatrices[m].m[1][0];
            tempValue[5] += (double)inputMatrices[m].m[1][1];
            tempValue[6] += (double)inputMatrices[m].m[1][2];
            tempValue[7] += (double)inputMatrices[m].m[1][3];
            tempValue[8] += (double)inputMatrices[m].m[2][0];
            tempValue[9] += (double)inputMatrices[m].m[2][1];
            tempValue[10] += (double)inputMatrices[m].m[2][2];
            tempValue[11] += (double)inputMatrices[m].m[2][3];
            tempValue[12] += (double)inputMatrices[m].m[3][0];
            tempValue[13] += (double)inputMatrices[m].m[3][1];
            tempValue[14] += (double)inputMatrices[m].m[3][2];
            tempValue[15] += (double)inputMatrices[m].m[3][3];
        }
        // Average matrix is computed
        tempValue[0] /= (double)matrixNumber;
        tempValue[1] /= (double)matrixNumber;
        tempValue[2] /= (double)matrixNumber;
        tempValue[3] /= (double)matrixNumber;
        tempValue[4] /= (double)matrixNumber;
        tempValue[5] /= (double)matrixNumber;
        tempValue[6] /= (double)matrixNumber;
        tempValue[7] /= (double)matrixNumber;
        tempValue[8] /= (double)matrixNumber;
        tempValue[9] /= (double)matrixNumber;
        tempValue[10] /= (double)matrixNumber;
        tempValue[11] /= (double)matrixNumber;
        tempValue[12] /= (double)matrixNumber;
        tempValue[13] /= (double)matrixNumber;
        tempValue[14] /= (double)matrixNumber;
        tempValue[15] /= (double)matrixNumber;
        // The final matrix is exponentiated
        mat44 outputMatrix;
        outputMatrix.m[0][0] = (float)tempValue[0];
        outputMatrix.m[0][1] = (float)tempValue[1];
        outputMatrix.m[0][2] = (float)tempValue[2];
        outputMatrix.m[0][3] = (float)tempValue[3];
        outputMatrix.m[1][0] = (float)tempValue[4];
        outputMatrix.m[1][1] = (float)tempValue[5];
        outputMatrix.m[1][2] = (float)tempValue[6];
        outputMatrix.m[1][3] = (float)tempValue[7];
        outputMatrix.m[2][0] = (float)tempValue[8];
        outputMatrix.m[2][1] = (float)tempValue[9];
        outputMatrix.m[2][2] = (float)tempValue[10];
        outputMatrix.m[2][3] = (float)tempValue[11];
        outputMatrix.m[3][0] = (float)tempValue[12];
        outputMatrix.m[3][1] = (float)tempValue[13];
        outputMatrix.m[3][2] = (float)tempValue[14];
        outputMatrix.m[3][3] = (float)tempValue[15];
        outputMatrix = reg_mat44_expm(&outputMatrix);
        // Free the array containing the input matrices
        free(inputMatrices);
        // The final matrix is saved
        reg_tool_WriteAffineFile(&outputMatrix, arguments[1]);
    }

    // Free variables
    for (int i = 0; i < nbargs; i++) {
        free(arguments[i]);
    }
    free(arguments);
    emit averageCompleted();

    return EXIT_SUCCESS;
}

int milxQtRegistrationAlgos::cpp2def(milxQtRegistrationParams params)
{

    char referenceName[FILENAME_MAX + 1],
         defOutputName[FILENAME_MAX + 1],
         cppOutputName[FILENAME_MAX + 1];

    strncpy(referenceName, params.referenceName.toLatin1().constData(), FILENAME_MAX);
    strncpy(defOutputName, params.defOutputName.toLatin1().constData(), FILENAME_MAX);
    strncpy(cppOutputName, params.cppOutputName.toLatin1().constData(), FILENAME_MAX);


    /* Read the reference image */
    nifti_image *referenceImage = reg_io_ReadImageHeader(referenceName);
    if (referenceImage == NULL) {
        emit error("cpp2def()", "Error while while computing deformation field, the reference image (\"" + QString(referenceName) + "\") can not be read.");
        return 1;
    }
    reg_checkAndCorrectDimension(referenceImage);


    // Read the control point image
    nifti_image *controlPointImage = reg_io_ReadImageFile(cppOutputName);
    if (controlPointImage == NULL) {
        emit error("cpp2def()", "Error while while computing deformation field, the control point image (\"" + QString(cppOutputName) + "\") can not be read.");
        return 1;
    }
    reg_checkAndCorrectDimension(controlPointImage);
    // Allocate the deformation field
    nifti_image *deformationFieldImage = nifti_copy_nim_info(referenceImage);
    deformationFieldImage->dim[0] = deformationFieldImage->ndim = 5;
    deformationFieldImage->dim[1] = deformationFieldImage->nx = referenceImage->nx;
    deformationFieldImage->dim[2] = deformationFieldImage->ny = referenceImage->ny;
    deformationFieldImage->dim[3] = deformationFieldImage->nz = referenceImage->nz;
    deformationFieldImage->dim[4] = deformationFieldImage->nt = 1;
    deformationFieldImage->pixdim[4] = deformationFieldImage->dt = 1.0;
    if (referenceImage->nz>1) deformationFieldImage->dim[5] = deformationFieldImage->nu = 3;
    else deformationFieldImage->dim[5] = deformationFieldImage->nu = 2;
    deformationFieldImage->pixdim[5] = deformationFieldImage->du = 1.0;
    deformationFieldImage->dim[6] = deformationFieldImage->nv = 1;
    deformationFieldImage->pixdim[6] = deformationFieldImage->dv = 1.0;
    deformationFieldImage->dim[7] = deformationFieldImage->nw = 1;
    deformationFieldImage->pixdim[7] = deformationFieldImage->dw = 1.0;
    deformationFieldImage->nvox = deformationFieldImage->nx*deformationFieldImage->ny*deformationFieldImage->nz*deformationFieldImage->nt*deformationFieldImage->nu;
    deformationFieldImage->datatype = controlPointImage->datatype;
    deformationFieldImage->nbyper = controlPointImage->nbyper;
    deformationFieldImage->data = (void *)calloc(deformationFieldImage->nvox, deformationFieldImage->nbyper);
    //Computation of the deformation field
    if (controlPointImage->intent_code == NIFTI_INTENT_VECTOR &&
            strcmp(controlPointImage->intent_name, "NREG_VEL_STEP") == 0)
        reg_bspline_getDeformationFieldFromVelocityGrid(controlPointImage,
                deformationFieldImage
                                                       );
    else
        reg_spline_getDeformationField(controlPointImage,
                                       referenceImage,
                                       deformationFieldImage,
                                       NULL,
                                       false, //composition
                                       true // bspline
                                      );
    nifti_image_free(controlPointImage);
    // Ouput the deformation field image
    reg_io_WriteImageFile(deformationFieldImage, defOutputName);
    nifti_image_free(deformationFieldImage);


    nifti_image_free(referenceImage);
    emit cpp2defCompleted();
    return 0;
}

int milxQtRegistrationAlgos::f3d(milxQtRegistrationParams params)
{
    time_t start;
    time(&start);

    char referenceName[FILENAME_MAX + 1],
         floatingName[FILENAME_MAX + 1],
         outputName[FILENAME_MAX + 1],
         cppOutputName[FILENAME_MAX + 1];

    strncpy(referenceName, params.referenceName.toLatin1().constData(), FILENAME_MAX);
    strncpy(floatingName, params.floatingName.toLatin1().constData(), FILENAME_MAX);
    strncpy(outputName, params.outputName.toLatin1().constData(), FILENAME_MAX);
    strncpy(cppOutputName, params.cppOutputName.toLatin1().constData(), FILENAME_MAX);


    int maxiterationNumber = params.maxit;


    PrecisionTYPE spacing[3];
    spacing[0] = params.spacing[0];
    spacing[1] = params.spacing[1];
    spacing[2] = params.spacing[2];

    unsigned int levelNumber = params.ln;
    unsigned int levelToPerform = params.lp;
    bool noPyramid = params.nopy;
    bool useSym = params.useSym;

    char *referenceMaskName = NULL;
    char *inputControlPointGridName = NULL;
    char *affineTransformationName = NULL;
    bool flirtAffine = false;
    PrecisionTYPE bendingEnergyWeight = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE linearEnergyWeight0 = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE linearEnergyWeight1 = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE L2NormWeight = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE jacobianLogWeight = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    bool jacobianLogApproximation = true;
    PrecisionTYPE referenceSmoothingSigma = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE floatingSmoothingSigma = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE referenceThresholdUp[10];
    PrecisionTYPE referenceThresholdLow[10];
    PrecisionTYPE floatingThresholdUp[10];
    PrecisionTYPE floatingThresholdLow[10];
    unsigned int referenceBinNumber[10];
    unsigned int floatingBinNumber[10];
    for (int i = 0; i<10; i++) {
        referenceThresholdUp[i] = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
        referenceThresholdLow[i] = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
        floatingThresholdUp[i] = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
        floatingThresholdLow[i] = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
        referenceBinNumber[i] = 0;
        floatingBinNumber[i] = 0;
    }
    bool parzenWindowApproximation = true;
    PrecisionTYPE warpedPaddingValue = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    PrecisionTYPE gradientSmoothingSigma = std::numeric_limits<PrecisionTYPE>::quiet_NaN();
    bool verbose=true;
    bool useConjugate=true;
    bool useSSD=false;
    bool useKLD = false;
    int interpolation=1;
    bool xOptimisation=true;
    bool yOptimisation=true;
    bool zOptimisation=true;
    bool gridRefinement = true;

    bool additiveNMI = false;
    char *floatingMaskName = NULL;
    float inverseConsistencyWeight = std::numeric_limits<PrecisionTYPE>::quiet_NaN();

#ifdef _BUILD_NR_DEV
    int stepNumber = -1;
    bool useVel = false;
#endif

#ifdef _USE_CUDA
    bool useGPU = false;
    bool checkMem = false;
    int cardNumber = -1;
#endif


    if (referenceName[0] == '\0' || floatingName[0] == '\0') {
        //fprintf(stderr, "Err:\tThe reference and the floating image have to be defined.\n");
        emit error("f3d()", "F3D registration error: The reference and the floating image have to be defined.");
        return 1;
    }

#ifndef NDEBUG
    /*
    printf("[NiftyReg DEBUG] *******************************************\n");
    printf("[NiftyReg DEBUG] *******************************************\n");
    printf("[NiftyReg DEBUG] NiftyReg has been compiled in DEBUG mode\n");
    printf("[NiftyReg DEBUG] Please rerun cmake to set the variable\n");
    printf("[NiftyReg DEBUG] CMAKE_BUILD_TYPE to \"Release\" if required\n");
    printf("[NiftyReg DEBUG] *******************************************\n");
    printf("[NiftyReg DEBUG] *******************************************\n");
    */
#endif

    // Output the command line
#ifdef NDEBUG
    if (verbose == true) {
#endif
        /*
        printf("\n[NiftyReg F3DNifti] Command line:\n\t");
        for (int i = 0; i<argc; i++)
        	printf(" %s", argv[i]);
        printf("\n\n");
        */
#ifdef NDEBUG
    }
#endif

    // Read the reference image
    if (referenceName[0] == '\0') {
        //fprintf(stderr, "Error. No reference image has been defined\n");
        emit error("f3d()", "F3D registration error: No reference image has been defined.");
        return 1;
    }
    nifti_image *referenceImage = reg_io_ReadImageFile(referenceName);
    if (referenceImage == NULL) {
        //fprintf(stderr, "Error when reading the reference image %s\n", referenceName);
        emit error("f3d()", "F3D registration error: Error when reading the reference image \"" + QString(referenceName) + "\".");
        return 1;
    }

    // Read the floating image
    if (floatingName[0] == '\0') {
        //fprintf(stderr, "Error. No floating image has been defined\n");
        emit error("f3d()", "F3D registration error: No floating image has been defined.");
        return 1;
    }
    nifti_image *floatingImage = reg_io_ReadImageFile(floatingName);
    if (floatingImage == NULL) {
        //fprintf(stderr, "Error when reading the floating image %s\n", floatingName);
        emit error("f3d()", "F3D registration error: Error when reading the floating image \""+QString(floatingName)+"\".");
        return 1;
    }

    // Read the mask images
    nifti_image *referenceMaskImage = NULL;
    if (referenceMaskName != NULL) {
        referenceMaskImage = reg_io_ReadImageFile(referenceMaskName);
        if (referenceMaskImage == NULL) {
            //fprintf(stderr, "Error when reading the reference mask image %s\n", referenceMaskName);
            emit error("f3d()", "F3D registration error: Error when reading the reference mask image \"" + QString(referenceMaskName) + "\".");
            return 1;
        }
    }
    nifti_image *floatingMaskImage = NULL;
    if (floatingMaskName != NULL) {
        floatingMaskImage = reg_io_ReadImageFile(floatingMaskName);
        if (floatingMaskImage == NULL) {
            //fprintf(stderr, "Error when reading the reference mask image %s\n", floatingMaskName);
            emit error("f3d()", "F3D registration error: Error when reading the reference mask image \"" + QString(floatingMaskName) + "\".");
            return 1;
        }
    }

    // Read the input control point grid image
    nifti_image *controlPointGridImage = NULL;
    if (inputControlPointGridName != NULL) {
        controlPointGridImage = reg_io_ReadImageFile(inputControlPointGridName);
        if (controlPointGridImage == NULL) {
            //fprintf(stderr, "Error when reading the input control point grid image %s\n", inputControlPointGridName);
            emit error("f3d()", "F3D registration error: Error when reading the input control point grid image \"" + QString(inputControlPointGridName) + "\".");
            return 1;
        }
#ifdef _BUILD_NR_DEV
        if (controlPointGridImage->intent_code == NIFTI_INTENT_VECTOR &&
                strcmp(controlPointGridImage->intent_name, "NREG_VEL_STEP") == 0 &&
                fabs(controlPointGridImage->intent_p1)>1)
            useVel = true;
#endif
    }

    // Read the affine transformation
    mat44 *affineTransformation = NULL;
    if (affineTransformationName != NULL) {
        affineTransformation = (mat44 *)malloc(sizeof(mat44));
        // Check first if the specified affine file exist
        if (FILE *aff = fopen(affineTransformationName, "r")) {
            fclose(aff);
        }
        else {
            //fprintf(stderr, "The specified input affine file (%s) can not be read\n", affineTransformationName);
            emit error("f3d()", "F3D registration error: The specified input affine file (\"" + QString(affineTransformationName) + "\") can not be read.");
            return 1;
        }
        reg_tool_ReadAffineFile(affineTransformation,
                                referenceImage,
                                floatingImage,
                                affineTransformationName,
                                flirtAffine);
    }

    // Create the reg_f3d object
    reg_f3d<PrecisionTYPE> *REG = NULL;
#ifdef _USE_CUDA
    CUdevice dev;
    CUcontext ctx;
    if (useGPU) {

        if (linearEnergyWeight0 == linearEnergyWeight0 ||
                linearEnergyWeight1 == linearEnergyWeight1 ||
                L2NormWeight == L2NormWeight) {
            //fprintf(stderr, "NiftyReg ERROR CUDA] The linear elasticity has not been implemented with CUDA yet. Exit.\n");
            emit error("f3d()", "F3D registration error: NiftyReg ERROR CUDA, The linear elasticity has not been implemented with CUDA yet.");
            return 1;
        }

        if (useSym) {
            //fprintf(stderr, "\n[NiftyReg ERROR CUDA] GPU implementation of the symmetric registration is not available yet. Exit\n");
            emit error("f3d()", "F3D registration error: NiftyReg ERROR CUDA, GPU implementation of the symmetric registration is not available yet.");
            return 1;
        }
#ifdef _BUILD_NR_DEV
        if (useVel) {
            //fprintf(stderr, "\n[NiftyReg ERROR CUDA] GPU implementation of velocity field parametrisartion is not available yet. Exit\n");
            emit error("f3d()", "F3D registration error: NiftyReg ERROR CUDA, GPU implementation of velocity field parametrisartion is not available yet.");
            return 1;
        }
#endif

        if ((referenceImage->dim[4] == 1 && floatingImage->dim[4] == 1) || (referenceImage->dim[4] == 2 && floatingImage->dim[4] == 2)) {

            // The CUDA card is setup
            cuInit(0);
            struct cudaDeviceProp deviceProp;
            int device_count = 0;
            cudaGetDeviceCount(&device_count);
            int device = cardNumber;
            if(cardNumber==-1) {
                // following code is from cutGetMaxGflopsDeviceId()
                int max_gflops_device = 0;
                int max_gflops = 0;
                int current_device = 0;
                while( current_device < device_count ) {
                    cudaGetDeviceProperties( &deviceProp, current_device );
                    int gflops = deviceProp.multiProcessorCount * deviceProp.clockRate;
                    if (gflops > max_gflops) {
                        max_gflops = gflops;
                        max_gflops_device = current_device;
                    }
                    ++current_device;
                }
                device = max_gflops_device;
            }
            NR_CUDA_SAFE_CALL(cudaSetDevice( device ));
            NR_CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, device ));
            cuDeviceGet(&dev,device);
            cuCtxCreate(&ctx, 0, dev);
            if (deviceProp.major < 1) {
                //printf("[NiftyReg ERROR CUDA] The specified graphical card does not exist.\n");
                emit error("f3d()", "F3D registration error: NiftyReg ERROR CUDA, The specified graphical card does not exist.");
                return 1;
            }
            REG = new reg_f3d_gpu<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
            if (verbose == true) {
#endif
                //printf("\n[NiftyReg F3DNifti] GPU implementation is used\n");
#ifdef NDEBUG
            }
#endif
        }
        else {
            //fprintf(stderr, "[NiftyReg ERROR] The GPU implementation only handle 1 to 1 or 2 to 2 image(s) registration\n");
            emit error("f3d()", "F3D registration error: NiftyReg ERROR CUDA, The GPU implementation only handle 1 to 1 or 2 to 2 image(s) registration.");
            return 1;
        }
    }
    else
#endif
    {
        if(useSym) {
            REG = new reg_f3d_sym<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
            if (verbose == true) {
#endif // NDEBUG
                //printf("\n[NiftyReg F3DNifti SYM] CPU implementation is used\n");
#ifdef NDEBUG
            }
#endif // NDEBUG
        }
#ifdef _BUILD_NR_DEV
        else if (useVel) {
            REG = new reg_f3d2<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
            if(verbose==true) {
#endif
                //printf("\n[NiftyReg F3DNifti2] CPU implementation is used\n");
#ifdef NDEBUG
            }
#endif
        }
#endif // _BUILD_NR_DEV
        else {
            REG = new reg_f3d<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
            if (verbose == true) {
#endif // NDEBUG
                //printf("\n[NiftyReg F3DNifti] CPU implementation is used\n");
#ifdef NDEBUG
            }
#endif // NDEBUG
        }
    }
#ifdef _OPENMP
    int maxThreadNumber = omp_get_max_threads();
#ifdef NDEBUG
    if(verbose==true) {
#endif // NDEBUG
        printf("[NiftyReg F3DNifti] OpenMP is used with %i thread(s)\n", maxThreadNumber);
#ifdef NDEBUG
    }
#endif // NDEBUG
#endif // _OPENMP

    // Set the reg_f3d parameters

    REG->SetReferenceImage(referenceImage);

    REG->SetFloatingImage(floatingImage);

    if(verbose==false) REG->DoNotPrintOutInformation();
    else REG->PrintOutInformation();

    if (referenceMaskImage != NULL)
        REG->SetReferenceMask(referenceMaskImage);

    if(controlPointGridImage!=NULL)
        REG->SetControlPointGridImage(controlPointGridImage);

    if (affineTransformation != NULL)
        REG->SetAffineTransformation(affineTransformation);

    if(bendingEnergyWeight==bendingEnergyWeight)
        REG->SetBendingEnergyWeight(bendingEnergyWeight);

    if (linearEnergyWeight0 == linearEnergyWeight0 || linearEnergyWeight1 == linearEnergyWeight1) {
        if (linearEnergyWeight0 != linearEnergyWeight0) linearEnergyWeight0 = 0.0;
        if (linearEnergyWeight1 != linearEnergyWeight1) linearEnergyWeight1 = 0.0;
        REG->SetLinearEnergyWeights(linearEnergyWeight0, linearEnergyWeight1);
    }
    if (L2NormWeight == L2NormWeight)
        REG->SetL2NormDisplacementWeight(L2NormWeight);

    if(jacobianLogWeight==jacobianLogWeight)
        REG->SetJacobianLogWeight(jacobianLogWeight);

    if (jacobianLogApproximation)
        REG->ApproximateJacobianLog();
    else REG->DoNotApproximateJacobianLog();

    if (parzenWindowApproximation)
        REG->ApproximateParzenWindow();
    else REG->DoNotApproximateParzenWindow();

    if (maxiterationNumber>-1)
        REG->SetMaximalIterationNumber(maxiterationNumber);

    if (referenceSmoothingSigma == referenceSmoothingSigma)
        REG->SetReferenceSmoothingSigma(referenceSmoothingSigma);

    if(floatingSmoothingSigma==floatingSmoothingSigma)
        REG->SetFloatingSmoothingSigma(floatingSmoothingSigma);

    for (unsigned int t = 0; t<(unsigned int)referenceImage->nt; t++)
        if (referenceThresholdUp[t] == referenceThresholdUp[t])
            REG->SetReferenceThresholdUp(t, referenceThresholdUp[t]);

    for(unsigned int t=0; t<(unsigned int)referenceImage->nt; t++)
        if (referenceThresholdLow[t] == referenceThresholdLow[t])
            REG->SetReferenceThresholdLow(t, referenceThresholdLow[t]);

    for(unsigned int t=0; t<(unsigned int)floatingImage->nt; t++)
        if(floatingThresholdUp[t]==floatingThresholdUp[t])
            REG->SetFloatingThresholdUp(t, floatingThresholdUp[t]);

    for(unsigned int t=0; t<(unsigned int)floatingImage->nt; t++)
        if (floatingThresholdLow[t] == floatingThresholdLow[t])
            REG->SetFloatingThresholdLow(t, floatingThresholdLow[t]);

    for (unsigned int t = 0; t<(unsigned int)referenceImage->nt; t++)
        if (referenceBinNumber[t]>0)
            REG->SetReferenceBinNumber(t, referenceBinNumber[t]);

    for (unsigned int t = 0; t<(unsigned int)floatingImage->nt; t++)
        if (floatingBinNumber[t]>0)
            REG->SetFloatingBinNumber(t, floatingBinNumber[t]);

    if (warpedPaddingValue == warpedPaddingValue)
        REG->SetWarpedPaddingValue(warpedPaddingValue);

    for (unsigned int s = 0; s<3; s++)
        if (spacing[s] == spacing[s])
            REG->SetSpacing(s, spacing[s]);

    if (levelNumber>0)
        REG->SetLevelNumber(levelNumber);

    if (levelToPerform>0)
        REG->SetLevelToPerform(levelToPerform);

    if (gradientSmoothingSigma == gradientSmoothingSigma)
        REG->SetGradientSmoothingSigma(gradientSmoothingSigma);

    if (useSSD)
        REG->UseSSD();
    else REG->DoNotUseSSD();

    if (useKLD)
        REG->UseKLDivergence();
    else REG->DoNotUseKLDivergence();

    if (useConjugate == true)
        REG->UseConjugateGradient();
    else REG->DoNotUseConjugateGradient();

    if (noPyramid == 1)
        REG->DoNotUsePyramidalApproach();

    switch (interpolation) {
    case 0:
        REG->UseNeareatNeighborInterpolation();
        break;
    case 1:
        REG->UseLinearInterpolation();
        break;
    case 3:
        REG->UseCubicSplineInterpolation();
        break;
    default:
        REG->UseLinearInterpolation();
        break;
    }

    if (xOptimisation == false)
        REG->NoOptimisationAlongX();

    if (yOptimisation == false)
        REG->NoOptimisationAlongY();

    if (zOptimisation == false)
        REG->NoOptimisationAlongZ();

    if (gridRefinement == false)
        REG->NoGridRefinement();

    if (additiveNMI) REG->SetAdditiveMC();

    // F3DNifti SYM arguments
    if (floatingMaskImage != NULL) {
        if (useSym)
            REG->SetFloatingMask(floatingMaskImage);
#ifdef _BUILD_NR_DEV
        else if (useVel)
            REG->SetFloatingMask(floatingMaskImage);
#endif
        else {
            //fprintf(stderr, "[NiftyReg WARNING] The floating mask image is only used for symmetric or velocity field parametrisation\n");
            //fprintf(stderr, "[NiftyReg WARNING] The floating mask image is ignored\n");
        }
    }

    if (inverseConsistencyWeight == inverseConsistencyWeight)
        REG->SetInverseConsistencyWeight(inverseConsistencyWeight);

#ifdef _BUILD_NR_DEV
    // F3DNifti2 arguments
    if (stepNumber>-1)
        REG->SetCompositionStepNumber(stepNumber);
#endif

    // Run the registration
#ifdef _USE_CUDA
    if (useGPU && checkMem) {
        size_t free, total, requiredMemory = REG->CheckMemoryMB_f3d();
        cuMemGetInfo(&free, &total);
        //printf("[NiftyReg CUDA] The required memory to run the registration is %lu Mb\n",
        (unsigned long int)requiredMemory);
        //printf("[NiftyReg CUDA] The GPU card has %lu Mb from which %lu Mb are currenlty free\n",
        (unsigned long int)total / (1024 * 1024), (unsigned long int)free / (1024 * 1024));
    }
    else {
#endif
        REG->Run_f3d();

        // Save the control point result
        nifti_image *outputControlPointGridImage = REG->GetControlPointPositionImage();
        if (cppOutputName[0] == '\0') strcpy(cppOutputName, "outputCPP.nii");
        memset(outputControlPointGridImage->descrip, 0, 80);
        strcpy(outputControlPointGridImage->descrip, "Control point position from NiftyReg (reg_f3d)");
#ifdef _BUILD_NR_DEV
        if (useVel)
            strcpy(outputControlPointGridImage->descrip, "Velocity field grid from NiftyReg (reg_f3d2)");
#endif
        reg_io_WriteImageFile(outputControlPointGridImage, cppOutputName);
        nifti_image_free(outputControlPointGridImage);
        outputControlPointGridImage = NULL;

        // Save the backward control point result
#ifdef _BUILD_NR_DEV
        if (useSym || useVel) {
#else
        if (useSym) {
#endif
            // _backward is added to the forward control point grid image name
            std::string b(cppOutputName);
            if (b.find(".nii.gz") != std::string::npos)
                b.replace(b.find(".nii.gz"), 7, "_backward.nii.gz");
            else if (b.find(".nii") != std::string::npos)
                b.replace(b.find(".nii"), 4, "_backward.nii");
            else if (b.find(".hdr") != std::string::npos)
                b.replace(b.find(".hdr"), 4, "_backward.hdr");
            else if (b.find(".img.gz") != std::string::npos)
                b.replace(b.find(".img.gz"), 7, "_backward.img.gz");
            else if (b.find(".img") != std::string::npos)
                b.replace(b.find(".img"), 4, "_backward.img");
#ifdef _USE_NR_PNG
            else if (b.find(".png") != std::string::npos)
                b.replace(b.find(".png"), 4, "_backward.png");
#endif
#ifdef _USE_NR_NRRD
            else if (b.find(".nrrd") != std::string::npos)
                b.replace(b.find(".nrrd"), 5, "_backward.nrrd");
#endif
            else b.append("_backward.nii");
            nifti_image *outputBackwardControlPointGridImage = REG->GetBackwardControlPointPositionImage();
            memset(outputBackwardControlPointGridImage->descrip, 0, 80);
            strcpy(outputBackwardControlPointGridImage->descrip, "Backward Control point position from NiftyReg (reg_f3d)");
#ifdef _BUILD_NR_DEV
            if (useVel)
                strcpy(outputBackwardControlPointGridImage->descrip, "Backward velocity field grid from NiftyReg (reg_f3d2)");
#endif
            reg_io_WriteImageFile(outputBackwardControlPointGridImage, b.c_str());
            nifti_image_free(outputBackwardControlPointGridImage);
            outputBackwardControlPointGridImage = NULL;
        }

        // Save the warped image result(s)
        nifti_image **outputWarpedImage = (nifti_image **)malloc(2 * sizeof(nifti_image *));
        outputWarpedImage[0] = outputWarpedImage[1] = NULL;
        outputWarpedImage = REG->GetWarpedImage();
        if (outputName[0] == '\0') strcpy(outputName, "outputResult.nii");
        memset(outputWarpedImage[0]->descrip, 0, 80);
        strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d)");
        if (useSym) {
            strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
            strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
#ifdef _BUILD_NR_DEV
        }
        if (useVel) {
            strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d2)");
            strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d2)");
        }
        if (useSym || useVel) {
#endif
            if (outputWarpedImage[1] != NULL) {
                std::string b(outputName);
                if (b.find(".nii.gz") != std::string::npos)
                    b.replace(b.find(".nii.gz"), 7, "_backward.nii.gz");
                else if (b.find(".nii") != std::string::npos)
                    b.replace(b.find(".nii"), 4, "_backward.nii");
                else if (b.find(".hdr") != std::string::npos)
                    b.replace(b.find(".hdr"), 4, "_backward.hdr");
                else if (b.find(".img.gz") != std::string::npos)
                    b.replace(b.find(".img.gz"), 7, "_backward.img.gz");
                else if (b.find(".img") != std::string::npos)
                    b.replace(b.find(".img"), 4, "_backward.img");
#ifdef _USE_NR_PNG
                else if (b.find(".png") != std::string::npos)
                    b.replace(b.find(".png"), 4, "_backward.png");
#endif
#ifdef _USE_NR_NRRD
                else if (b.find(".nrrd") != std::string::npos)
                    b.replace(b.find(".nrrd"), 5, "_backward.nrrd");
#endif
                else b.append("_backward.nii");
                if (useSym)
                    strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
#ifdef _BUILD_NR_DEV
                if (useVel)
                    strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d2)");
#endif
                reg_io_WriteImageFile(outputWarpedImage[1], b.c_str());
                nifti_image_free(outputWarpedImage[1]);
                outputWarpedImage[1] = NULL;
            }
        }
        reg_io_WriteImageFile(outputWarpedImage[0], outputName);
        nifti_image_free(outputWarpedImage[0]);
        outputWarpedImage[0] = NULL;
        free(outputWarpedImage);
#ifdef _USE_CUDA
    }
    cuCtxDetach(ctx);
#endif
    // Erase the registration object
    delete REG;

    // Clean the allocated images
    if (referenceImage != NULL) nifti_image_free(referenceImage);
    if (floatingImage != NULL) nifti_image_free(floatingImage);
    if (controlPointGridImage != NULL) nifti_image_free(controlPointGridImage);
    if (affineTransformation != NULL) free(affineTransformation);
    if (referenceMaskImage != NULL) nifti_image_free(referenceMaskImage);
    if (floatingMaskImage != NULL) nifti_image_free(floatingMaskImage);

#ifdef NDEBUG
    if (verbose) {
#endif
        time_t end;
        time(&end);
        int minutes = (int)floorf(float(end - start) / 60.0f);
        int seconds = (int)(end - start - 60 * minutes);

#ifdef _USE_CUDA
        if (!checkMem) {
#endif
            //printf("[NiftyReg F3DNifti] Registration Performed in %i min %i sec\n", minutes, seconds);
            //printf("[NiftyReg F3DNifti] Have a good day !\n");
#ifdef _USE_CUDA
        }
#endif
#ifdef NDEBUG
    }
#endif

    emit registrationCompleted();
    return 0;
}

int milxQtRegistrationAlgos::aladin(milxQtRegistrationParams params)
{
    time_t start;
    time(&start);

    char referenceName[FILENAME_MAX + 1],
         floatingName[FILENAME_MAX + 1],
         outputName[FILENAME_MAX + 1];

    strncpy(referenceName, params.referenceName.toLatin1().constData(), FILENAME_MAX);
    strncpy(floatingName, params.floatingName.toLatin1().constData(), FILENAME_MAX);
    strncpy(outputName, params.outputName.toLatin1().constData(), FILENAME_MAX);

    int referenceImageFlag = 1;
    int floatingImageFlag = 1;
    int outputResultFlag = 1;

    int maxIter = params.maxit;
    int symFlag = params.useSym;
    int nLevels =  params.ln;
    int levelsToPerform = params.lp;
    float blockPercentage = params.percentBlock;
    int affineFlag = 1;
    int rigidFlag = 1;

    if (params.rigOnly)
    {
        rigidFlag = 1;
        affineFlag = 0;
    }

    if (params.affineDirect)
    {
        rigidFlag = 0;
        affineFlag = 1;
    }

    char *outputAffineName = NULL;
    int outputAffineFlag = 0;

    char *inputAffineName = NULL;
    int inputAffineFlag = 0;
    int flirtAffineFlag = 0;

    char *referenceMaskName = NULL;
    int referenceMaskFlag = 0;

    char *floatingMaskName = NULL;
    int floatingMaskFlag = 0;

    float inlierLts=50.0f;
    int alignCentre=1;
    int interpolation = 1;
    float floatingSigma = 0.0;
    float referenceSigma=0.0;

    if (!referenceImageFlag || !floatingImageFlag) {
        //fprintf(stderr, "Err:\tThe reference and the floating image have to be defined.\n");
        emit error("aladin()", "Aladin registration error: The reference and the floating image have to be defined");
        return 1;
    }

    reg_aladin<PrecisionTYPE> *REG;
#ifdef _BUILD_NR_DEV
    if (symFlag)
    {
        REG = new reg_aladin_sym<PrecisionTYPE>;
        if ((referenceMaskFlag && !floatingMaskName) || (!referenceMaskFlag && floatingMaskName))
        {
            //fprintf(stderr, "[NiftyReg Warning] You have one image mask option turned on but not the other.\n");
            //fprintf(stderr, "[NiftyReg Warning] This will affect the degree of symmetry achieved.\n");
        }
    }
    else
    {
#endif
        REG = new reg_aladin<PrecisionTYPE>;
#ifdef _BUILD_NR_DEV
        if (floatingMaskFlag)
        {
            //fprintf(stderr, "Note: Floating mask flag only used in symmetric method. Ignoring this option\n");
        }
    }
#endif
    REG->SetMaxIterations(maxIter);
    REG->SetNumberOfLevels(nLevels);
    REG->SetLevelsToPerform(levelsToPerform);
    REG->SetReferenceSigma(referenceSigma);
    REG->SetFloatingSigma(floatingSigma);
    REG->SetAlignCentre(alignCentre);
    REG->SetPerformAffine(affineFlag);
    REG->SetPerformRigid(rigidFlag);
    REG->SetBlockPercentage(blockPercentage);
    REG->SetInlierLts(inlierLts);
    REG->SetInterpolation(interpolation);

    if (REG->GetLevelsToPerform() > REG->GetNumberOfLevels())
        REG->SetLevelsToPerform(REG->GetNumberOfLevels());

    /* Read the reference image and check its dimension */
    nifti_image *referenceHeader = reg_io_ReadImageFile(referenceName);
    if (referenceHeader == NULL) {
        //fprintf(stderr, "* ERROR Error when reading the reference  image: %s\n", referenceImageName);
        emit error("aladin()", "Aladin registration error: Error when reading the reference  image (\"" + QString(referenceName) + "\").");
        return 1;
    }

    /* Read teh floating image and check its dimension */
    nifti_image *floatingHeader = reg_io_ReadImageFile(floatingName);
    if (floatingHeader == NULL) {
        //fprintf(stderr, "* ERROR Error when reading the floating image: %s\n", floatingImageName);
        emit error("aladin()", "Aladin registration error: Error when reading the floating image (\"" + QString(floatingName) + "\").");
        return 1;
    }

    // Set the reference and floating image
    REG->SetInputReference(referenceHeader);
    REG->SetInputFloating(floatingHeader);

    // Set the input affine transformation if defined
    if (inputAffineFlag == 1)
        REG->SetInputTransform(inputAffineName, flirtAffineFlag);

    /* read the reference mask image */
    nifti_image *referenceMaskImage = NULL;
    if (referenceMaskFlag) {
        referenceMaskImage = reg_io_ReadImageFile(referenceMaskName);
        if (referenceMaskImage == NULL) {
            //fprintf(stderr, "* ERROR Error when reading the reference mask image: %s\n", referenceMaskName);
            emit error("aladin()", "Aladin registration error: Error when reading the reference mask image (\"" + QString(referenceMaskName) + "\").");
            return 1;
        }
        /* check the dimension */
        for (int i = 1; i <= referenceHeader->dim[0]; i++) {
            if (referenceHeader->dim[i] != referenceMaskImage->dim[i]) {
                //fprintf(stderr, "* ERROR The reference image and its mask do not have the same dimension\n");
                emit error("aladin()", "Aladin registration error: The reference image and its mask do not have the same dimension.");
                return 1;
            }
        }
        REG->SetInputMask(referenceMaskImage);
    }
#ifdef _BUILD_NR_DEV
    nifti_image *floatingMaskImage = NULL;
    if (floatingMaskFlag && symFlag) {
        floatingMaskImage = reg_io_ReadImageFile(floatingMaskName);
        if (floatingMaskImage == NULL) {
            //fprintf(stderr, "* ERROR Error when reading the floating mask image: %s\n", referenceMaskName);
            emit error("aladin()", "Aladin registration error: Error when reading the floating mask image (\"" + QString(floatingMaskName) + "\").");
            return 1;
        }
        /* check the dimension */
        for (int i = 1; i <= floatingHeader->dim[0]; i++) {
            if (floatingHeader->dim[i] != floatingMaskImage->dim[i]) {
                //fprintf(stderr, "* ERROR The floating image and its mask do not have the same dimension\n");
                emit error("aladin()", "Aladin registration error: Error the floating image and its mask do not have the same dimension.");
                return 1;
            }
        }
        REG->SetInputFloatingMask(floatingMaskImage);
    }
#endif
    REG->Run();

    // The warped image is saved
    nifti_image *outputResultImage = REG->GetFinalWarpedImage();
    if (!outputResultFlag) strcpy(outputName, "outputResult.nii");
    reg_io_WriteImageFile(outputResultImage, outputName);
    nifti_image_free(outputResultImage);

    /* The affine transformation is saved */
    if (outputAffineFlag)
        reg_tool_WriteAffineFile(REG->GetTransformationMatrix(), outputAffineName);
    else reg_tool_WriteAffineFile(REG->GetTransformationMatrix(), (char *)"outputAffine.txt");

    nifti_image_free(referenceHeader);
    nifti_image_free(floatingHeader);

    delete REG;
    time_t end;
    time(&end);
    int minutes = (int)floorf((end - start) / 60.0f);
    int seconds = (int)(end - start - 60 * minutes);


    emit registrationCompleted();
    return 0;
}

int milxQtRegistrationAlgos::similarities(milxQtRegistration * image)
{
    QString program = "reg_measure";
    milxQtSimilarities similarities;

    // We compute the similarity before the registration (i == 0) and after the registration (i == 1)
    for (int i = 0; i < 2; i++)
    {
        QStringList arguments;
        arguments << "-ref" << image->reference->getPath();

        if (i == 0)
            arguments << "-flo" << image->getPath();
        else if (i == 1)
            arguments << "-flo" << image->getOutputPath();

        QString similarityOutput = image->createFile(QDir::tempPath() + QDir::separator() + "similarity_XXXXXX.txt");
        arguments << "-out" << similarityOutput;

        arguments << "-ncc" << "-lncc" << "-nmi" << "-ssd";

        /*
        QString test2 = "";
        for (int j = 0; j < arguments.size(); j++)
        {
            QString test = arguments[j];
            test2 = test2 + test + " ";
        }
        */

        // Call reg_measure
        QProcess *myProcess = new QProcess(this);
        myProcess->start(program, arguments);
        myProcess->waitForFinished();

        // Read the results
        QFile inputFile(similarityOutput);
        if (inputFile.open(QIODevice::ReadOnly))
        {
            int lineNumber = 0;
            QTextStream in(&inputFile);
            while (!in.atEnd())
            {
                QString line = in.readLine();
                bool ok = false;
                double value = line.toDouble(&ok);

                if (ok && lineNumber == 0) {
                    similarities.ncc = value;
                }
                else if (ok && lineNumber == 1) {
                    similarities.lncc = value;
                }
                else if (ok && lineNumber == 2) {
                    similarities.nmi = value;
                }
                else if (ok && lineNumber == 3) {
                    similarities.ssd = value;
                }

                lineNumber++;
            }

            if (i == 0)
                image->similarities_before = similarities;
            else if (i == 1)
                image->similarities_after = similarities;


            inputFile.close();
        }

        // Remove the file
        if (QFile::exists(similarityOutput))
        {
            QFile::remove(similarityOutput);
        }

        delete myProcess;
    }


    emit similaritiesComputed();
    return 0;
}

#endif

#ifdef USE_ELASTIX

void milxQtRegistrationAlgos::elastix_async(milxQtRegistrationParams params)
{
    future = QtConcurrent::run(this, &milxQtRegistrationAlgos::elastix, params);
}

int milxQtRegistrationAlgos::elastix(milxQtRegistrationParams params)
{
    // Create the arguments from the parameters
    QStringList args;
    args.append("smilx.exe");

    args.append("-f");
    args.append(params.referenceName);

    args.append("-m");
    args.append(params.floatingName);

    args.append("-out");
    args.append(params.outputFolder);

    args.append("-p");
    args.append(params.parameterFile);

    // Some typedef's.
    typedef elx::ElastixMain                            ElastixMainType;
    typedef ElastixMainType::Pointer                    ElastixMainPointer;
    typedef std::vector< ElastixMainPointer >           ElastixMainVectorType;
    typedef ElastixMainType::ObjectPointer              ObjectPointer;
    typedef ElastixMainType::DataObjectContainerPointer DataObjectContainerPointer;
    typedef ElastixMainType::FlatDirectionCosinesType   FlatDirectionCosinesType;

    typedef ElastixMainType::ArgumentMapType ArgumentMapType;
    typedef ArgumentMapType::value_type      ArgumentMapEntryType;

    typedef std::pair< std::string, std::string > ArgPairType;
    typedef std::queue< ArgPairType >             ParameterFileListType;
    typedef ParameterFileListType::value_type     ParameterFileListEntryType;

    // Support Mevis Dicom Tiff (if selected in cmake)
    RegisterMevisDicomTiff();

    // Some declarations and initialisations.
    ElastixMainVectorType elastices;

    ObjectPointer              transform = 0;
    DataObjectContainerPointer fixedImageContainer = 0;
    DataObjectContainerPointer movingImageContainer = 0;
    DataObjectContainerPointer fixedMaskContainer = 0;
    DataObjectContainerPointer movingMaskContainer = 0;
    FlatDirectionCosinesType   fixedImageOriginalDirection;
    int                        returndummy = 0;
    unsigned long              nrOfParameterFiles = 0;
    ArgumentMapType            argMap;
    ParameterFileListType      parameterFileList;
    bool                       outFolderPresent = false;
    std::string                outFolder = "";
    std::string                logFileName = "";

    // Put command line parameters into parameterFileList.
    for (unsigned int i = 1; static_cast<long>(i) < (args.count() - 1); i += 2)
    {
        std::string key(args[i].toStdString());
        std::string value(args[i + 1].toStdString());

        if (key == "-p")
        {
            // Queue the ParameterFileNames.
            nrOfParameterFiles++;
            parameterFileList.push(
                ParameterFileListEntryType(key.c_str(), value.c_str()));
            // The different '-p' are stored in the argMap, with
            // keys p(1), p(2), etc.
            std::ostringstream tempPname("");
            tempPname << "-p(" << nrOfParameterFiles << ")";
            std::string tempPName = tempPname.str();
            argMap.insert(ArgumentMapEntryType(tempPName.c_str(), value.c_str()));
        }
        else
        {
            if (key == "-out")
            {
                // Make sure that last character of the output folder equals a '/' or '\'.
                const char last = value[value.size() - 1];
                if (last != '/' && last != '\\') {
                    value.append("/");
                }
                value = itksys::SystemTools::ConvertToOutputPath(value.c_str());

                // Note that on Windows, in case the output folder contains a space,
                // the path name is double quoted by ConvertToOutputPath, which is undesirable.
                // So, we remove these quotes again.

                if (itksys::SystemTools::StringStartsWith(value.c_str(), "\"")
                        && itksys::SystemTools::StringEndsWith(value.c_str(), "\""))
                {
                    value = value.substr(1, value.length() - 2);
                }

                // Save this information.
                outFolderPresent = true;
                outFolder = value;

            } // end if key == "-out"

            // Attempt to save the arguments in the ArgumentMap.
            if (argMap.count(key.c_str()) == 0)
            {
                argMap.insert(ArgumentMapEntryType(key.c_str(), value.c_str()));
            }
            else
            {
                // Duplicate arguments.
                std::cerr << "WARNING!" << std::endl;
                std::cerr << "Argument " << key.c_str() << "is only required once." << std::endl;
                std::cerr << "Arguments " << key.c_str() << " " << value.c_str() << "are ignored" << std::endl;
            }

        } // end else (so, if key does not equal "-p")

    } // end for loop

    // The argv0 argument, required for finding the component.dll/so's.
    argMap.insert(ArgumentMapEntryType("-argv0", args[0].toStdString()));

    // Check if at least once the option "-p" is given.
    if (nrOfParameterFiles == 0)
    {
        std::cerr << "ERROR: No CommandLine option \"-p\" given!" << std::endl;
        returndummy |= -1;
    }

    // Check if the -out option is given.
    if (outFolderPresent)
    {
        // Check if the output directory exists.
        bool outFolderExists = itksys::SystemTools::FileIsDirectory(outFolder.c_str());
        if (!outFolderExists)
        {
            std::cerr << "ERROR: the output directory \"" << outFolder << "\" does not exist." << std::endl;
            std::cerr << "You are responsible for creating it." << std::endl;
            returndummy |= -2;
        }
        else
        {
            // Setup xout.
            // CODE REMOVE FROM HERE, AND ADDED IN MILX REGISTRATION WINDOW, WE ONLY EXECUTE THE SETUP ONCE
            /*
            logFileName = outFolder + "elastix.log";
            int returndummy2 = elx::xoutSetup(logFileName.c_str(), false, false);
            if (returndummy2)
            {
                std::cerr << "ERROR while setting up xout." << std::endl;
            }
            returndummy |= returndummy2;
            */
        }
    }
    else
    {
        returndummy = -2;
        std::cerr << "ERROR: No CommandLine option \"-out\" given!" << std::endl;
    }

    // Stop if some fatal errors occurred.
    if (returndummy)
    {
        emit error("elastix()", "Elastix registration error, error code: " + QString::number(returndummy) + ".");
        return returndummy;
    }

    elxout << std::endl;

    // Declare a timer, start it and print the start time.
    tmr::Timer::Pointer totaltimer = tmr::Timer::New();
    totaltimer->StartTimer();
    elxout << "elastix is started at " << totaltimer->PrintStartTime()
           << ".\n" << std::endl;

    // Print where elastix was run.
    //	elxout << "which elastix:   " << argv[0] << std::endl;
    itksys::SystemInformation info;
    info.RunCPUCheck();
    info.RunOSCheck();
    info.RunMemoryCheck();
    elxout << "elastix runs at: " << info.GetHostname() << std::endl;
    elxout << "  " << info.GetOSName() << " "
           << info.GetOSRelease() << (info.Is64Bits() ? " (x64), " : ", ")
           << info.GetOSVersion() << std::endl;
    elxout << "  with " << info.GetTotalPhysicalMemory() << " MB memory, and "
           << info.GetNumberOfPhysicalCPU() << " cores @ "
           << static_cast< unsigned int >(info.GetProcessorClockFrequency())
           << " MHz." << std::endl;

    // ********************* START REGISTRATION *********************
    //
    // Do the (possibly multiple) registration(s).

    for (unsigned int i = 0; i < nrOfParameterFiles; i++)
    {
        // Create another instance of ElastixMain.
        elastices.push_back(ElastixMainType::New());

        // Set stuff we get from a former registration.
        elastices[i]->SetInitialTransform(transform);
        elastices[i]->SetFixedImageContainer(fixedImageContainer);
        elastices[i]->SetMovingImageContainer(movingImageContainer);
        elastices[i]->SetFixedMaskContainer(fixedMaskContainer);
        elastices[i]->SetMovingMaskContainer(movingMaskContainer);
        elastices[i]->SetOriginalFixedImageDirectionFlat(fixedImageOriginalDirection);

        // Set the current elastix-level.
        elastices[i]->SetElastixLevel(i);
        elastices[i]->SetTotalNumberOfElastixLevels(nrOfParameterFiles);

        // Delete the previous ParameterFileName.
        if (argMap.count("-p"))
        {
            argMap.erase("-p");
        }

        // Read the first parameterFileName in the queue.
        ArgPairType argPair = parameterFileList.front();
        parameterFileList.pop();

        // Put it in the ArgumentMap.
        argMap.insert(ArgumentMapEntryType(argPair.first, argPair.second));

        // Print a start message.
        elxout << "-------------------------------------------------------------------------" << "\n" << std::endl;
        elxout << "Running elastix with parameter file " << i
               << ": \"" << argMap["-p"] << "\".\n" << std::endl;

        // Declare a timer, start it and print the start time.
        tmr::Timer::Pointer timer = tmr::Timer::New();
        timer->StartTimer();
        elxout << "Current time: " << timer->PrintStartTime() << "." << std::endl;

        // Start registration.
        returndummy = elastices[i]->Run(argMap);

        // Check for errors.
        if (returndummy != 0)
        {
            xl::xout["error"] << "Errors occurred!" << std::endl;
            emit error("elastix()", "Elastix registration error, error code: " + QString::number(returndummy) + ".");
            return returndummy;
        }

        // Get the transform, the fixedImage and the movingImage
        // in order to put it in the (possibly) next registration.

        transform = elastices[i]->GetFinalTransform();
        fixedImageContainer = elastices[i]->GetFixedImageContainer();
        movingImageContainer = elastices[i]->GetMovingImageContainer();
        fixedMaskContainer = elastices[i]->GetFixedMaskContainer();
        movingMaskContainer = elastices[i]->GetMovingMaskContainer();
        fixedImageOriginalDirection = elastices[i]->GetOriginalFixedImageDirectionFlat();

        // Print a finish message.
        elxout << "Running elastix with parameter file " << i
               << ": \"" << argMap["-p"] << "\", has finished.\n" << std::endl;

        // Stop timer and print it.
        timer->StopTimer();
        elxout << "\nCurrent time: " << timer->PrintStopTime() << "." << std::endl;
        elxout << "Time used for running elastix with this parameter file: "
               << timer->PrintElapsedTimeDHMS() << ".\n" << std::endl;

        // Try to release some memory.
        elastices[i] = 0;

    } // end loop over registrations

    elxout << "-------------------------------------------------------------------------" << "\n" << std::endl;

    // Stop totaltimer and print it.
    totaltimer->StopTimer();
    elxout << "Total time elapsed: " << totaltimer->PrintElapsedTimeDHMS() << ".\n" << std::endl;

    // Make sure all the components that are defined in a Module (.DLL/.so)
    // are deleted before the modules are closed.

    for (unsigned int i = 0; i < nrOfParameterFiles; i++)
    {
        elastices[i] = 0;
    }

    transform = 0;
    fixedImageContainer = 0;
    movingImageContainer = 0;
    fixedMaskContainer = 0;
    movingMaskContainer = 0;

    // Close the modules.
    ElastixMainType::UnloadComponents();

    // Remove all unecessary files
    elastixClean(params);

    // Exit and return the error code.
    emit registrationCompleted();
    return returndummy;
}

void milxQtRegistrationAlgos::elastixClean(milxQtRegistrationParams params)
{
    // copy output to the correct path
    // we look for the output
    QString src_basepath = params.outputFolder + "/" + "result.";
    QString src_path;

    int i = -1;
    do
    {
        i++;
        src_path = src_basepath + QString::number(i) +".nii";
    } while (QFile::exists(src_path));

    i--;

    if (i != -1)
    {
        src_path = src_basepath + QString::number(i) + ".nii";
        if (QFile::exists(params.outputName))
        {
            QFile::remove(params.outputName);
        }
        QFile::copy(src_path, params.outputName);
        QFile::remove(src_path);
    }

    // remove file
    QDir directory(params.outputFolder);
    QStringList nameFilter;
    nameFilter << "elastix.log" << "IterationInfo.*" << "TransformParameters.*" << "result.*";
    QStringList fileList = directory.entryList(nameFilter);

    for (i = 0; i < fileList.size(); i++)
    {
        QString filepath = params.outputFolder + "/" + fileList[i];
        if (QFile::exists(filepath))
        {
            QFile::remove(filepath);
        }
    }
}
#endif