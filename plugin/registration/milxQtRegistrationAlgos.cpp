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

#ifdef USE_NIFTI_REG
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
            for (int j = 0; j < nbargs; j++) {
                free(arguments[j]);
            }
            free(arguments);
            emit error("average()", "Error while computing average, the following image can not be read: \"" + QString(arguments[2]) + "\".");

            return EXIT_FAILURE;
        }
        reg_checkAndCorrectDimension(tempImage);

        // Create the average image
        nifti_image *averageImage = nifti_copy_nim_info(tempImage);
        averageImage->scl_slope = 1.f;
        averageImage->scl_inter = 0.f;
        nifti_image_free(tempImage);
        tempImage = NULL;
        averageImage->datatype = NIFTI_TYPE_FLOAT32;
        if (sizeof(PrecisionTYPE) == sizeof(double))
            averageImage->datatype = NIFTI_TYPE_FLOAT64;
        averageImage->nbyper = sizeof(PrecisionTYPE);
        averageImage->data = (void *)calloc(averageImage->nvox, averageImage->nbyper);

        int imageTotalNumber = 0;
        for (int i = 2; i<nbargs; ++i)
        {
            nifti_image *tempImage = reg_io_ReadImageFile(arguments[i]);
            if (tempImage == NULL) {
                //fprintf(stderr, "[!] The following image can not be read: %s\n", arguments[i]);
                for (int j = 0; j < nbargs; j++) {
                    free(arguments[j]);
                }
                free(arguments);
                emit error("average()", "Error while computing average, the following image can not be read: \"" + QString(arguments[i]) + "\".");
                return EXIT_FAILURE;
            }
            reg_checkAndCorrectDimension(tempImage);
            if (sizeof(PrecisionTYPE) == sizeof(double))
                reg_tools_changeDatatype<double>(tempImage);
            else reg_tools_changeDatatype<float>(tempImage);
            if (averageImage->nvox != tempImage->nvox) {
                //fprintf(stderr, "[!] All images must have the same size. Error when processing: %s\n", arguments[i]);
                for (int j = 0; j < nbargs; j++) {
                    free(arguments[j]);
                }
                free(arguments);
                emit error("average()", "Error while computing average, all images must have the same size. Error when processing: \"" + QString(arguments[i]) + "\".");
                return EXIT_FAILURE;
            }
            //            if(sizeof(PrecisionTYPE)==sizeof(double))
            //               average_norm_intensity<double>(tempImage);
            //            else average_norm_intensity<float>(tempImage);
            reg_tools_addImageToImage(averageImage, tempImage, averageImage);
            imageTotalNumber++;
            nifti_image_free(tempImage);
            tempImage = NULL;
        }
        reg_tools_divideValueToImage(averageImage, averageImage, (float)imageTotalNumber);

        reg_io_WriteImageFile(averageImage, arguments[1]);
        nifti_image_free(averageImage);
    }
    else
    {
        // Free variables
        for (int j = 0; j < nbargs; j++) {
            free(arguments[j]);
        }
        free(arguments);
        emit error("average()", "Error while computing average, invalid filetype input (check accepted file extensions: .nii, .nii.gz, .hdr, .img, .img.gz)");
        return EXIT_FAILURE;
    }

    // Free variables
    for (int j = 0; j < nbargs; j++) {
        free(arguments[j]);
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

    // Create some variables
    mat44 *affineTransformation = NULL;
    nifti_image *referenceImage = NULL;
    nifti_image *inputTransformationImage = NULL;
    nifti_image *outputTransformationImage = NULL;
    // First check if the input filename is an image
    if (reg_isAnImageFileName(cppOutputName))
    {
        inputTransformationImage = reg_io_ReadImageFile(cppOutputName);
        if (inputTransformationImage == NULL)
        {
            fprintf(stderr, "[NiftyReg ERROR] Error when reading the provided transformation: %s\n",
                    cppOutputName);
            return 1;
        }
        reg_checkAndCorrectDimension(inputTransformationImage);
        // If the input transformation is a grid, check that the reference image has been specified
        if (inputTransformationImage->intent_p1 == SPLINE_GRID ||
                inputTransformationImage->intent_p1 == SPLINE_VEL_GRID)
        {
            referenceImage = reg_io_ReadImageHeader(referenceName);
            if (referenceImage == NULL)
            {
                fprintf(stderr, "[NiftyReg ERROR] Error when reading the reference image: %s\n",
                        referenceName);
                return 1;
            }
            reg_checkAndCorrectDimension(referenceImage);
        }
    }
    else
    {
        // Read the affine transformation
        affineTransformation = (mat44 *)malloc(sizeof(mat44));
        reg_tool_ReadAffineFile(affineTransformation, cppOutputName);
        referenceImage = reg_io_ReadImageHeader(referenceName);
        if (referenceImage == NULL)
        {
            emit error("cpp2def()", "Error while while computing deformation field, the reference image (\"" + QString(referenceName) + "\") can not be read.");
            return 1;
        }
        reg_checkAndCorrectDimension(referenceImage);
    }
    // Create a dense field
    if (affineTransformation != NULL ||
            inputTransformationImage->intent_p1 == SPLINE_GRID ||
            inputTransformationImage->intent_p1 == SPLINE_VEL_GRID)
    {
        // Create a field image from the reference image
        outputTransformationImage = nifti_copy_nim_info(referenceImage);
        outputTransformationImage->ndim = outputTransformationImage->dim[0] = 5;
        outputTransformationImage->nt = outputTransformationImage->dim[4] = 1;
        outputTransformationImage->nu = outputTransformationImage->dim[5] = outputTransformationImage->nz>1 ? 3 : 2;
        outputTransformationImage->nvox = (size_t)outputTransformationImage->nx *
                                          outputTransformationImage->ny * outputTransformationImage->nz *
                                          outputTransformationImage->nt * outputTransformationImage->nu;
        outputTransformationImage->nbyper = sizeof(float);
        outputTransformationImage->datatype = NIFTI_TYPE_FLOAT32;
        outputTransformationImage->intent_code = NIFTI_INTENT_VECTOR;
        memset(outputTransformationImage->intent_name, 0, 16);
        strcpy(outputTransformationImage->intent_name, "NREG_TRANS");
        outputTransformationImage->scl_slope = 1.f;
        outputTransformationImage->scl_inter = 0.f;
    }
    else
    {
        // Create a deformation field from in the input transformation
        outputTransformationImage = nifti_copy_nim_info(inputTransformationImage);
    }
    // Allocate the output field data array
    outputTransformationImage->data = (void *)malloc
                                      (outputTransformationImage->nvox*outputTransformationImage->nbyper);

    // Create a deformation or displacement field
    if (affineTransformation != NULL)
    {
        reg_affine_getDeformationField(affineTransformation, outputTransformationImage);
    }
    else
    {
        switch (static_cast<int>(reg_round(inputTransformationImage->intent_p1)))
        {
        case DEF_FIELD:
            //printf("[NiftyReg] The specified transformation is a deformation field:\n[NiftyReg] %s\n", inputTransformationImage->fname);
            // the current in transformation is copied
            memcpy(outputTransformationImage->data, inputTransformationImage->data,
                   outputTransformationImage->nvox*outputTransformationImage->nbyper);
            break;
        case DISP_FIELD:
            //printf("[NiftyReg] The specified transformation is a displacement field:\n[NiftyReg] %s\n", inputTransformationImage->fname);
            // the current in transformation is copied and converted
            memcpy(outputTransformationImage->data, inputTransformationImage->data,
                   outputTransformationImage->nvox*outputTransformationImage->nbyper);
            reg_getDeformationFromDisplacement(outputTransformationImage);
            break;
        case SPLINE_GRID:
            // printf("[NiftyReg] The specified transformation is a spline parametrisation:\n[NiftyReg] %s\n", inputTransformationImage->fname);
            // The output field is filled with an identity deformation field
            memset(outputTransformationImage->data,
                   0,
                   outputTransformationImage->nvox*outputTransformationImage->nbyper);
            reg_getDeformationFromDisplacement(outputTransformationImage);
            // The spline transformation is composed with the identity field
            reg_spline_getDeformationField(inputTransformationImage,
                                           outputTransformationImage,
                                           NULL, // no mask
                                           true, // composition is used,
                                           true // b-spline are used
                                          );
            break;
        case DEF_VEL_FIELD:
            //printf("[NiftyReg] The specified transformation is a deformation velocity field:\n[NiftyReg] %s\n", inputTransformationImage->fname);
            // The flow field is exponentiated
            reg_defField_getDeformationFieldFromFlowField(inputTransformationImage,
                    outputTransformationImage,
                    false // step number is not updated
                                                         );
            break;
        case DISP_VEL_FIELD:
            //printf("[NiftyReg] The specified transformation is a displacement velocity field:\n[NiftyReg] %s\n", inputTransformationImage->fname);
            // The input transformation is converted into a def flow
            reg_getDeformationFromDisplacement(outputTransformationImage);
            // The flow field is exponentiated
            reg_defField_getDeformationFieldFromFlowField(inputTransformationImage,
                    outputTransformationImage,
                    false // step number is not updated
                                                         );
            break;
        case SPLINE_VEL_GRID:
            //printf("[NiftyReg] The specified transformation is a spline velocity parametrisation:\n[NiftyReg] %s\n", inputTransformationImage->fname);
            // The spline parametrisation is converted into a dense flow and exponentiated
            reg_spline_getDefFieldFromVelocityGrid(inputTransformationImage,
                                                   outputTransformationImage,
                                                   false // step number is not updated
                                                  );
            break;
        default:
            emit error("cpp2def()", "[NiftyReg ERROR] Unknown input transformation type");
            return 1;
        }
    }
    outputTransformationImage->intent_p1 = DEF_FIELD;
    outputTransformationImage->intent_p2 = 0;


    // Save the generated transformation
    reg_io_WriteImageFile(outputTransformationImage, defOutputName);

    // Free the allocated images and arrays
    if (affineTransformation != NULL) free(affineTransformation);
    if (referenceImage != NULL) nifti_image_free(referenceImage);
    if (inputTransformationImage != NULL) nifti_image_free(inputTransformationImage);
    if (outputTransformationImage != NULL) nifti_image_free(outputTransformationImage);
    emit cpp2defCompleted();
    return 0;


}

int milxQtRegistrationAlgos::f3d(milxQtRegistrationParams params)
{
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

    time_t start;
    time(&start);
    int verbose = true;

    nifti_image *referenceImage = NULL;
    nifti_image *floatingImage = NULL;

    referenceImage = reg_io_ReadImageFile(referenceName);
    if (referenceImage == NULL)
    {
        emit error("f3d()", "F3D registration error: Error when reading the reference image \"" + QString(referenceName) + "\".");
        return 1;
    }

    floatingImage = reg_io_ReadImageFile(floatingName);
    if (floatingImage == NULL)
    {
        emit error("f3d()", "F3D registration error: Error when reading the floating image \"" + QString(floatingName) + "\".");
        return 1;
    }

    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    // Check the type of registration object to create
#ifdef _USE_CUDA
    CUcontext ctx;
#endif // _USE_CUDA
    reg_f3d<PrecisionTYPE> *REG = NULL;

    if (useSym)
    {
        REG = new reg_f3d_sym<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
    }

    if (REG == NULL)
        REG = new reg_f3d<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
    REG->SetReferenceImage(referenceImage);
    REG->SetFloatingImage(floatingImage);

    // Create some pointers that could be used
    mat44 affineMatrix;
    nifti_image *inputCCPImage = NULL;
    nifti_image *referenceMaskImage = NULL;
    nifti_image *floatingMaskImage = NULL;

    int refBinNumber = 0;
    int floBinNumber = 0;

    /* read the input parameter */
    // Verbose off
    REG->DoNotPrintOutInformation();

    REG->SetMaximalIterationNumber(maxiterationNumber);
    REG->SetSpacing(0, spacing[0]);
    REG->SetSpacing(1, spacing[1]);
    REG->SetSpacing(2, spacing[2]);
    REG->SetLevelNumber(levelNumber);
    REG->SetLevelToPerform(levelToPerform);
    if (noPyramid)
    {
        REG->DoNotUsePyramidalApproach();
    }

    // Run the registration
    REG->Run();

    // Save the control point result
    nifti_image *outputControlPointGridImage = REG->GetControlPointPositionImage();
    if (cppOutputName[0] == '\0') strcpy(cppOutputName, "outputCPP.nii.gz");
    memset(outputControlPointGridImage->descrip, 0, 80);
    strcpy(outputControlPointGridImage->descrip, "Control point position from NiftyReg (reg_f3d)");
    if (strcmp("NiftyReg F3D2", REG->GetExecutableName()) == 0)
        strcpy(outputControlPointGridImage->descrip, "Velocity field grid from NiftyReg (reg_f3d2)");
    reg_io_WriteImageFile(outputControlPointGridImage, cppOutputName);
    nifti_image_free(outputControlPointGridImage);
    outputControlPointGridImage = NULL;

    // Save the backward control point result
    if (REG->GetSymmetricStatus())
    {
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
        else if (b.find(".png") != std::string::npos)
            b.replace(b.find(".png"), 4, "_backward.png");
        else if (b.find(".nrrd") != std::string::npos)
            b.replace(b.find(".nrrd"), 5, "_backward.nrrd");
        else b.append("_backward.nii");
        nifti_image *outputBackwardControlPointGridImage = REG->GetBackwardControlPointPositionImage();
        memset(outputBackwardControlPointGridImage->descrip, 0, 80);
        strcpy(outputBackwardControlPointGridImage->descrip, "Backward Control point position from NiftyReg (reg_f3d)");
        if (strcmp("NiftyReg F3D2", REG->GetExecutableName()) == 0)
            strcpy(outputBackwardControlPointGridImage->descrip, "Backward velocity field grid from NiftyReg (reg_f3d2)");
        reg_io_WriteImageFile(outputBackwardControlPointGridImage, b.c_str());
        nifti_image_free(outputBackwardControlPointGridImage);
        outputBackwardControlPointGridImage = NULL;
    }

    // Save the warped image result(s)
    nifti_image **outputWarpedImage = (nifti_image **)malloc(2 * sizeof(nifti_image *));
    outputWarpedImage[0] = NULL;
    outputWarpedImage[1] = NULL;
    outputWarpedImage = REG->GetWarpedImage();
    if (outputName[0] == '\0') strcpy(outputName, "outputResult.nii");

    memset(outputWarpedImage[0]->descrip, 0, 80);
    strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d)");
    if (strcmp("NiftyReg F3D SYM", REG->GetExecutableName()) == 0)
    {
        strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
        strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
    }
    if (strcmp("NiftyReg F3D2", REG->GetExecutableName()) == 0)
    {
        strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d2)");
        strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d2)");
    }
    if (REG->GetSymmetricStatus())
    {
        if (outputWarpedImage[1] != NULL)
        {
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
            else if (b.find(".png") != std::string::npos)
                b.replace(b.find(".png"), 4, "_backward.png");
            else if (b.find(".nrrd") != std::string::npos)
                b.replace(b.find(".nrrd"), 5, "_backward.nrrd");
            else b.append("_backward.nii");
            reg_io_WriteImageFile(outputWarpedImage[1], b.c_str());
        }
    }
    reg_io_WriteImageFile(outputWarpedImage[0], outputName);
    if (outputWarpedImage[0] != NULL)
        nifti_image_free(outputWarpedImage[0]);
    outputWarpedImage[0] = NULL;
    if (outputWarpedImage[1] != NULL)
        nifti_image_free(outputWarpedImage[1]);
    outputWarpedImage[1] = NULL;
    free(outputWarpedImage);
    outputWarpedImage = NULL;
#ifdef _USE_CUDA
    cudaCommon_unsetCUDACard(&ctx);
#endif
    // Erase the registration object
    delete REG;

    // Clean the allocated images
    if (referenceImage != NULL) nifti_image_free(referenceImage);
    if (floatingImage != NULL) nifti_image_free(floatingImage);

    time_t end;
    time(&end);
    int minutes = (int)floorf(float(end - start) / 60.0f);
    int seconds = (int)(end - start - 60 * minutes);

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
    if (!outputResultFlag) strcpy(outputName, "outputResult.nii.gz");
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

    milxQtSimilarities similarities;

    // We compute the similarity before the registration (i == 0) and after the registration (i == 1)
    for (int i = 0; i < 2; i++)
    {
        QString similarityOutput = image->createFile(QDir::tempPath() + QDir::separator() + "similarity_XXXXXX.txt");
        char referenceName[FILENAME_MAX + 1],
             floatingName[FILENAME_MAX + 1],
             outputName[FILENAME_MAX + 1];

        strncpy(referenceName, image->reference->getPath().toLatin1().constData(), FILENAME_MAX);

        if (i == 0)
            strncpy(floatingName, image->getPath().toLatin1().constData(), FILENAME_MAX);
        else if (i == 1)
            strncpy(floatingName, image->getOutputPath().toLatin1().constData(), FILENAME_MAX);

        strncpy(outputName, similarityOutput.toLatin1().constData(), FILENAME_MAX);


        /* Read the reference image */
        nifti_image *refImage = reg_io_ReadImageFile(referenceName);
        if (refImage == NULL)
        {
            emit error("similarity()", "Similarity measure error: Error when reading the reference image (\"" + QString(referenceName) + "\").");
            return 1;
        }
        reg_checkAndCorrectDimension(refImage);
        reg_tools_changeDatatype<float>(refImage);

        /* Read the floating image */
        nifti_image *floImage = reg_io_ReadImageFile(floatingName);
        if (floImage == NULL)
        {
            emit error("similarity()", "Similarity measure error: Error when reading the floating image (\"" + QString(floatingName) + "\").");
            return 1;
        }
        reg_checkAndCorrectDimension(floImage);
        reg_tools_changeDatatype<float>(floImage);

        /* Read and create the mask array */
        int *refMask = NULL;
        int refMaskVoxNumber = refImage->nx*refImage->ny*refImage->nz;
        refMask = (int *)calloc(refMaskVoxNumber, sizeof(int));
        for (int j = 0; j<refMaskVoxNumber; ++j) refMask[j] = j;


        /* Create the warped floating image */
        nifti_image *warpedFloImage = nifti_copy_nim_info(refImage);
        warpedFloImage->ndim = warpedFloImage->dim[0] = floImage->ndim;
        warpedFloImage->nt = warpedFloImage->dim[4] = floImage->nt;
        warpedFloImage->nu = warpedFloImage->dim[5] = floImage->nu;
        warpedFloImage->nvox = (size_t)warpedFloImage->nx * warpedFloImage->ny *
                               warpedFloImage->nz * warpedFloImage->nt * warpedFloImage->nu;
        warpedFloImage->cal_min = floImage->cal_min;
        warpedFloImage->cal_max = floImage->cal_max;
        warpedFloImage->scl_inter = floImage->scl_inter;
        warpedFloImage->scl_slope = floImage->scl_slope;
        warpedFloImage->datatype = floImage->datatype;
        warpedFloImage->nbyper = floImage->nbyper;
        warpedFloImage->data = (void *)malloc(warpedFloImage->nvox*warpedFloImage->nbyper);

        /* Create the deformation field */
        nifti_image *defField = nifti_copy_nim_info(refImage);
        defField->ndim = defField->dim[0] = 5;
        defField->nt = defField->dim[4] = 1;
        defField->nu = defField->dim[5] = refImage->nz>1 ? 3 : 2;
        defField->nvox = (size_t)defField->nx * defField->ny *
                         defField->nz * defField->nt * defField->nu;
        defField->datatype = NIFTI_TYPE_FLOAT32;
        defField->nbyper = sizeof(float);
        defField->data = (void *)calloc(defField->nvox, defField->nbyper);
        defField->scl_slope = 1.f;
        defField->scl_inter = 0.f;
        reg_tools_multiplyValueToImage(defField, defField, 0.f);
        defField->intent_p1 = DISP_FIELD;
        reg_getDeformationFromDisplacement(defField);

        /* Warp the floating image */
        reg_resampleImage(floImage,
                          warpedFloImage,
                          defField,
                          refMask,
                          3,
                          std::numeric_limits<float>::quiet_NaN());
        nifti_image_free(defField);

        FILE *outFile = NULL;
        outFile = fopen(outputName, "w");

        /* Compute the NCC if required */
        {
            float *refPtr = static_cast<float *>(refImage->data);
            float *warPtr = static_cast<float *>(warpedFloImage->data);
            double refMeanValue = 0.;
            double warMeanValue = 0.;
            refMaskVoxNumber = 0;
            for (size_t i = 0; i<refImage->nvox; ++i) {
                if (refMask[i]>-1 && refPtr[i] == refPtr[i] && warPtr[i] == warPtr[i]) {
                    refMeanValue += refPtr[i];
                    warMeanValue += warPtr[i];
                    ++refMaskVoxNumber;
                }
            }
            if (refMaskVoxNumber == 0)
                fprintf(stderr, "No active voxel\n");
            refMeanValue /= (double)refMaskVoxNumber;
            warMeanValue /= (double)refMaskVoxNumber;
            double refSTDValue = 0.;
            double warSTDValue = 0.;
            double measure = 0.;
            for (size_t i = 0; i<refImage->nvox; ++i) {
                if (refMask[i]>-1 && refPtr[i] == refPtr[i] && warPtr[i] == warPtr[i]) {
                    refSTDValue += reg_pow2((double)refPtr[i] - refMeanValue);
                    warSTDValue += reg_pow2((double)warPtr[i] - warMeanValue);
                    measure += ((double)refPtr[i] - refMeanValue) *
                               ((double)warPtr[i] - warMeanValue);
                }
            }
            refSTDValue /= (double)refMaskVoxNumber;
            warSTDValue /= (double)refMaskVoxNumber;
            measure /= sqrt(refSTDValue)*sqrt(warSTDValue)*
                       (double)refMaskVoxNumber;
            if (outFile != NULL)
                fprintf(outFile, "%g\n", measure);
            else printf("NCC: %g\n", measure);
        }
        /* Compute the LNCC if required */
        {
            reg_lncc *lncc_object = new reg_lncc();
            for (int j = 0; j < (refImage->nt < warpedFloImage->nt ? refImage->nt : warpedFloImage->nt); ++j)
                lncc_object->SetActiveTimepoint(j);
            lncc_object->InitialiseMeasure(refImage,
                                           warpedFloImage,
                                           refMask,
                                           warpedFloImage,
                                           NULL,
                                           NULL);
            double measure = lncc_object->GetSimilarityMeasureValue();
            if (outFile != NULL)
                fprintf(outFile, "%g\n", measure);
            else printf("LNCC: %g\n", measure);
            delete lncc_object;
        }
        /* Compute the NMI if required */
        {
            reg_nmi *nmi_object = new reg_nmi();
            for (int j = 0; j < (refImage->nt < warpedFloImage->nt ? refImage->nt : warpedFloImage->nt); ++j)
                nmi_object->SetActiveTimepoint(i);
            nmi_object->InitialiseMeasure(refImage,
                                          warpedFloImage,
                                          refMask,
                                          warpedFloImage,
                                          NULL,
                                          NULL);
            double measure = nmi_object->GetSimilarityMeasureValue();
            if (outFile != NULL)
                fprintf(outFile, "%g\n", measure);
            else printf("NMI: %g\n", measure);
            delete nmi_object;
        }
        /* Compute the SSD if required */
        {
            reg_ssd *ssd_object = new reg_ssd();
            for (int j = 0; j < (refImage->nt < warpedFloImage->nt ? refImage->nt : warpedFloImage->nt); ++j)
                ssd_object->SetActiveTimepoint(i);
            ssd_object->InitialiseMeasure(refImage,
                                          warpedFloImage,
                                          refMask,
                                          warpedFloImage,
                                          NULL,
                                          NULL);
            double measure = ssd_object->GetSimilarityMeasureValue();
            if (outFile != NULL)
                fprintf(outFile, "%g\n", measure);
            else printf("SSD: %g\n", measure);
            delete ssd_object;
        }
        // Close the output file if required
        if (outFile != NULL)
            fclose(outFile);

        // Free the allocated images
        nifti_image_free(refImage);
        nifti_image_free(floImage);
        free(refMask);

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
