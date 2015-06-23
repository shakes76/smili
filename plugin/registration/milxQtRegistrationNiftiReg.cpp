#include "milxQtRegistrationNiftiReg.h"


milxQtRegistrationNifti::milxQtRegistrationNifti()
{

}

milxQtRegistrationNifti::~milxQtRegistrationNifti()
{

}

void milxQtRegistrationNifti::cpp2def_async(PARAMSCPP2DEF params)
{
	future = QtConcurrent::run(this, &milxQtRegistrationNifti::cpp2def, params);
}

void milxQtRegistrationNifti::f3d_async(PARAMSF3D params)
{
	future = QtConcurrent::run(this, &milxQtRegistrationNifti::f3d, params);
}

void milxQtRegistrationNifti::aladin_async(PARAMSALADIN params)
{
	future = QtConcurrent::run(this, &milxQtRegistrationNifti::aladin, params);
}

int milxQtRegistrationNifti::cpp2def(PARAMSCPP2DEF params)
{
	/* Read the reference image */
	nifti_image *referenceImage = reg_io_ReadImageHeader(params.referenceImageName);
	if (referenceImage == NULL){
		emit cpp2defFinished();
		return 1;
	}
	reg_checkAndCorrectDimension(referenceImage);


	// Read the control point image
	nifti_image *controlPointImage = reg_io_ReadImageFile(params.cpp2defInputName);
    if (controlPointImage == NULL){
		emit cpp2defFinished();
		return 1;
	}
	reg_checkAndCorrectDimension(controlPointImage);
	// Allocate the deformation field
	nifti_image *deformationFieldImage = nifti_copy_nim_info(referenceImage);
	deformationFieldImage->dim[0] = deformationFieldImage->ndim = 5;
	deformationFieldImage->dim[1] = deformationFieldImage->nx = referenceImage->nx;
	deformationFieldImage->dim[2] = deformationFieldImage->ny = referenceImage->ny;
	deformationFieldImage->dim[3] = deformationFieldImage->nz = referenceImage->nz;
	deformationFieldImage->dim[4] = deformationFieldImage->nt = 1; deformationFieldImage->pixdim[4] = deformationFieldImage->dt = 1.0;
	if (referenceImage->nz>1) deformationFieldImage->dim[5] = deformationFieldImage->nu = 3;
	else deformationFieldImage->dim[5] = deformationFieldImage->nu = 2;
	deformationFieldImage->pixdim[5] = deformationFieldImage->du = 1.0;
	deformationFieldImage->dim[6] = deformationFieldImage->nv = 1; deformationFieldImage->pixdim[6] = deformationFieldImage->dv = 1.0;
	deformationFieldImage->dim[7] = deformationFieldImage->nw = 1; deformationFieldImage->pixdim[7] = deformationFieldImage->dw = 1.0;
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
	reg_io_WriteImageFile(deformationFieldImage, params.cpp2defOutputName);
	nifti_image_free(deformationFieldImage);


	nifti_image_free(referenceImage);
	emit cpp2defFinished();
	return 0;
}

int milxQtRegistrationNifti::f3d(PARAMSF3D params)
{
	time_t start; time(&start);
	
	char *referenceName = params.referenceName;
	char *floatingName = params.floatingName;
	char * outputControlPointGridName = params.outputControlPointGridName;
	char * outputWarpedName = params.outputWarpedName;
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
	for (int i = 0; i<10; i++){
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


	if (referenceName == NULL || floatingName == NULL){
		//fprintf(stderr, "Err:\tThe reference and the floating image have to be defined.\n");
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
	if (verbose == true){
#endif
		/*
		printf("\n[NiftyReg F3D] Command line:\n\t");
		for (int i = 0; i<argc; i++)
			printf(" %s", argv[i]);
		printf("\n\n");
		*/
#ifdef NDEBUG
	}
#endif

	// Read the reference image
	if (referenceName == NULL){
		//fprintf(stderr, "Error. No reference image has been defined\n");
		return 1;
	}
	nifti_image *referenceImage = reg_io_ReadImageFile(referenceName);
	if (referenceImage == NULL){
		//fprintf(stderr, "Error when reading the reference image %s\n", referenceName);
		return 1;
	}

	// Read the floating image
	if (floatingName == NULL){
		//fprintf(stderr, "Error. No floating image has been defined\n");
		return 1;
	}
	nifti_image *floatingImage = reg_io_ReadImageFile(floatingName);
	if (floatingImage == NULL){
		//fprintf(stderr, "Error when reading the floating image %s\n", floatingName);
		return 1;
	}

	// Read the mask images
	nifti_image *referenceMaskImage = NULL;
	if (referenceMaskName != NULL){
		referenceMaskImage = reg_io_ReadImageFile(referenceMaskName);
		if (referenceMaskImage == NULL){
			//fprintf(stderr, "Error when reading the reference mask image %s\n", referenceMaskName);
			return 1;
		}
	}
	nifti_image *floatingMaskImage = NULL;
	if (floatingMaskName != NULL){
		floatingMaskImage = reg_io_ReadImageFile(floatingMaskName);
		if (floatingMaskImage == NULL){
			//fprintf(stderr, "Error when reading the reference mask image %s\n", floatingMaskName);
			return 1;
		}
	}

	// Read the input control point grid image
	nifti_image *controlPointGridImage = NULL;
	if (inputControlPointGridName != NULL){
		controlPointGridImage = reg_io_ReadImageFile(inputControlPointGridName);
		if (controlPointGridImage == NULL){
			//fprintf(stderr, "Error when reading the input control point grid image %s\n", inputControlPointGridName);
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
	if (affineTransformationName != NULL){
		affineTransformation = (mat44 *)malloc(sizeof(mat44));
		// Check first if the specified affine file exist
		if (FILE *aff = fopen(affineTransformationName, "r")){
			fclose(aff);
		}
		else{
			//fprintf(stderr, "The specified input affine file (%s) can not be read\n", affineTransformationName);
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
	if (useGPU){

		if (linearEnergyWeight0 == linearEnergyWeight0 ||
			linearEnergyWeight1 == linearEnergyWeight1 ||
			L2NormWeight == L2NormWeight){
			//fprintf(stderr, "NiftyReg ERROR CUDA] The linear elasticity has not been implemented with CUDA yet. Exit.\n");
			exit(0);
		}

		if (useSym){
			//fprintf(stderr, "\n[NiftyReg ERROR CUDA] GPU implementation of the symmetric registration is not available yet. Exit\n");
			exit(0);
		}
#ifdef _BUILD_NR_DEV
		if (useVel){
			//fprintf(stderr, "\n[NiftyReg ERROR CUDA] GPU implementation of velocity field parametrisartion is not available yet. Exit\n");
			exit(0);
		}
#endif

		if ((referenceImage->dim[4] == 1 && floatingImage->dim[4] == 1) || (referenceImage->dim[4] == 2 && floatingImage->dim[4] == 2)){

			// The CUDA card is setup
			cuInit(0);
			struct cudaDeviceProp deviceProp;
			int device_count = 0;
			cudaGetDeviceCount(&device_count);
			int device = cardNumber;
			if(cardNumber==-1){
				// following code is from cutGetMaxGflopsDeviceId()
				int max_gflops_device = 0;
				int max_gflops = 0;
				int current_device = 0;
				while( current_device < device_count ){
					cudaGetDeviceProperties( &deviceProp, current_device );
					int gflops = deviceProp.multiProcessorCount * deviceProp.clockRate;
					if (gflops > max_gflops){
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
			if (deviceProp.major < 1){
				//printf("[NiftyReg ERROR CUDA] The specified graphical card does not exist.\n");
				return 1;
			}
			REG = new reg_f3d_gpu<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
			if (verbose == true){
#endif
				//printf("\n[NiftyReg F3D] GPU implementation is used\n");
#ifdef NDEBUG
			}
#endif
}
		else{
			//fprintf(stderr, "[NiftyReg ERROR] The GPU implementation only handle 1 to 1 or 2 to 2 image(s) registration\n");
			exit(1);
		}
	}
	else
#endif
	{
		if(useSym){
			REG = new reg_f3d_sym<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
			if (verbose == true){
#endif // NDEBUG
				//printf("\n[NiftyReg F3D SYM] CPU implementation is used\n");
#ifdef NDEBUG
			}
#endif // NDEBUG
	}
#ifdef _BUILD_NR_DEV
		else if (useVel){
			REG = new reg_f3d2<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
			if(verbose==true){
#endif
				//printf("\n[NiftyReg F3D2] CPU implementation is used\n");
#ifdef NDEBUG
		}
#endif
		}
#endif // _BUILD_NR_DEV
		else{
			REG = new reg_f3d<PrecisionTYPE>(referenceImage->nt, floatingImage->nt);
#ifdef NDEBUG
			if (verbose == true){
#endif // NDEBUG
				//printf("\n[NiftyReg F3D] CPU implementation is used\n");
#ifdef NDEBUG
			}
#endif // NDEBUG
		}
	}
#ifdef _OPENMP
	int maxThreadNumber = omp_get_max_threads();
#ifdef NDEBUG
	if(verbose==true){
#endif // NDEBUG
		printf("[NiftyReg F3D] OpenMP is used with %i thread(s)\n", maxThreadNumber);
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

	if (linearEnergyWeight0 == linearEnergyWeight0 || linearEnergyWeight1 == linearEnergyWeight1){
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

	for(unsigned int t=0;t<(unsigned int)referenceImage->nt;t++)
		if (referenceThresholdLow[t] == referenceThresholdLow[t])
			REG->SetReferenceThresholdLow(t, referenceThresholdLow[t]);

	for(unsigned int t=0;t<(unsigned int)floatingImage->nt;t++)
		if(floatingThresholdUp[t]==floatingThresholdUp[t])
			REG->SetFloatingThresholdUp(t, floatingThresholdUp[t]);

	for(unsigned int t=0;t<(unsigned int)floatingImage->nt;t++)
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

	switch (interpolation){
	case 0:REG->UseNeareatNeighborInterpolation();
		break;
	case 1:REG->UseLinearInterpolation();
		break;
	case 3:REG->UseCubicSplineInterpolation();
		break;
	default:REG->UseLinearInterpolation();
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

	// F3D SYM arguments
	if (floatingMaskImage != NULL){
		if (useSym)
			REG->SetFloatingMask(floatingMaskImage);
#ifdef _BUILD_NR_DEV
		else if (useVel)
			REG->SetFloatingMask(floatingMaskImage);
#endif
		else{
			//fprintf(stderr, "[NiftyReg WARNING] The floating mask image is only used for symmetric or velocity field parametrisation\n");
			//fprintf(stderr, "[NiftyReg WARNING] The floating mask image is ignored\n");
		}
	}

	if (inverseConsistencyWeight == inverseConsistencyWeight)
		REG->SetInverseConsistencyWeight(inverseConsistencyWeight);

#ifdef _BUILD_NR_DEV
	// F3D2 arguments
	if (stepNumber>-1)
		REG->SetCompositionStepNumber(stepNumber);
#endif

	// Run the registration
#ifdef _USE_CUDA
	if (useGPU && checkMem){
		size_t free, total, requiredMemory = REG->CheckMemoryMB_f3d();
		cuMemGetInfo(&free, &total);
		//printf("[NiftyReg CUDA] The required memory to run the registration is %lu Mb\n",
			(unsigned long int)requiredMemory);
		//printf("[NiftyReg CUDA] The GPU card has %lu Mb from which %lu Mb are currenlty free\n",
			(unsigned long int)total / (1024 * 1024), (unsigned long int)free / (1024 * 1024));
	}
	else{
#endif
		REG->Run_f3d();

		// Save the control point result
		nifti_image *outputControlPointGridImage = REG->GetControlPointPositionImage();
		if (outputControlPointGridName == NULL) outputControlPointGridName = (char *)"outputCPP.nii";
		memset(outputControlPointGridImage->descrip, 0, 80);
		strcpy(outputControlPointGridImage->descrip, "Control point position from NiftyReg (reg_f3d)");
#ifdef _BUILD_NR_DEV
		if (useVel)
			strcpy(outputControlPointGridImage->descrip, "Velocity field grid from NiftyReg (reg_f3d2)");
#endif
		reg_io_WriteImageFile(outputControlPointGridImage, outputControlPointGridName);
		nifti_image_free(outputControlPointGridImage); outputControlPointGridImage = NULL;

		// Save the backward control point result
#ifdef _BUILD_NR_DEV
		if (useSym || useVel){
#else
		if (useSym){
#endif
			// _backward is added to the forward control point grid image name
			std::string b(outputControlPointGridName);
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
			nifti_image_free(outputBackwardControlPointGridImage); outputBackwardControlPointGridImage = NULL;
		}

		// Save the warped image result(s)
		nifti_image **outputWarpedImage = (nifti_image **)malloc(2 * sizeof(nifti_image *));
		outputWarpedImage[0] = outputWarpedImage[1] = NULL;
		outputWarpedImage = REG->GetWarpedImage();
		if (outputWarpedName == NULL) outputWarpedName = (char *)"outputResult.nii";
		memset(outputWarpedImage[0]->descrip, 0, 80);
		strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d)");
		if (useSym){
			strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
			strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d_sym)");
#ifdef _BUILD_NR_DEV
		}
		if (useVel){
			strcpy(outputWarpedImage[0]->descrip, "Warped image using NiftyReg (reg_f3d2)");
			strcpy(outputWarpedImage[1]->descrip, "Warped image using NiftyReg (reg_f3d2)");
		}
		if (useSym || useVel){
#endif
			if (outputWarpedImage[1] != NULL){
				std::string b(outputWarpedName);
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
				nifti_image_free(outputWarpedImage[1]); outputWarpedImage[1] = NULL;
			}
		}
		reg_io_WriteImageFile(outputWarpedImage[0], outputWarpedName);
		nifti_image_free(outputWarpedImage[0]); outputWarpedImage[0] = NULL;
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
	if (verbose){
#endif
		time_t end; time(&end);
		int minutes = (int)floorf(float(end - start) / 60.0f);
		int seconds = (int)(end - start - 60 * minutes);

#ifdef _USE_CUDA
		if (!checkMem){
#endif
			//printf("[NiftyReg F3D] Registration Performed in %i min %i sec\n", minutes, seconds);
			//printf("[NiftyReg F3D] Have a good day !\n");
#ifdef _USE_CUDA
		}
#endif
#ifdef NDEBUG
	}
#endif

emit registrationFinished();
return 0;
}

int milxQtRegistrationNifti::aladin(PARAMSALADIN params)
{
	time_t start; time(&start);

	char *referenceImageName = params.referenceName;
	int referenceImageFlag = 1;

	char *floatingImageName = params.floatingName;
	int floatingImageFlag = 1;

	char *outputResultName = params.outputResultName;
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

	if (params.aF3Direct)
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

	if (!referenceImageFlag || !floatingImageFlag){
		//fprintf(stderr, "Err:\tThe reference and the floating image have to be defined.\n");
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
	nifti_image *referenceHeader = reg_io_ReadImageFile(referenceImageName);
	if (referenceHeader == NULL){
		//fprintf(stderr, "* ERROR Error when reading the reference  image: %s\n", referenceImageName);
		return 1;
	}

	/* Read teh floating image and check its dimension */
	nifti_image *floatingHeader = reg_io_ReadImageFile(floatingImageName);
	if (floatingHeader == NULL){
		//fprintf(stderr, "* ERROR Error when reading the floating image: %s\n", floatingImageName);
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
	if (referenceMaskFlag){
		referenceMaskImage = reg_io_ReadImageFile(referenceMaskName);
		if (referenceMaskImage == NULL){
			//fprintf(stderr, "* ERROR Error when reading the reference mask image: %s\n", referenceMaskName);
			return 1;
		}
		/* check the dimension */
		for (int i = 1; i <= referenceHeader->dim[0]; i++){
			if (referenceHeader->dim[i] != referenceMaskImage->dim[i]){
				//fprintf(stderr, "* ERROR The reference image and its mask do not have the same dimension\n");
				return 1;
			}
		}
		REG->SetInputMask(referenceMaskImage);
	}
#ifdef _BUILD_NR_DEV
	nifti_image *floatingMaskImage = NULL;
	if (floatingMaskFlag && symFlag){
		floatingMaskImage = reg_io_ReadImageFile(floatingMaskName);
		if (floatingMaskImage == NULL){
			//fprintf(stderr, "* ERROR Error when reading the floating mask image: %s\n", referenceMaskName);
			return 1;
		}
		/* check the dimension */
		for (int i = 1; i <= floatingHeader->dim[0]; i++){
			if (floatingHeader->dim[i] != floatingMaskImage->dim[i]){
				//fprintf(stderr, "* ERROR The floating image and its mask do not have the same dimension\n");
				return 1;
			}
		}
		REG->SetInputFloatingMask(floatingMaskImage);
	}
#endif
	REG->Run();

	// The warped image is saved
	nifti_image *outputResultImage = REG->GetFinalWarpedImage();
	if (!outputResultFlag) outputResultName = (char *)"outputResult.nii";
	reg_io_WriteImageFile(outputResultImage, outputResultName);
	nifti_image_free(outputResultImage);

	/* The affine transformation is saved */
	if (outputAffineFlag)
		reg_tool_WriteAffineFile(REG->GetTransformationMatrix(), outputAffineName);
	else reg_tool_WriteAffineFile(REG->GetTransformationMatrix(), (char *)"outputAffine.txt");

	nifti_image_free(referenceHeader);
	nifti_image_free(floatingHeader);

	delete REG;
	time_t end; time(&end);
	int minutes = (int)floorf((end - start) / 60.0f);
	int seconds = (int)(end - start - 60 * minutes);


	emit registrationFinished();
	return 0;
}
