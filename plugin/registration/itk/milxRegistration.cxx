/*=========================================================================
Program: MILX MixView
Module: milxRegistration.cxx
Author: Jurgen Fripp
Modified by:
Language: C++
Created: Tue 24 Mar 2009 16:25:48 EST

Copyright: (c) 2009 CSIRO, Australia.

This software is protected by international copyright laws.
Any unauthorised copying, distribution or reverse engineering is prohibited.

Licence:
All rights in this Software are reserved to CSIRO. You are only permitted
to have this Software in your possession and to make use of it if you have
agreed to a Software License with CSIRO.

BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

// #define USE_BOOST
//TODO ND: this should really be set by cmake

//TODO Set Affine transform to initialise NR registration (method missing in ITK 4.0.0)

#include "milxRegistration.h"

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkPermuteAxesImageFilter.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include "itkHistogramMatchingImageFilter.h"
#include "itkMultiResolutionPDEDeformableRegistration.h"

#include <itkImageRegionIteratorWithIndex.h>

#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>

#include <itkIdentityTransform.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCenteredTransformInitializer.h>
#include <itkVersorRigid3DTransform.h>
#include <itkVersorRigid3DTransformOptimizer.h>
#include <itkAffineTransform.h>
#include <itkEuler3DTransform.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineResampleImageFunction.h>
#include <itkBSplineDecompositionImageFilter.h>
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include "itkImageMaskSpatialObject.h"
#include <itkThresholdImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>

#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFileWriter.h>

#if ITK_VERSION_MAJOR < 4
#include <itkOrientedImage.h>
#include <itkBSplineDeformableTransform.h>
#include <itkCompose3DVectorImageFilter.h>
#else
#include <itkImage.h>
#include <itkBSplineTransform.h>
#include <itkComposeImageFilter.h>
#endif

#include <cmath>

#include <itkLinearInterpolateImageFunction.h>
#include <itkImageRegistrationMethod.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresHistogramImageToImageMetric.h"
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkCenteredVersorTransformInitializer.h>
#include <itkSpatialObjectToImageFilter.h>

#include <fstream>
#include <iostream>

//!ND
//#define USE_BOOST
#ifdef USE_BOOST
#include "itkAliBabaPyramidWrapper.h"
#endif

namespace milx {


/* STAPLE Filter files */
#define USE_EXTERNAL_CALL2VOTE 1
#define STAPLE_EXEC "itkMultiLabelSTAPLEImageFilter"
#define STAPLE_TYPE "VOTE"
#define STAPLE_TYPE2 "MULTISTAPLE"

#define MILXBSPLINE "nrr"

/********************** DEFAULT PARAMETERS ************************/

const unsigned int   RIGID_ITERATIONS                      = 1200;
const          float RIGID_MAX_STEP                        = 2;
const          float RIGID_MIN_STEP                        = 0.01;

const unsigned int   RIGID_AFFINE_ITERATIONS               = 1200;
const          float RIGID_AFFINE_MAX_STEP                 = 2;
const          float RIGID_AFFINE_MIN_STEP                 = 0.2;

const unsigned int   AFFINE_ITERATIONS                     = 800;
const          float AFFINE_MAX_STEP                       = 2;
const          float AFFINE_MIN_STEP                       = 0.01;

const unsigned int   NRR_ITERATIONS_ONE                    = 800;
const unsigned int   NRR_ITERATIONS_TWO                    = 400;
const unsigned int   NRR_ITERATIONS_THREE                  = 200;
const          float NRR_MAX_STEP_ONE                      = 10;
const          float NRR_MAX_STEP_TWO                      = 8.6;
const          float NRR_MAX_STEP_THREE                    = 7;
const          float NRR_MIN_STEP                          = 0.00001;
const          float NRR_RELAXATION_FACTOR                 = 0.7;
const          float NRR_GRADIENT_MAGNITUDE_TOLERANCE      = 0.0000001;
const unsigned int   PYRAMID_RESOLUTION_ONE                = 2;
const unsigned int   PYRAMID_RESOLUTION_TWO                = 2;
const unsigned int   PYRAMID_RESOLUTION_THREE              = 2;
const unsigned int   GRID_POINT_SPACING_ONE                = 20;
const unsigned int   GRID_POINT_SPACING_TWO                = 10;
const unsigned int   GRID_POINT_SPACING_THREE              = 5;

const    unsigned int    SplineOrder = 3;
#define ImageDimension 3
#define SpaceDimension 3
/*******************************************************************/

#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
    typedef  CommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    CommandIterationUpdate() {};
public:
    typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;

    typedef const OptimizerType                       *  OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
        OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }
        if(optimizer->GetCurrentIteration()%50 == 0)
        {
            std::cout << optimizer->GetCurrentIteration() << "   ";
            std::cout << optimizer->GetValue() << "   ";
            std::cout << std::endl;
        }
    }
};

class RegistrationIterationUpdate : public itk::Command
{
public:
    typedef  RegistrationIterationUpdate   Self;
    typedef  itk::Command                  Superclass;
    typedef itk::SmartPointer<Self>        Pointer;
    itkNewMacro( Self );

    typedef float                                    InputPixelType;
#if ITK_VERSION_MAJOR < 4
    typedef itk::OrientedImage< InputPixelType, 3>   FixedImageType;
    typedef itk::OrientedImage< InputPixelType, 3 >  MovingImageType;
#else
    typedef itk::Image< InputPixelType, 3>   FixedImageType;
    typedef itk::Image< InputPixelType, 3 >  MovingImageType;
#endif
    typedef itk::MultiResolutionImageRegistrationMethod< FixedImageType, MovingImageType > RegistrationType;

    typedef   const RegistrationType   *    RegistrationPointer;
    typedef   RegistrationType         *    Registration2Pointer;

    typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
    typedef OptimizerType                           *    Optimizer2Pointer;

    typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType >    MetricType;
    typedef MetricType                                                                      *    Metric2Pointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        //Execute( (const itk::Object *)caller, event);

        Registration2Pointer registration = dynamic_cast< Registration2Pointer >( caller );
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }

        Optimizer2Pointer optimizer = dynamic_cast< Optimizer2Pointer >(registration->GetOptimizer());
        int iterations = optimizer->GetNumberOfIterations();
        int current = optimizer->GetCurrentIteration();
        if(current > 0)
        {
            std::cout << "Iteration event for registration entered: Prev Iterations: " << current << " of " << iterations << std::endl;
            //if(iterations == current)
            optimizer->SetNumberOfIterations( (iterations+1)/2 );
            std::cout << "    Next level use " << optimizer->GetNumberOfIterations() << " iterations" << std::endl;

            // Update Max-Min Step lengths
            //std::cout << "    Next level use " << optimizer->GetCurrentStepLength() << " Max Step Length" << std::endl;
            optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength()/2.0 );
            optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength()/2.0 );
            std::cout << "    Next level use " << optimizer->GetMinimumStepLength() << " Min Step Length" << std::endl;

            Metric2Pointer metric = dynamic_cast< Metric2Pointer >(registration->GetMetric());
            int samples = metric->GetNumberOfSpatialSamples();
            //std::cout << "    Metric samples previously used " << samples << std::endl;
            if(metric->GetUseAllPixels() == false)
            {
                metric->SetNumberOfSpatialSamples( static_cast<int>(4*samples) );
                //metric->SetFixedImageSamplesIntensityThreshold(-0.6);
                //metric->Initialize();
                std::cout << "    Will now use metric sampling: " << metric->GetNumberOfSpatialSamples() << std::endl;
            }
        }
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
        RegistrationPointer( dynamic_cast< RegistrationPointer >( object ) ); //ND
        if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
            return;
        }
    }

protected:
    RegistrationIterationUpdate() {};

};

Registration
::Registration()
{
    m_NearestNeighborInterpolator = NearestNeighborInterpolatorType::New();
    m_LinearInterpolator = LinearInterpolatorType::New();
    m_BsplineInterpolator = BsplineInterpolatorType::New();
    m_BsplineInterpolator->SetSplineOrder(3);

    m_NearestNeighborOrientatedInterpolator = NearestNeighborOrientatedInterpolatorType::New();
    m_LinearOrientatedInterpolator = LinearOrientatedInterpolatorType::New();
    m_BsplineOrientatedInterpolator = BsplineOrientatedInterpolatorType::New();
    m_BsplineOrientatedInterpolator->SetSplineOrder(3);

    m_ConsensusMode = VOTING;
    m_RegistrationMode = ALADIN_DIFFDEMONS;
    m_RegistrationInterpolation = INTERP_LINEAR;
    m_AtlasInterpolation = INTERP_BSPLINE;
    m_LabelledImageInterpolation = INTERP_NEAREST;

    m_UseConsensus = false;

    m_SaveFiles = 5;
}

Registration
::~Registration()
{

}

bool
Registration
::PropagateAffine(std::string dataFile,
                  std::string outFile,
                  std::string targetFile,
                  double * array, InterpolateMode mode,
                  itk::SpatialOrientation::ValidCoordinateOrientationFlags &/*flag*/,
                  itk::SpatialOrientationAdapter::DirectionType &dir)
{
    typedef itk::AffineTransform< double, 3 >  AffineTransformType;
    AffineTransformType::Pointer affineTransform = AffineTransformType::New();

    typedef AffineTransformType::ParametersType ParametersType;
    ParametersType params( 12 );

    for(int i = 0; i < 12; i++)
    {
        params[i] = array[i];
    }
    affineTransform->SetParameters(params);

    itk::ImageFileReader<InputImageType>::Pointer reader;
    reader = itk::ImageFileReader<InputImageType>::New();
    reader->SetFileName(dataFile);
    reader->Update();

    //InputImageType::SpacingType inputSpacing = reader->GetOutput()->GetSpacing();
    //InputImageType::SizeType inputSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

    itk::ImageFileReader<OutputImageType>::Pointer reader2;
    reader2 = itk::ImageFileReader<OutputImageType>::New();
    reader2->SetFileName(targetFile);
    reader2->Update();

    OutputImageType::SpacingType outputSpacing = reader2->GetOutput()->GetSpacing();
    OutputImageType::SizeType outputSize = reader2->GetOutput()->GetLargestPossibleRegion().GetSize();
    OutputImageType::PointType outputOrigin = reader2->GetOutput()->GetOrigin();

    // We want to resample based on affine transform and given spacing!
    typedef itk::ResampleImageFilter< InputImageType, OutputImageType >  ResampleFilterType;
    ResampleFilterType::Pointer resampleImageFilter = ResampleFilterType::New();
    if( mode == INTERP_NEAREST)
        resampleImageFilter->SetInterpolator( m_NearestNeighborInterpolator );
    else if(mode == INTERP_LINEAR)
        resampleImageFilter->SetInterpolator( m_LinearInterpolator );
    else if(mode == INTERP_BSPLINE)
    {
        m_BsplineInterpolator->SetSplineOrder(3);
        resampleImageFilter->SetInterpolator( m_BsplineInterpolator);
    }
    else
    {
        m_BsplineInterpolator->SetSplineOrder(5);
        resampleImageFilter->SetInterpolator( m_BsplineInterpolator);
    }
    resampleImageFilter->SetTransform( affineTransform );
    resampleImageFilter->SetDefaultPixelValue( 0 );

    // As Aladin transforms are in voxel space
    OutputImageType::SpacingType dummySpacing;
    dummySpacing[0] = dummySpacing[1] = dummySpacing[2] = 1.0;
    reader->GetOutput()->SetSpacing(dummySpacing);
    resampleImageFilter->SetOutputSpacing( dummySpacing );
    resampleImageFilter->SetSize( outputSize );
    OutputImageType::PointType dummyOrigin;
    dummyOrigin[0] = dummyOrigin[1] = dummyOrigin[2] = 0.0;
    resampleImageFilter->SetOutputOrigin( dummyOrigin );
    resampleImageFilter->SetInput(reader->GetOutput());
    resampleImageFilter->Update();

    OutputImageType::Pointer image = OutputImageType::New();
    OutputImageType::IndexType imageStart;
    imageStart[0] = 0;
    imageStart[1] = 0;
    imageStart[2] = 0;
    OutputImageType::RegionType imageRegion;
    imageRegion.SetSize(outputSize);
    imageRegion.SetIndex(imageStart);
    image->SetRegions(imageRegion);
    image->Allocate();
    //outputSpacing[0] = outputSpacing[0];
    //outputSpacing[1] = outputSpacing[1];
    //outputSpacing[2] = outputSpacing[2];
    image->SetSpacing(outputSpacing);
    image->SetOrigin(outputOrigin);

    typedef itk::ImageRegionIteratorWithIndex<OutputImageType> ImageRegionIteratorWithIndexType;
    ImageRegionIteratorWithIndexType it( resampleImageFilter->GetOutput(), resampleImageFilter->GetOutput()->GetRequestedRegion() );
    it.GoToBegin();

    OutputImageType::IndexType index;
    while ( ! it.IsAtEnd() )
    {
        index = it.GetIndex();
        image->SetPixel(index, it.Value());
        ++it;
    }
    image->SetDirection(dir);

    // Save resampled image
    itk::ImageFileWriter<OutputImageType>::Pointer writer;
    writer = itk::ImageFileWriter<OutputImageType>::New();
    writer->SetInput(image);
    writer->SetFileName(outFile);
    writer->Write();

    return true;
}

// Simple function to convert a string into an array of char
// which can be passed as an input for the nifty-reg executables
// which have been converted in libraries but still take (int argc, char *argv[]) as input
unsigned int
Registration
::StringToCharArray(std::string & cmdline, char ** argument )
{
    unsigned int count=0;

    std::string word;

    // transform the string into a stream
    std::stringstream stream(cmdline);

    // split the string
    while( getline(stream, word, ' ') )
    {
        char *tmp_word;
        // allocate a new char array
        tmp_word = new char[word.length() + 1];
        // copy the substring into the char array
        std::strcpy(tmp_word,word.c_str());
        // assign the pointer of the char array to the argument array
        argument[count++] = tmp_word;
    }
    // return argument count
    return count;
}

bool
Registration
::RegisterMilxNonRigid(std::string dataFile, std::string atlasFile,
                       std::string outputName, bool /*Invert*/)
{

    if (dataFile.length() == 0 || atlasFile.length() == 0) {
        return false;
    }

    /*
    milx::Orientation::ValidCoordinateOrientationFlagsType inputFlag;
    milx::Orientation::DirectionType inputDirection;
    milx::Orientation::ValidCoordinateOrientationFlagsType atlasFlag;
    milx::Orientation::DirectionType atlasDirection;
    milx::Orientation *imageOrientation = new milx::Orientation();
    bool result = imageOrientation->GetSpatialOrientation(dataFile.c_str(), inputFlag, inputDirection);
    if (result == false) {
        delete imageOrientation;
        return false;
    }
    result = imageOrientation->GetSpatialOrientation(atlasFile.c_str(), atlasFlag, atlasDirection);
    if (result == false) {
        delete imageOrientation;
        return false;
    }
    // Need to apply orientation to each output file.
    delete imageOrientation;
    */

    itk::ImageFileReader<InputImageType>::Pointer reader = itk::ImageFileReader<InputImageType>::New();
    reader->SetFileName(dataFile);
    try
    {
        reader->Update();
    }
    catch (itk::ExceptionObject & ex)
    {
        std::cerr << "Failed reading data file" << std::endl;
        std::cerr << ex.GetDescription() << std::endl;
        return false;
    }

    itk::ImageFileReader<InputImageType>::Pointer atlasReader = itk::ImageFileReader<InputImageType>::New();
    atlasReader->SetFileName(atlasFile);
    try
    {
        atlasReader->Update();
    }
    catch (itk::ExceptionObject & ex)
    {
        std::cerr << "Failed reading atlas file" << std::endl;
        std::cerr << ex.GetDescription() << std::endl;
        return false;
    }

    InputImageType::Pointer atlas = atlasReader->GetOutput();
    InputImageType::Pointer mriData = reader->GetOutput();

    typedef itk::HistogramMatchingImageFilter <InputImageType, InputImageType> MatchingFilterType;
    MatchingFilterType::Pointer matcher = MatchingFilterType::New();
    matcher->SetInput(atlas);
    matcher->SetReferenceImage(mriData);
    matcher->SetNumberOfHistogramLevels(1024);
    matcher->SetNumberOfMatchPoints(7);
    matcher->ThresholdAtMeanIntensityOn();

    typedef itk::PDEDeformableRegistrationFilter < InputImageType, InputImageType, DemonsDeformationFieldType> BaseRegistrationFilterType;
    BaseRegistrationFilterType::Pointer filter;

    // s <- s o exp(u) (Diffeomorphic demons)
    typedef itk::DiffeomorphicDemonsRegistrationFilter < InputImageType, InputImageType, DemonsDeformationFieldType> ActualRegistrationFilterType;
    typedef ActualRegistrationFilterType::GradientType GradientType;

    ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();

    actualfilter->SetMaximumUpdateStepLength(2.0);
    actualfilter->SetUseGradientType(static_cast<GradientType>(0));
    filter = actualfilter;
#if ITK_VERSION_MAJOR < 4
    filter->SmoothDeformationFieldOn();
#else
    filter->SetSmoothDisplacementField(true);
#endif
    filter->SetStandardDeviations(3.0);


    typedef itk::MultiResolutionPDEDeformableRegistration< InputImageType, InputImageType, DemonsDeformationFieldType, float >   MultiResRegistrationFilterType;
    MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

    int numLevels = 3;
    std::vector<unsigned int> numIterations;
    numIterations = std::vector<unsigned int>(numLevels, 10u);

    multires->SetRegistrationFilter(filter);
    multires->SetNumberOfLevels(numLevels);
    multires->SetNumberOfIterations(&numIterations[0]);

    multires->SetFixedImage(mriData);
    multires->SetMovingImage(matcher->GetOutput());

    // Compute the deformation field
    try
    {
        multires->UpdateLargestPossibleRegion();
    }
    catch (itk::ExceptionObject& err)
    {
        std::cout << "Unexpected error." << std::endl;
        std::cout << err << std::endl;
        return false;
    }

    // The outputs
    DemonsDeformationFieldType::Pointer defField = multires->GetOutput();

    // We want to save deformation field
    // and propogated aal

    itk::ImageFileWriter<DemonsDeformationFieldType>::Pointer defWriter = itk::ImageFileWriter<DemonsDeformationFieldType>::New();
    defWriter->SetInput(defField);
    std::string filename1 = outputName;
    defWriter->SetFileName(filename1);
    try
    {
        defWriter->Write();
    }
    catch (itk::ExceptionObject & ex)
    {
        std::cerr << "Failed writing file" << std::endl;
        std::cerr << ex.GetDescription() << std::endl;
        return false;
    }

    return true;
}

void
Registration
::LoadTRSF(std::string filename1, std::string filename2, double * array, bool invert)
{
    std::cout << "Loading file with name " << filename1 << std::endl;
    std::cout << "Loading file with name " << filename2 << std::endl;
    char * buffer1 = new char[512];
    char * buffer2 = new char[512];
    gzFile fin1 = ::gzopen( filename1.c_str(), "rb" );
    gzFile fin2 = ::gzopen( filename2.c_str(), "rb" );
    if((fin1 == NULL) || (fin2 == NULL))
    {
        std::cerr << "Cannot read " << filename1 << std::endl;
        delete [] buffer1;
        delete [] buffer2;
    }
    buffer1 = ::gzgets (fin1, buffer1, 512);
    buffer2 = ::gzgets (fin2, buffer2, 512);
    if((buffer1[0] != '(') || (buffer2[0] != '('))
    {
        std::cerr << "File is not a transform file " << buffer1[0] << std::endl;
    }
    buffer1 = ::gzgets (fin1, buffer1, 512);
    buffer2 = ::gzgets (fin2, buffer2, 512);
    if(((buffer1[0] != 0) && (buffer1[1] != '8')) || (((buffer2[0] != 0) && (buffer2[1] != '8'))))
    {
        // Failed Magic number test ;)
        std::cerr << buffer1 << std::endl;
        std::cerr << "File is not a transform file" << buffer1[1] << std::endl;
    }
    char * str1 = new char [128];
    char * str2 = new char [128];
    char * str3 = new char [128];
    char * str4 = new char [128];

    vnl_matrix<double> matrix1(4,4);
    vnl_matrix<double> matrix2(4,4);
    for (unsigned long i = 0; i < 4; i++)
    {
        buffer1 = ::gzgets (fin1, buffer1, 512);
        sscanf(buffer1,"%s %s %s %s", str1, str2, str3, str4);
        matrix1(i, 0) = atof(str1);
        matrix1(i, 1) = atof(str2);
        matrix1(i, 2) = atof(str3);
        matrix1(i, 3) = atof(str4);

        buffer2 = ::gzgets (fin2, buffer2, 512);
        sscanf(buffer2,"%s %s %s %s", str1, str2, str3, str4);
        matrix2(i, 0) = atof(str1);
        matrix2(i, 1) = atof(str2);
        matrix2(i, 2) = atof(str3);
        matrix2(i, 3) = atof(str4);
        //std::cout << "(" << str1 << "," << str2 << "," << str3 << "," << str4 << ")" << std::endl;
    }

    vnl_matrix<double> matrixComposite = matrix2*matrix1;
    //std::cout << "Transform Matrix 1 is " << std::endl;
    //std::cout << matrixComposite << std::endl;

    if(invert == true)
    {
        matrixComposite = vnl_matrix_inverse<double>(matrixComposite);
        //std::cout << "Inverse Matrix is " << std::endl;
        //std::cout << matrixComposite << std::endl;
    }

    //matrix.set_identity();
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
            array[3*i+j] = matrixComposite(i,j);
        array[9+i] = matrixComposite(i, 3);
    }

    // Get last line
    buffer1 = ::gzgets (fin1, buffer1, 512);
    if(buffer1[0] != ')')
    {
        std::cerr << buffer1 << std::endl;
        std::cerr << "File is not a valid trsf transform file ?" << std::endl;
    }
    delete [] buffer1;
    delete [] buffer2;
    delete [] str1;
    delete [] str2;
    delete [] str3;
    delete [] str4;
    gzclose(fin1);
    gzclose(fin2);
}

void
Registration
::LoadTRSF(std::string filename, double * array, bool invert)
{
    std::cout << "Loading file with name " << filename << std::endl;
    char * buffer = new char[512];
    gzFile fin = ::gzopen( filename.c_str(), "rb" );
    if(fin == NULL)
    {
        std::cerr << "Cannot read " << filename << std::endl;
        delete [] buffer;
    }
    buffer = ::gzgets (fin, buffer, 512);
    if(buffer[0] != '(')
    {
        std::cerr << "File is not a transform file " << buffer[0] << std::endl;
    }
    buffer = ::gzgets (fin, buffer, 512);
    if((buffer[0] != 0) && (buffer[1] != '8'))
    {
        // Failed Magic number test ;)
        std::cerr << buffer << std::endl;
        std::cerr << "File is not a transform file" << buffer[1] << std::endl;
    }
    char * str1 = new char [128];
    char * str2 = new char [128];
    char * str3 = new char [128];
    char * str4 = new char [128];

    vnl_matrix<double> matrix(4,4);
    for (unsigned long i = 0; i < 4; i++)
    {
        buffer = ::gzgets (fin, buffer, 512);
        sscanf(buffer,"%s %s %s %s", str1, str2, str3, str4);
        matrix(i, 0) = atof(str1);
        matrix(i, 1) = atof(str2);
        matrix(i, 2) = atof(str3);
        matrix(i, 3) = atof(str4);
        //std::cout << "(" << str1 << "," << str2 << "," << str3 << "," << str4 << ")" << std::endl;
    }
    //std::cout << "Transform Matrix is " << std::endl;
    //std::cout << matrix << std::endl;

    if(invert == true)
    {
        matrix = vnl_matrix_inverse<double>(matrix);
        //std::cout << "Inverse Matrix is " << std::endl;
        //std::cout << matrix << std::endl;
    }

    //matrix.set_identity();
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
            array[3*i+j] = matrix(i,j);
        array[9+i] = matrix(i, 3);
    }

    // Get last line
    buffer = ::gzgets (fin, buffer, 512);
    if(buffer[0] != ')')
    {
        std::cerr << buffer << std::endl;
        std::cerr << "File is not a valid trsf transform file ?" << std::endl;
    }
    delete [] buffer;
    delete [] str1;
    delete [] str2;
    delete [] str3;
    delete [] str4;
    gzclose(fin);
}

bool
Registration
::RegisterRigidITK(std::string dataFile, std::string atlasFile,
                   std::string outputName, bool /*Invert*/,
                   std::string dataMaskFile, std::string atlasMaskFile)
{
    bool writeFiles = false;
    char outputFileName[512];
    typedef double CoordinateRepType;

    typedef itk::VersorRigid3DTransform< CoordinateRepType > RigidTransformType;
    typedef itk::CenteredTransformInitializer< RigidTransformType, FixedImageType, MovingImageType > TransformInitializerType;
    typedef itk::LinearInterpolateImageFunction< MovingImageType, double> InterpolatorType;
    typedef itk::MultiResolutionImageRegistrationMethod< FixedImageType, MovingImageType > RegistrationType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType >    MetricType;
    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >  ResampleFilterType;
    typedef itk::RecursiveMultiResolutionPyramidImageFilter< FixedImageType, FixedImageType > PyramidType;

    OptimizerType::Pointer      optimizer              = OptimizerType::New();
    InterpolatorType::Pointer   interpolator           = InterpolatorType::New();
    RegistrationType::Pointer   registration           = RegistrationType::New();
    MetricType::Pointer         metric                  = MetricType::New();
    TransformInitializerType::Pointer initializer           = TransformInitializerType::New();

    // Enable the use of Cubic spline interpolation in registration scheme

    registration->SetMetric( metric );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );

    RegistrationIterationUpdate::Pointer observer1 = RegistrationIterationUpdate::New();
    registration->AddObserver( itk::IterationEvent(), observer1 );

    // Auxiliary identity transform.
    typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
    IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

    // Read in the files and set them to the registration object
    MovingImageType::Pointer movingImage = this->readFile3D( atlasFile );
    FixedImageType::Pointer fixedImage   = this->readFile3D( dataFile );

    // This is another copy of the moving image to be used for transforming (re-orientating) during optimisation
    registration->SetFixedImage(  fixedImage   );
    registration->SetMovingImage(   movingImage  );

    // Add a time and memory probes collector for profiling the computation time of every stage.
    //itkProbesCreate();

    /*********************************************************************************
     *                                                                               *
     *                          Perform Rigid Registration                           *
     *                                                                               *
     *********************************************************************************/

    RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
    initializer->SetTransform( rigidTransform );
    initializer->SetFixedImage(  fixedImage );
    initializer->SetMovingImage( movingImage );
    //initializer->MomentsOn();
    initializer->GeometryOn();
    initializer->InitializeTransform();
    // Display initial transform
    std::cout << "Initial Rigid Parameters " << rigidTransform->GetParameters() << std::endl;


    FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
    registration->SetFixedImageRegion( fixedRegion );
    registration->SetInitialTransformParameters( rigidTransform->GetParameters() );
    registration->SetTransform( rigidTransform );

    //itk::Array2D<unsigned int> schedule;
    RegistrationType::ScheduleType fixedSchedule(3,3);
    fixedSchedule[0][0] = 8;
    fixedSchedule[0][1] = 8;
    fixedSchedule[0][2] = 8;
    fixedSchedule[1][0] = 4;
    fixedSchedule[1][1] = 4;
    fixedSchedule[1][2] = 4;
    fixedSchedule[2][0] = 2;
    fixedSchedule[2][1] = 2;
    fixedSchedule[2][2] = 2;

    RegistrationType::ScheduleType movingSchedule(3,3);
    movingSchedule[0][0] = 8;
    movingSchedule[0][1] = 8;
    movingSchedule[0][2] = 8;
    movingSchedule[1][0] = 4;
    movingSchedule[1][1] = 4;
    movingSchedule[1][2] = 4;
    movingSchedule[2][0] = 2;
    movingSchedule[2][1] = 2;
    movingSchedule[2][2] = 2;
    registration->SetSchedules(fixedSchedule, movingSchedule);
    //registration->SetNumberOfLevels(3);

    // Add mask
    typedef itk::ImageMaskSpatialObject< SpaceDimension >   MaskType;
    typedef itk::Image< unsigned char, SpaceDimension >     ImageMaskType;
    typedef itk::ImageFileReader< ImageMaskType >           MaskReaderType;

    typedef itk::LabelStatisticsImageFilter<FixedImageType, ImageMaskType> LabelStatisticsImageFilterType;
    LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();

    int numberOfVoxelsInMask = 0;
    bool useFixedMask = false;

    if ( dataMaskFile.length() != 0) {
        useFixedMask = true;
    }


    MaskReaderType::Pointer fixedImageMaskReader  = MaskReaderType::New();
    MaskType::Pointer spatialObjectMask = MaskType::New();
    if(useFixedMask)
    {
        std::cout << "RegisterRigidITK: Using Fixed mask " << dataMaskFile <<  std::endl;
        fixedImageMaskReader->SetFileName( dataMaskFile );
        try {
            fixedImageMaskReader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        typedef itk::BinaryThresholdImageFilter<ImageMaskType,ImageMaskType> BinaryThresholdImageFilterType;
        BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
        threshold->SetInput(fixedImageMaskReader->GetOutput());
        threshold->SetLowerThreshold(1);
        threshold->SetUpperThreshold(255);
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(255);
        threshold->Update();

        spatialObjectMask->SetImage( threshold->GetOutput() );
        useFixedMask = true;

        labelStatisticsImageFilter->SetInput(fixedImage);
        labelStatisticsImageFilter->SetLabelInput(threshold->GetOutput());
        labelStatisticsImageFilter->Update();
        //for(int i = 50; i < 300; i++)
        numberOfVoxelsInMask = labelStatisticsImageFilter->GetCount(255);
        std::cout << "We have " << numberOfVoxelsInMask << " in fixed mask region" << std::endl;

    }
    bool useMovingMask = false;
    if ( atlasMaskFile.length() != 0)
        useMovingMask = true;

    MaskReaderType::Pointer movingImageMaskReader  = MaskReaderType::New();
    MaskType::Pointer movingSpatialObjectMask = MaskType::New();
    if(useMovingMask)
    {
        std::cout << "RegisterRigidITK: Using Moving mask: " << atlasMaskFile << std::endl;
        movingImageMaskReader->SetFileName( atlasMaskFile );
        try {
            movingImageMaskReader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        typedef itk::BinaryThresholdImageFilter<ImageMaskType,ImageMaskType> BinaryThresholdImageFilterType;
        BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
        threshold->SetInput(movingImageMaskReader->GetOutput());
        threshold->SetLowerThreshold(1);
        threshold->SetUpperThreshold(255);
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(255);
        threshold->Update();

        movingSpatialObjectMask->SetImage( threshold->GetOutput() );
        useMovingMask = true;

        labelStatisticsImageFilter->SetInput(movingImage);
        labelStatisticsImageFilter->SetLabelInput(threshold->GetOutput());
        labelStatisticsImageFilter->Update();
        //for(int i = 50; i < 300; i++)
        numberOfVoxelsInMask = labelStatisticsImageFilter->GetCount(255);
        std::cout << "We have " << numberOfVoxelsInMask << " in moving mask region" << std::endl;
    }

    const unsigned int numberOfPixels = fixedRegion.GetNumberOfPixels();
    metric->SetNumberOfHistogramBins( 128 );
    metric->SetUseExplicitPDFDerivatives( true );
    metric->SetUseCachingOfBSplineWeights(false);
    metric->ReinitializeSeed( 76926294 );

    // Define optimizer normalization to compensate for different dynamic range of rotations and translations
    typedef OptimizerType::ScalesType       OptimizerScalesType;
    OptimizerScalesType optimizerScales( rigidTransform->GetNumberOfParameters() );
    const double translationScale = 1.0 / 4000.0;
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = translationScale;
    optimizerScales[4] = translationScale;
    optimizerScales[5] = translationScale;
    optimizer->SetScales( optimizerScales );

    optimizer->SetMaximumStepLength( RIGID_MAX_STEP  );
    optimizer->SetMinimumStepLength( RIGID_MIN_STEP );
    optimizer->SetNumberOfIterations( RIGID_ITERATIONS );

    // Create the Command observer and register it with the optimizer.
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );

    std::cout << "Number of Pixels " << numberOfPixels << std::endl;

    if(useFixedMask)
    {
        metric->SetFixedImageMask( spatialObjectMask );
        if(numberOfVoxelsInMask >  numberOfPixels/4.0)
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4*16)) );
        else
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4)) );
    }
    else
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / (256)) );

    std::cout << "Initialize metric with " << metric->GetNumberOfSpatialSamples() << " samples " << std::endl;

    std::cout << "Starting Rigid Registration " << std::endl;

    try {
        //itkProbesStart( "Rigid Registration" );
#if ITK_VERSION_MAJOR < 4
        registration->StartRegistration();
#else
        registration->Update();
#endif    //itkProbesStop( "Rigid Registration" );
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Rigid Registration completed" << std::endl;
    std::cout << std::endl;

    rigidTransform->SetParameters( registration->GetLastTransformParameters() );

    /*********************************************************************************
     *                                                                               *
     *                    Output the Rigid Parameters                               *
     *                                                                               *
     *********************************************************************************/

    //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"rigid.txt");
    sprintf(outputFileName, "%s%s", outputName.c_str(), "_rigidtransform.txt");

    itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
    transformWriter->AddTransform(rigidTransform);
    transformWriter->SetFileName(outputFileName);
    try {
        transformWriter->Update();
    } catch (itk::ExceptionObject & excp) {
        std::cerr << excp << std::endl;
        return false;
    }

    /*********************************************************************************
     *                                                                               *
     *                    Output the Rigid Transformed Image                        *
     *                                                                               *
     *********************************************************************************/

    if (writeFiles == true) {
        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"rigid_orig.nii.gz");
        sprintf(outputFileName, "%s", outputName.c_str());
        ResampleFilterType::Pointer resampleAffine = ResampleFilterType::New();
        resampleAffine->SetTransform( rigidTransform );
        resampleAffine->SetInput( movingImage );
        resampleAffine->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleAffine->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleAffine->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleAffine->SetOutputDirection( fixedImage->GetDirection() );
        resampleAffine->SetDefaultPixelValue(0.0);
        resampleAffine->Update();

        std::cout << std::endl << "Writing  Rigid resampled moving image..." << std::flush;
        writeFile3D( resampleAffine->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl << std::endl;
    }

    return true;
}

bool
Registration
::RegisterAffineITK(std::string dataFile, std::string atlasFile,
                    std::string outputName, bool /*Invert*/, std::string dataMaskFile,
                    std::string atlasMaskFile)
{
    bool writeFiles = true;
    char outputFileName[512];
    typedef double CoordinateRepType;

    typedef itk::VersorRigid3DTransform< CoordinateRepType > RigidTransformType;
    typedef itk::AffineTransform< CoordinateRepType, SpaceDimension > AffineTransformType;
    //typedef itk::BSplineDeformableTransform< CoordinateRepType, SpaceDimension, SplineOrder > DeformableTransformType;
    typedef itk::CenteredTransformInitializer< RigidTransformType, FixedImageType, MovingImageType > TransformInitializerType;
    //typedef itk::CenteredVersorTransformInitializer< MovingImageType, MovingImageType > TransformInitializerType;
    typedef itk::LinearInterpolateImageFunction< MovingImageType, double> InterpolatorType;
    //typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType > ImageRegistrationType;
    typedef itk::MultiResolutionImageRegistrationMethod< FixedImageType, MovingImageType > RegistrationType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    //typedef itk::VersorRigid3DTransformOptimizer VersorOptimizerType;
    typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType >    MetricType;
    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >  ResampleFilterType;
    typedef itk::RecursiveMultiResolutionPyramidImageFilter< FixedImageType, FixedImageType > PyramidType;

    OptimizerType::Pointer      optimizer              = OptimizerType::New();
    //VersorOptimizerType::Pointer      VersorOptimizer              = VersorOptimizerType::New();
    InterpolatorType::Pointer   interpolator           = InterpolatorType::New();
    RegistrationType::Pointer   registration           = RegistrationType::New();
    MetricType::Pointer         metric                  = MetricType::New();
    TransformInitializerType::Pointer initializer           = TransformInitializerType::New();

    registration->SetMetric( metric );
    //registration->SetOptimizer(     VersorOptimizer     );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );

    RegistrationIterationUpdate::Pointer observer1 = RegistrationIterationUpdate::New();
    registration->AddObserver( itk::IterationEvent(), observer1 );

    // Auxiliary identity transform.
    //typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
    //IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

    // Read in the files and set them to the registration object
    MovingImageType::Pointer movingImage = this->readFile3D( atlasFile );
    FixedImageType::Pointer fixedImage   = this->readFile3D( dataFile );

    // This is another copy of the moving image to be used for transforming (re-orientating) during optimisation
    registration->SetFixedImage(  fixedImage   );
    registration->SetMovingImage(   movingImage  );

    FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
    FixedImageType::SizeType fixedImageSize       = fixedRegion.GetSize();
    FixedImageType::SpacingType fixedImageSpacing = fixedImage->GetSpacing();
    std::cout << "ImageSize: " << fixedImageSize[0] << ", " << fixedImageSize[1] << ", " << fixedImageSize[2] << std::endl;
    std::cout << "ImageSpacing: " << fixedImageSpacing[0] << ", " <<  fixedImageSpacing[1] << ", " <<  fixedImageSpacing[2] << std::endl;
    std::cout << "ImageOrigin: " << fixedImage->GetOrigin() << std::endl;

    // Add a time and memory probes collector for profiling the computation time of every stage.
    //itkProbesCreate();

    /*********************************************************************************
     *                                                                               *
     *                          Perform Rigid Registration                           *
     *                                                                               *
     *********************************************************************************/

    RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
    initializer->SetTransform( rigidTransform );
    initializer->SetFixedImage(  fixedImage );
    initializer->SetMovingImage( movingImage );
    //initializer->MomentsOn();
    initializer->GeometryOn();
    initializer->InitializeTransform();

    registration->SetFixedImageRegion( fixedRegion );
    registration->SetInitialTransformParameters( rigidTransform->GetParameters() );
    registration->SetTransform( rigidTransform );

    RegistrationType::ScheduleType fixedSchedule(3,3);
    fixedSchedule[0][0] = 8;
    fixedSchedule[0][1] = 8;
    fixedSchedule[0][2] = 8;
    fixedSchedule[1][0] = 4;
    fixedSchedule[1][1] = 4;
    fixedSchedule[1][2] = 4;
    fixedSchedule[2][0] = 2;
    fixedSchedule[2][1] = 2;
    fixedSchedule[2][2] = 2;

    RegistrationType::ScheduleType movingSchedule(3,3);
    movingSchedule[0][0] = 8;
    movingSchedule[0][1] = 8;
    movingSchedule[0][2] = 8;
    movingSchedule[1][0] = 4;
    movingSchedule[1][1] = 4;
    movingSchedule[1][2] = 4;
    movingSchedule[2][0] = 2;
    movingSchedule[2][1] = 2;
    movingSchedule[2][2] = 2;
    registration->SetSchedules(fixedSchedule, movingSchedule);
    //registration->SetNumberOfLevels(3);

    // Add mask
    typedef itk::ImageMaskSpatialObject< SpaceDimension >         MaskType;
#if ITK_VERSION_MAJOR < 4
    typedef itk::OrientedImage< unsigned char, SpaceDimension >   ImageMaskType;
#else
    typedef itk::Image< unsigned char, SpaceDimension >   ImageMaskType;
#endif
    typedef itk::ImageFileReader< FixedImageType >                MaskReaderType;

    int numberOfVoxelsInMask = 0;
    bool useFixedMask = false;

    if ( dataMaskFile.length() != 0) {
        useFixedMask = true;
    }

    MaskReaderType::Pointer fixedImageMaskReader  = MaskReaderType::New();
    MaskType::Pointer spatialObjectMask = MaskType::New();

    typedef itk::LabelStatisticsImageFilter<FixedImageType, ImageMaskType> LabelStatisticsImageFilterType;
    LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();

    if(useFixedMask)
    {
        std::cout << "RegisterAffineITK: Using Fixed mask " << dataMaskFile <<  std::endl;
        fixedImageMaskReader->SetFileName( dataMaskFile );
        try {
            fixedImageMaskReader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return false;
        }
        typedef itk::BinaryThresholdImageFilter<FixedImageType,ImageMaskType> BinaryThresholdImageFilterType;
        BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
        threshold->SetInput(fixedImageMaskReader->GetOutput());
        threshold->SetLowerThreshold(0.1);
        threshold->SetUpperThreshold(255);
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(255);
        threshold->Update();

        std::cout << "Direction " << threshold->GetOutput()->GetDirection() << std::endl;
        std::cout << "Origin " << threshold->GetOutput()->GetOrigin() << std::endl;

        spatialObjectMask->SetImage( threshold->GetOutput() );
        useFixedMask = true;


        labelStatisticsImageFilter->SetInput(fixedImage);
        labelStatisticsImageFilter->SetLabelInput(threshold->GetOutput());
        labelStatisticsImageFilter->Update();
        //for(int i = 50; i < 300; i++)
        numberOfVoxelsInMask = labelStatisticsImageFilter->GetCount(255);
        std::cout << "We have " << numberOfVoxelsInMask << " in fixed mask region" << std::endl;
    }

    bool useMovingMask = false;
    if ( atlasMaskFile.length() != 0) {
        useMovingMask = true;
    }

    MaskReaderType::Pointer movingImageMaskReader  = MaskReaderType::New();
    MaskType::Pointer movingSpatialObjectMask = MaskType::New();
    if(useMovingMask)
    {
        std::cout << "RegisterAffineITK: Using Moving mask: " << atlasMaskFile << std::endl;
        movingImageMaskReader->SetFileName( atlasMaskFile );
        try {
            movingImageMaskReader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }

        FixedImageType::Pointer image = FixedImageType::New();
        FixedImageType::IndexType imageStart;
        imageStart[0] = 0;
        imageStart[1] = 0;
        imageStart[2] = 0;
        FixedImageType::RegionType imageRegion;
        imageRegion.SetSize(movingImageMaskReader->GetOutput()->GetLargestPossibleRegion().GetSize());
        imageRegion.SetIndex(imageStart);
        image->SetRegions(imageRegion);
        image->Allocate();

        image->SetSpacing(movingImageMaskReader->GetOutput()->GetSpacing());
        image->SetOrigin(movingImageMaskReader->GetOutput()->GetOrigin());
        //image->SetDirection(reader->GetOutput()->GetDirection());

        // Apply Direction Cosine to mask as ITK doesn't support this properly.
        // Need to apply Direction and resample
        //typedef itk::IdentityTransform< InputPixelType, 3 >  IdentityTransformType;
        //IdentityTransformType::Pointer identityInterpolator = IdentityTransformType::New();
        typedef itk::ResampleImageFilter< FixedImageType, FixedImageType >  ResampleFilterType2;
        ResampleFilterType2::Pointer resampleAffine = ResampleFilterType2::New();
        //resampleAffine->SetInterpolator( identityInterpolator );
        //resampleAffine->SetTransform( affineTransform );
        resampleAffine->SetInput( movingImageMaskReader->GetOutput() );
        resampleAffine->SetReferenceImage(image);
        resampleAffine->SetUseReferenceImage(true);
        // By default 1,1,1 ?
        //resampleAffine->SetOutputDirection( fixedImage->GetDirection() );
        resampleAffine->SetDefaultPixelValue(0.0);
        try {
            resampleAffine->Update();
        } catch (itk::ExceptionObject & err) {
            std::cerr << "ExceptionObject caught:" << std::endl;
            std::cerr << err << std::endl;
            exit(1);
        }

        typedef itk::BinaryThresholdImageFilter<FixedImageType,ImageMaskType> BinaryThresholdImageFilterType;
        BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
        threshold->SetInput(resampleAffine->GetOutput());
        threshold->SetLowerThreshold(0.1);
        threshold->SetUpperThreshold(255);
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(255);
        threshold->Update();

        std::cout << "Direction " << threshold->GetOutput()->GetDirection() << std::endl;
        std::cout << "Origin " << threshold->GetOutput()->GetOrigin() << std::endl;

        //movingSpatialObjectMask->DebugOn();
        movingSpatialObjectMask->SetImage( threshold->GetOutput() );
        std::cout << movingSpatialObjectMask->GetIndexToWorldTransform()->GetMatrix() << std::endl;
        std::cout << movingSpatialObjectMask->GetIndexToObjectTransform()->GetMatrix()  << std::endl;
        useMovingMask = true;

        // Save Spatial Object
        typedef itk::SpatialObjectToImageFilter< MaskType, ImageMaskType > SpatialObjectToImageFilterType;
        SpatialObjectToImageFilterType::Pointer spat = SpatialObjectToImageFilterType::New();
        spat->SetInput(movingSpatialObjectMask);
        //spat->SetDirection(threshold->GetOutput()->GetDirection());
        spat->SetOrigin(threshold->GetOutput()->GetOrigin());
        spat->SetSpacing(threshold->GetOutput()->GetSpacing());
        spat->SetSize(threshold->GetOutput()->GetLargestPossibleRegion().GetSize());

        spat->Update();

        typedef itk::ImageFileWriter< ImageMaskType >  WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput(spat->GetOutput());
        writer->SetFileName("testMask.nii.gz");
        try {
            writer->Update();
        } catch (itk::ExceptionObject & err) {
            std::cerr << "ExceptionObject caught:" << std::endl;
            std::cerr << err << std::endl;
            exit(1);
        }

    }

    const unsigned int numberOfPixels = fixedRegion.GetNumberOfPixels();
    metric->SetNumberOfHistogramBins( 128 );
    metric->SetUseExplicitPDFDerivatives( true );
    metric->SetUseCachingOfBSplineWeights(false);
    metric->ReinitializeSeed( 76926294 );

    // Define optimizer normalization to compensate for different dynamic range of rotations and translations
    typedef OptimizerType::ScalesType       OptimizerScalesType;
    OptimizerScalesType optimizerScales( rigidTransform->GetNumberOfParameters() );
    const double translationScale = 1.0 / 4000.0;
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = translationScale;
    optimizerScales[4] = translationScale;
    optimizerScales[5] = translationScale;
    optimizer->SetScales( optimizerScales );

    optimizer->SetMaximumStepLength( RIGID_AFFINE_MAX_STEP  );
    optimizer->SetMinimumStepLength( RIGID_AFFINE_MIN_STEP );
    optimizer->SetNumberOfIterations( RIGID_AFFINE_ITERATIONS );

    // Create the Command observer and register it with the optimizer.
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );

    //CommandIterationUpdate::Pointer VersorObserver = CommandIterationUpdate::New();
    //VersorOptimizer->AddObserver( itk::IterationEvent(), VersorObserver );

    if(useFixedMask)
    {
        metric->SetFixedImageMask( spatialObjectMask );
        if(numberOfVoxelsInMask >  numberOfPixels/4.0)
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4*16)) );
        else
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4)) );
    }
    else if(useMovingMask)
    {
        metric->DebugOn();
        metric->SetMovingImageMask( movingSpatialObjectMask );
        if(numberOfVoxelsInMask >  numberOfPixels/4.0)
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4*16)) );
        else
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4)) );
    }
    else
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / (256)) );

    std::cout << "Initialize metric with " << metric->GetNumberOfSpatialSamples() << std::endl;

    std::cout << "Starting Rigid Registration " << std::endl;

    try {
        //itkProbesStart( "Rigid Registration" );
#if ITK_VERSION_MAJOR < 4
        registration->StartRegistration();
#else
        registration->Update();
#endif
        //itkProbesStop( "Rigid Registration" );
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return false;
    }

    std::cout << "Rigid Registration completed" << std::endl;
    std::cout << std::endl;

    rigidTransform->SetParameters( registration->GetLastTransformParameters() );

    /*********************************************************************************
     *                                                                               *
     *                    Output the Rigid Parameters                               *
     *                                                                               *
     *********************************************************************************/

    /*
    if (writeFiles == true) {
        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"_rigid.txt");
        sprintf(outputFileName, "%s%s", outputName.c_str(), "_rigidtransform.txt");

        itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
        transformWriter->AddTransform(rigidTransform);
        transformWriter->SetFileName(outputFileName);
        try {
            transformWriter->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return false;
        }
    }
    */
    /*********************************************************************************
     *                                                                               *
     *                    Output the Rigid Transformed Image                        *
     *                                                                               *
     *********************************************************************************/

    writeFiles = false;
    if (writeFiles == true) {
        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"rigid.nii.gz");
        sprintf(outputFileName, "%s", outputName.c_str());
        ResampleFilterType::Pointer resampleAffine = ResampleFilterType::New();
        resampleAffine->SetTransform( rigidTransform );
        resampleAffine->SetInput( movingImage );
        resampleAffine->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleAffine->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleAffine->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleAffine->SetOutputDirection( fixedImage->GetDirection() );
        resampleAffine->SetDefaultPixelValue(0.0);
        resampleAffine->Update();

        std::cout << std::endl << "Writing  Rigid resampled moving image..." << std::flush;
        writeFile3D( resampleAffine->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl << std::endl;
    }
    writeFiles = true;


    /*********************************************************************************
     *                                                                               *
     *                         Perform Affine Registration                           *
     *                                                                               *
     *********************************************************************************/

    registration->SetOptimizer(     optimizer     );

    // Setup the affine transform based on the result of the rigid
    AffineTransformType::Pointer  affineTransform = AffineTransformType::New();
    affineTransform->SetCenter( rigidTransform->GetCenter() );
    affineTransform->SetTranslation( rigidTransform->GetTranslation() );
    affineTransform->SetMatrix( rigidTransform->GetMatrix() );
    registration->SetTransform( affineTransform );
    registration->SetInitialTransformParameters( affineTransform->GetParameters() );

    RegistrationType::ScheduleType fixedScheduleAffine(2,3);
    fixedScheduleAffine[0][0] = 4;
    fixedScheduleAffine[0][1] = 4;
    fixedScheduleAffine[0][2] = 4;
    fixedScheduleAffine[1][0] = 2;
    fixedScheduleAffine[1][1] = 2;
    fixedScheduleAffine[1][2] = 2;
    RegistrationType::ScheduleType movingScheduleAffine(2,3);
    movingScheduleAffine[0][0] = 4;
    movingScheduleAffine[0][1] = 4;
    movingScheduleAffine[0][2] = 4;
    movingScheduleAffine[1][0] = 2;
    movingScheduleAffine[1][1] = 2;
    movingScheduleAffine[1][2] = 2;

    registration->SetSchedules(fixedScheduleAffine, movingScheduleAffine);
    //registration->SetNumberOfLevels(2);

    // Define optimizer normaliztion to compensate for different dynamic range of rotations and translations.
    optimizerScales = OptimizerScalesType( affineTransform->GetNumberOfParameters() );
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = 1.0;
    optimizerScales[4] = 1.0;
    optimizerScales[5] = 1.0;
    optimizerScales[6] = 1.0;
    optimizerScales[7] = 1.0;
    optimizerScales[8] = 1.0;
    optimizerScales[9]  = translationScale;
    optimizerScales[10] = translationScale;
    optimizerScales[11] = translationScale;
    optimizer->SetScales( optimizerScales );

    optimizer->SetMaximumStepLength( AFFINE_MAX_STEP );
    optimizer->SetMinimumStepLength( AFFINE_MIN_STEP );
    //optimizer->SetMaximumStepLength( optimizer->GetCurrentStepLength() );
    //optimizer->SetMinimumStepLength( 10*optimizer->GetMinimumStepLength() );
    optimizer->SetNumberOfIterations( 2*AFFINE_ITERATIONS );

    //#ifdef MATTES
    // Start with N/(2*64) then inc by 8 (3 level pyramid
    //metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / 16) );
    // Hack as bug makes it run obs 1st it

    // Use mask region on Fixed Image
    if(useFixedMask)
    {
        metric->SetFixedImageMask( spatialObjectMask );
        if(numberOfVoxelsInMask >  numberOfPixels/4.0)
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4*16)) );
        else
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4)) );
    }
    else if(useMovingMask)
    {
        metric->SetMovingImageMask( spatialObjectMask );
        if(numberOfVoxelsInMask >  numberOfPixels/4.0)
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4*16)) );
        else
            metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4)) );
    }
    else
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / (256)) );

    //#endif
    //  metric->UseAllPixelsOn();
    //  metric->SetUseAllPixels(true);

    std::cout << "Starting Affine Registration " << std::endl;

    try {
        //itkProbesStart( "Affine Registration" );
#if ITK_VERSION_MAJOR < 4
        registration->StartRegistration();
#else
        registration->Update();
#endif
        //itkProbesStop( "Affine Registration" );
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return false;
    }

    std::cout << "Affine Registration completed" << std::endl;
    std::cout << std::endl;

    affineTransform->SetParameters( registration->GetLastTransformParameters() );

    /*********************************************************************************
     *                                                                               *
     *                    Output the Affine Parameters                               *
     *                                                                               *
     *********************************************************************************/

    /*
    if (writeFiles == true) {
        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"_affine.txt");
        sprintf(outputFileName, "%s%s", outputName.c_str(), "_affinetransform.txt");

        itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
        transformWriter->AddTransform(affineTransform);
        transformWriter->SetFileName(outputFileName);
        try {
            transformWriter->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return false;
        }
    }*/

    /*********************************************************************************
     *                                                                               *
     *                    Output the Affine Transformed Image                        *
     *                                                                               *
     *********************************************************************************/

    if (writeFiles == true) {
        sprintf(outputFileName,"%s", outputName.c_str());
        ResampleFilterType::Pointer resampleAffine = ResampleFilterType::New();
        resampleAffine->SetTransform( affineTransform );
        resampleAffine->SetInput( movingImage );
        resampleAffine->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleAffine->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleAffine->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleAffine->SetOutputDirection( fixedImage->GetDirection() );
        resampleAffine->SetDefaultPixelValue(0.0);
        resampleAffine->Update();

        std::cout << std::endl << "Writing  Affine resampled moving image..." << std::flush;
        writeFile3D( resampleAffine->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl << std::endl;
    }
    return true;
}

bool
Registration
::RegisterNonRigidBsplineITK(std::string dataFile, std::string atlasFile,
                             std::string outputName, bool /*Invert*/, std::string dataMaskFile,
                             std::string atlasMaskFile)
{
    bool writeFiles = false;
    char outputFileName[512];
    typedef double CoordinateRepType;

    typedef itk::AffineTransform< CoordinateRepType, SpaceDimension >                         AffineTransformType;
#if ITK_VERSION_MAJOR < 4
    typedef itk::BSplineDeformableTransform< CoordinateRepType, SpaceDimension, SplineOrder > DeformableTransformType;
#else
    typedef itk::BSplineTransform< CoordinateRepType, SpaceDimension, SplineOrder > DeformableTransformType;
#endif

    typedef itk::LinearInterpolateImageFunction< MovingImageType, double>                     InterpolatorType;
    typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType >                   RegistrationType;
    typedef itk::RegularStepGradientDescentOptimizer                                          OptimizerType;
    typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType > MetricType;
    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >                       ResampleFilterType;
    typedef itk::RecursiveMultiResolutionPyramidImageFilter< FixedImageType, FixedImageType > PyramidType;

    OptimizerType::Pointer      optimizer              = OptimizerType::New();
    InterpolatorType::Pointer   interpolator           = InterpolatorType::New();
    RegistrationType::Pointer   image_registration     = RegistrationType::New();
    MetricType::Pointer         metric                 = MetricType::New();

    image_registration->SetMetric( metric );
    image_registration->SetOptimizer(     optimizer     );
    image_registration->SetInterpolator(  interpolator  );

    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );

    RegistrationIterationUpdate::Pointer observer1 = RegistrationIterationUpdate::New();
    image_registration->AddObserver( itk::IterationEvent(), observer1 );

    // Auxiliary identity transform.
    typedef itk::IdentityTransform<double,SpaceDimension> IdentityTransformType;
    IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

    // Read in the files and set them to the registration object
    std::cout << "Moving image " << atlasFile << std::endl;
    MovingImageType::Pointer movingImage = this->readFile3D( atlasFile );
    FixedImageType::Pointer fixedImage   = this->readFile3D( dataFile );

    // Read in ITK transform
    //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"affine.txt");
    sprintf(outputFileName, "%s%s", outputName.c_str(), "_affinetransform.txt");

    itk::TransformFileReader::Pointer affineReader = itk::TransformFileReader::New();
    affineReader->SetFileName(outputFileName);
    try {
        affineReader->Update();
    } catch(itk::ExceptionObject &ex)
    {
        std::cerr << "Failed reading Transform" << std::endl;
        // Output the error to the terminal
        std::cout << "Failed reading Transform... " << ex.GetDescription() << std::endl;
        return false;
    }

    typedef itk::AffineTransform< double, ImageDimension > AffineTransformType;
    AffineTransformType::Pointer affineTransform;

    /* copy the transform into a new AffineTransform object */
    if(!strcmp(affineReader->GetTransformList()->front()->GetNameOfClass(), "AffineTransform"))
    {
        affineTransform = static_cast<AffineTransformType*>(affineReader->GetTransformList()->front().GetPointer());
    }
    else
    {
        std::cerr << "Affine Transform Input is not of AffineTransformType" << std::endl;
        return false;
    }

    // Add mask
    typedef itk::ImageMaskSpatialObject< SpaceDimension >   MaskType;
    typedef itk::Image< unsigned char, SpaceDimension >     ImageMaskType;
    typedef itk::ImageFileReader< ImageMaskType >           MaskReaderType;

    int numberOfVoxelsInMask = 0;
    bool useFixedMask = false;
    if ( dataMaskFile.length() != 0)
        useFixedMask = true;


    MaskReaderType::Pointer fixedImageMaskReader  = MaskReaderType::New();
    MaskType::Pointer spatialObjectMask = MaskType::New();
    if(useFixedMask)
    {
        std::cout << "Using Fixed mask " << dataMaskFile <<  std::endl;
        fixedImageMaskReader->SetFileName( dataMaskFile );

        try {
            fixedImageMaskReader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        typedef itk::BinaryThresholdImageFilter<ImageMaskType,ImageMaskType> BinaryThresholdImageFilterType;
        BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
        threshold->SetInput(fixedImageMaskReader->GetOutput());
        threshold->SetLowerThreshold(1);
        threshold->SetUpperThreshold(255);
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(255);
        threshold->Update();

        spatialObjectMask->SetImage( threshold->GetOutput() );
        useFixedMask = true;

        typedef itk::LabelStatisticsImageFilter<FixedImageType, ImageMaskType> LabelStatisticsImageFilterType;
        LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
        labelStatisticsImageFilter->SetInput(fixedImage);
        labelStatisticsImageFilter->SetLabelInput(threshold->GetOutput());
        labelStatisticsImageFilter->Update();
        //for(int i = 50; i < 300; i++)
        numberOfVoxelsInMask = labelStatisticsImageFilter->GetCount(255);
        std::cout << "We have " << numberOfVoxelsInMask << " in fixed mask region" << std::endl;

    }
    bool useMovingMask = false;
    if ( atlasMaskFile.length() != 0)
        useMovingMask = true;
    MaskReaderType::Pointer movingImageMaskReader  = MaskReaderType::New();
    MaskType::Pointer movingSpatialObjectMask = MaskType::New();
    if(useMovingMask)
    {
        std::cout << "Using Moving mask: " << atlasMaskFile << std::endl;
        movingImageMaskReader->SetFileName( atlasMaskFile );
        try {
            movingImageMaskReader->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        typedef itk::BinaryThresholdImageFilter<ImageMaskType,ImageMaskType> BinaryThresholdImageFilterType;
        BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
        threshold->SetInput(movingImageMaskReader->GetOutput());
        threshold->SetLowerThreshold(1);
        threshold->SetUpperThreshold(255);
        threshold->SetOutsideValue(0);
        threshold->SetInsideValue(255);
        threshold->Update();

        movingSpatialObjectMask->SetImage( threshold->GetOutput() );
        useMovingMask = true;
    }

    /*********************************************************************************
     *                                                                               *
     *                     Perform Coarse Deformable Registration                    *
     *                                                                               *
     *********************************************************************************/

    DeformableTransformType::Pointer  bsplineTransformCoarse = DeformableTransformType::New();

    // Setup the resample pyramid schedule
    PyramidType::Pointer movingPyramid = PyramidType::New();
    PyramidType::Pointer fixedPyramid  = PyramidType::New();
    int levels = 3;
    fixedPyramid->SetNumberOfLevels( levels );
    movingPyramid->SetNumberOfLevels( levels );
    fixedPyramid->SetInput( readFile3D( dataFile ) );
    movingPyramid->SetInput( readFile3D( atlasFile ) );
    itk::Array2D<unsigned int> schedule = fixedPyramid->GetSchedule();
    //for (int levelCounter = 0; levelCounter < levels; levelCounter++) schedule[levelCounter][3] = 0;  // ensure that the 4th dimension is not rescaled
    //fixedPyramid->SetSchedule( schedule );
    //movingPyramid->SetSchedule( schedule );
    fixedPyramid->Update();
    movingPyramid->Update();

    fixedImage  = fixedPyramid->GetOutput( PYRAMID_RESOLUTION_ONE );
    movingImage = movingPyramid->GetOutput( PYRAMID_RESOLUTION_ONE );
    FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
    image_registration->SetFixedImage(  fixedImage   );
    image_registration->SetMovingImage(   movingImage  );
    //image_registration->SetNumberOfLevels(1);

    std::cout << "------ Multi level resolution one ------" << std::endl;
    FixedImageType::SizeType fixedImageSize       = fixedRegion.GetSize();
    FixedImageType::SpacingType fixedImageSpacing = fixedImage->GetSpacing();
    std::cout << "ImageSize: " << fixedImageSize[0] << ", " << fixedImageSize[1] << ", " << fixedImageSize[2] << std::endl;
    std::cout << "ImageSpacing: " << fixedImageSpacing[0] << ", " <<  fixedImageSpacing[1] << ", " <<  fixedImageSpacing[2] << std::endl;
    std::cout << "ImageOrigin: " << fixedImage->GetOrigin() << std::endl;

    unsigned int coarseSpacing = GRID_POINT_SPACING_ONE; // Set the initial gridpoint spacing to 10mm
    std::cout << "Coarse Grid Point Spacing (mm): " << coarseSpacing << std::endl;

#if ITK_VERSION_MAJOR < 4
    typedef DeformableTransformType::RegionType RegionType;
    RegionType bsplineRegion;
    RegionType::SizeType   gridSizeOnImage;
    RegionType::SizeType   gridBorderSize;
    RegionType::SizeType   totalGridSize;

    gridSizeOnImage[0] = (int)floor( fixedImageSize[0] * fixedImageSpacing[0] / coarseSpacing );
    gridSizeOnImage[1] = (int)floor( fixedImageSize[1] * fixedImageSpacing[1] / coarseSpacing );
    gridSizeOnImage[2] = (int)floor( fixedImageSize[2] * fixedImageSpacing[2] / coarseSpacing );
    std::cout << "Number of grid points: " << gridSizeOnImage[0] << ", " << gridSizeOnImage[1] << ", " << gridSizeOnImage[2] << ", " << std::endl;
    gridBorderSize.Fill( SplineOrder );    // Border for spline order = 3 ( 1 lower, 2 upper )
    totalGridSize = gridSizeOnImage + gridBorderSize;

    bsplineRegion.SetSize( totalGridSize );

    typedef DeformableTransformType::SpacingType SpacingType;
    SpacingType spacing = fixedImage->GetSpacing();

    typedef DeformableTransformType::OriginType OriginType;
    OriginType origin = fixedImage->GetOrigin();

    for(unsigned int r=0; r<ImageDimension; r++) {
        spacing[r] *= static_cast<double>( fixedImageSize[r] - 1)  /  static_cast<double>(gridSizeOnImage[r] - 1 );
    }

    FixedImageType::DirectionType gridDirection = fixedImage->GetDirection();
    SpacingType gridOriginOffset = gridDirection * spacing;
    OriginType gridOrigin = origin - gridOriginOffset;

    std::cout << "Direction " << gridDirection << std::endl;
    std::cout << "Spacing " << spacing << std::endl;
    std::cout << "GridOrigin " << gridOrigin << std::endl;

    bsplineTransformCoarse->SetGridSpacing( spacing );
    bsplineTransformCoarse->SetGridOrigin( gridOrigin );
    bsplineTransformCoarse->SetGridRegion( bsplineRegion );
    bsplineTransformCoarse->SetGridDirection( gridDirection );

    bsplineTransformCoarse->SetBulkTransform( affineTransform );

#else

    DeformableTransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
    DeformableTransformType::MeshSizeType             meshSize;
    for( unsigned int i=0; i < ImageDimension; i++ )
    {
        fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
                                     static_cast<double>(
                                         fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1 );
    }
//   unsigned int numberOfGridNodesInOneDimension = (int)floor( fixedImageSize[0] * fixedImageSpacing[0] / coarseSpacing );

    for( unsigned int i=0; i < ImageDimension; i++ )
    {
        meshSize[i] = (int)floor( fixedImageSize[i] * fixedImageSpacing[i] / coarseSpacing ) - SplineOrder;
    }
//   meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );
    bsplineTransformCoarse->SetTransformDomainOrigin( fixedImage->GetOrigin() );
    bsplineTransformCoarse->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
    bsplineTransformCoarse->SetTransformDomainMeshSize( meshSize );
    bsplineTransformCoarse->SetTransformDomainDirection( fixedImage->GetDirection() );
    //TODO Set Affine transform to initialise NR registration

#endif


    typedef DeformableTransformType::ParametersType  ParametersType;

    unsigned int numberOfBSplineParameters = bsplineTransformCoarse->GetNumberOfParameters();

    typedef OptimizerType::ScalesType       OptimizerScalesType;
    OptimizerScalesType optimizerScales = OptimizerScalesType( numberOfBSplineParameters );
    optimizerScales.Fill( 1.0 );
    optimizer->SetScales( optimizerScales );

    ParametersType initialDeformableTransformParameters( numberOfBSplineParameters );
    initialDeformableTransformParameters.Fill( 0.0 );
    bsplineTransformCoarse->SetParameters( initialDeformableTransformParameters );

    image_registration->SetInitialTransformParameters( bsplineTransformCoarse->GetParameters() );
    image_registration->SetTransform( bsplineTransformCoarse );

    //  Next we set the parameters of the RegularStepGradientDescentOptimizer object.
    optimizer->SetMaximumStepLength( NRR_MAX_STEP_ONE  );
    optimizer->SetMinimumStepLength(  NRR_MIN_STEP );
    optimizer->SetRelaxationFactor( NRR_RELAXATION_FACTOR);
    optimizer->SetNumberOfIterations( NRR_ITERATIONS_ONE );
    optimizer->SetGradientMagnitudeTolerance( NRR_GRADIENT_MAGNITUDE_TOLERANCE);

    //#ifdef MATTES
    const unsigned int numberOfPixels = fixedRegion.GetNumberOfPixels();
    if(useFixedMask)
    {
        metric->SetFixedImageMask( spatialObjectMask );
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (2*16)) );
    }
    else
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / (2*32)) );

    //#endif
    //metric->UseAllPixelsOn();
    std::cout << std::endl << "Starting Deformable Registration Coarse Grid" << std::endl;
    try {
        //itkProbesStart( "Deformable Registration Coarse" );
#if ITK_VERSION_MAJOR < 4
        image_registration->StartRegistration();
#else
        image_registration->Update();
#endif
        //itkProbesStop( "Deformable Registration Coarse" );
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return false;
    }
    std::cout << "Deformable Registration Coarse Grid completed" << std::endl << std::endl;

    OptimizerType::ParametersType finalParameters = image_registration->GetLastTransformParameters();
    bsplineTransformCoarse->SetParameters( finalParameters );

    if (writeFiles == true) {
        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"coarse_nrr.nii.gz");
        sprintf(outputFileName, "%s", outputName.c_str());

        ResampleFilterType::Pointer resampleNRR = ResampleFilterType::New();
        resampleNRR->SetTransform(bsplineTransformCoarse );
        resampleNRR->SetInput(  movingImage );
        resampleNRR->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleNRR->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleNRR->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleNRR->SetOutputDirection( fixedImage->GetDirection() );
        resampleNRR->SetDefaultPixelValue(0.0);
        std::cout << "Computing NRR output" << std::endl;
        resampleNRR->Update();

        std::cout << "Writing coarse NRR resampled moving image..." << std::flush;
        writeFile3D( resampleNRR->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl;

        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"coarse_Bsplines.txt");
        sprintf(outputFileName, "%s%s", outputName.c_str(), "_coarse_Bsplines_transform.txt");
        itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
        transformWriter->AddTransform(bsplineTransformCoarse);
        transformWriter->SetFileName(outputFileName);
        try {
            transformWriter->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return false;
        }

        //sprintf(outputFileName,"%s%s", outputPrefix.c_str(),"coarse_Bsplines.nii.gz");
        sprintf(outputFileName,"%s", outputName.c_str());
#if ITK_VERSION_MAJOR < 4
        typedef itk::Compose3DVectorImageFilter<DeformableTransformType::ImageType, FieldType> VectorComposerType;
        VectorComposerType::Pointer composer = VectorComposerType::New();
        composer->SetInput1(bsplineTransformCoarse->GetCoefficientImage()[0]);
        composer->SetInput2(bsplineTransformCoarse->GetCoefficientImage()[1]);
        composer->SetInput3(bsplineTransformCoarse->GetCoefficientImage()[2]);
        composer->Update();
        typedef itk::ImageFileWriter<VectorComposerType::OutputImageType> CoeffWriterType;
        CoeffWriterType::Pointer coeffWriter = CoeffWriterType::New();
        coeffWriter->SetInput(composer->GetOutput());
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#else
        typedef itk::TransformFileWriter TransformWriterType;
        TransformWriterType::Pointer coeffWriter = TransformWriterType::New();
        coeffWriter->SetInput(bsplineTransformCoarse);
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#endif

        sprintf(outputFileName,"%s%s", outputName.c_str(),"_coarse_DefField.nii.gz");

        FieldType::Pointer field = FieldType::New();
        field->SetRegions( fixedRegion );
        field->SetOrigin( fixedImage->GetOrigin() );
        field->SetSpacing( fixedImage->GetSpacing() );
        field->SetDirection( fixedImage->GetDirection() );
        field->Allocate();

        typedef itk::ImageRegionIterator< FieldType > FieldIterator;
        FieldIterator fi( field, fixedRegion );
        fi.GoToBegin();

        DeformableTransformType::InputPointType  fixedPoint;
        DeformableTransformType::OutputPointType movingPoint;
        FieldType::IndexType index;
        VectorType displacement;

        while( ! fi.IsAtEnd() ) {
            index = fi.GetIndex();
            field->TransformIndexToPhysicalPoint( index, fixedPoint );
            movingPoint = bsplineTransformCoarse->TransformPoint( fixedPoint );
            displacement = movingPoint - fixedPoint;
            fi.Set( displacement );
            ++fi;
        }

        std::cout << "Writing deformation field ..." << std::flush;
        //writeFile4D( convertFromVectorDFImage( field ), argv[6] );
        typedef itk::ImageFileWriter<FieldType> DFWriterType;
        DFWriterType::Pointer DFWriter = DFWriterType::New();
        DFWriter->SetInput(field);
        DFWriter->SetFileName( outputFileName );
        DFWriter->Update();
        std::cout << " Done!" << std::endl << std::endl;
    }

    /*********************************************************************************
     *                                                                               *
     *                     Perform Medium Deformable Registration                    *
     *                                                                               *
     *********************************************************************************/
    metric->SetUseExplicitPDFDerivatives( false );

    DeformableTransformType::Pointer  bsplineTransformMedium = DeformableTransformType::New();

    fixedImage  = fixedPyramid->GetOutput( PYRAMID_RESOLUTION_TWO );
    movingImage = movingPyramid->GetOutput( PYRAMID_RESOLUTION_TWO );
    fixedRegion = fixedImage->GetBufferedRegion();
    image_registration->SetFixedImage(  fixedImage   );
    image_registration->SetMovingImage(   movingImage  );
    //image_registration->SetNumberOfLevels(1);

    std::cout << "------ Multi level resolution two ------" << std::endl;
    fixedImageSize = fixedRegion.GetSize();
    fixedImageSpacing = fixedImage->GetSpacing();
    std::cout << "ImageSize: " << fixedImageSize[0] << ", " << fixedImageSize[1] << ", " << fixedImageSize[2] << std::endl;
    std::cout << "ImageSpacing: " << fixedImageSpacing[0] << ", " <<  fixedImageSpacing[1] << ", " <<  fixedImageSpacing[2] << std::endl;

    unsigned int mediumSpacing = GRID_POINT_SPACING_TWO;
    std::cout << "Medium Grid Point Spacing (mm): " << mediumSpacing << std::endl;

#if ITK_VERSION_MAJOR < 4

    RegionType::SizeType   gridMediumSizeOnImage;
    gridMediumSizeOnImage[0] = (int)floor(fixedImageSize[0] * fixedImageSpacing[0] / mediumSpacing);
    gridMediumSizeOnImage[1] = (int)floor(fixedImageSize[1] * fixedImageSpacing[1] / mediumSpacing);
    gridMediumSizeOnImage[2] = (int)floor(fixedImageSize[2] * fixedImageSpacing[2] / mediumSpacing);
    std::cout << "Number of grid points: " << gridMediumSizeOnImage[0] << ", " << gridMediumSizeOnImage[1] << ", " << gridMediumSizeOnImage[2] << ", " << std::endl;

    totalGridSize = gridMediumSizeOnImage + gridBorderSize;
    bsplineRegion.SetSize( totalGridSize );

    SpacingType spacingMedium = fixedImage->GetSpacing();
    OriginType  originMedium  = fixedImage->GetOrigin();

    for(unsigned int rh=0; rh<ImageDimension; rh++)
    {
        spacingMedium[rh] *= static_cast<double>(fixedImageSize[rh] - 1)  / static_cast<double>(gridMediumSizeOnImage[rh] - 1);
        originMedium[rh] -= spacingMedium[rh];
    }

    gridDirection = fixedImage->GetDirection();
    SpacingType gridOriginOffsetMedium = gridDirection * spacingMedium;
    OriginType gridOriginMedium = origin - gridOriginOffsetMedium;

    bsplineTransformMedium->SetGridSpacing( spacingMedium );
    bsplineTransformMedium->SetGridOrigin( gridOriginMedium );
    bsplineTransformMedium->SetGridRegion( bsplineRegion );
    bsplineTransformMedium->SetGridDirection( gridDirection );

    bsplineTransformMedium->SetBulkTransform( affineTransform );

#else

    for( unsigned int i=0; i < ImageDimension; i++ )
    {
        fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
                                     static_cast<double>(
                                         fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1 );
        meshSize[i] = fixedImage->GetSpacing()[i] * static_cast<double>(fixedImageSize[i] - 1)  / static_cast<double>((int)floor(fixedImageSize[i] * fixedImageSpacing[i] / mediumSpacing) - 1) - SplineOrder;
    }
    /*  meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );*/
    bsplineTransformMedium->SetTransformDomainOrigin( fixedImage->GetOrigin() );
    bsplineTransformMedium->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
    bsplineTransformMedium->SetTransformDomainMeshSize( meshSize );
    bsplineTransformMedium->SetTransformDomainDirection( fixedImage->GetDirection() );
    //TODO Set Affine transform to initialise NR registration

#endif


    numberOfBSplineParameters = bsplineTransformMedium->GetNumberOfParameters();
    ParametersType parametersMedium( numberOfBSplineParameters );
    parametersMedium.Fill( 0.0 );

    //  Now we need to initialize the BSpline coefficients of the higher resolution
    //  transform. This is done by first computing the actual deformation field
    //  at the higher resolution from the lower resolution BSpline coefficients.
    //  Then a BSpline decomposition is done to obtain the BSpline coefficient of
    //  the higher resolution transform.

#if ITK_VERSION_MAJOR < 4
    unsigned int counter = 0;

    for ( unsigned int k = 0; k < SpaceDimension; k++ ) {
        typedef DeformableTransformType::ImageType ParametersImageType;
        typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
        ResamplerType::Pointer upsampler = ResamplerType::New();

        typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
        FunctionType::Pointer function = FunctionType::New();

        upsampler->SetInput( bsplineTransformCoarse->GetCoefficientImage()[k] );
        upsampler->SetInterpolator( function );
        upsampler->SetTransform( identityTransform );
        upsampler->SetSize( bsplineTransformMedium->GetGridRegion().GetSize() );
        upsampler->SetOutputSpacing( bsplineTransformMedium->GetGridSpacing() );
        upsampler->SetOutputOrigin( bsplineTransformMedium->GetGridOrigin() );

        typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType> DecompositionType;
        DecompositionType::Pointer decomposition = DecompositionType::New();
        decomposition->SetSplineOrder( SplineOrder );
        decomposition->SetInput( upsampler->GetOutput() );
        decomposition->Update();

        ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

        // copy the coefficients into the parameter array
        typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
        Iterator it( newCoefficients, bsplineTransformMedium->GetGridRegion() );
        while ( !it.IsAtEnd() )
        {
            parametersMedium[ counter++ ] = it.Get();
            ++it;
        }
    }
    bsplineTransformMedium->SetParameters( parametersMedium );
#else


#endif

    optimizerScales = OptimizerScalesType( numberOfBSplineParameters );
    optimizerScales.Fill( 1.0 );
    optimizer->SetScales( optimizerScales );


    if (false) {
        sprintf(outputFileName,"%s%s", outputName.c_str(),"mediumMR_initial.nii.gz");

        ResampleFilterType::Pointer resampleNRR = ResampleFilterType::New();
        resampleNRR->SetTransform(bsplineTransformMedium);
        resampleNRR->SetInput(  movingImage );
        resampleNRR->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleNRR->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleNRR->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleNRR->SetOutputDirection( fixedImage->GetDirection() );
        resampleNRR->SetDefaultPixelValue(0.0);
        std::cout << "Computing NRR output" << std::endl;
        resampleNRR->Update();

        std::cout << "Writing medium NRR resampled moving image..." << std::flush;
        writeFile3D( resampleNRR->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl;

        sprintf(outputFileName, "%s%s", outputName.c_str(), "mediumBspline_initial.txt");

        itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
        transformWriter->AddTransform(bsplineTransformMedium);
        transformWriter->SetFileName(outputFileName);
        try {
            transformWriter->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }

        sprintf(outputFileName,"%s%s", outputName.c_str(),"mediumBspline_initial.nii.gz");

#if ITK_VERSION_MAJOR < 4
        typedef itk::Compose3DVectorImageFilter<DeformableTransformType::ImageType, FieldType> VectorComposerType;
        VectorComposerType::Pointer composer = VectorComposerType::New();
        composer->SetInput1(bsplineTransformMedium->GetCoefficientImage()[0]);
        composer->SetInput2(bsplineTransformMedium->GetCoefficientImage()[1]);
        composer->SetInput3(bsplineTransformMedium->GetCoefficientImage()[2]);
        composer->Update();
        typedef itk::ImageFileWriter<VectorComposerType::OutputImageType> CoeffWriterType;
        CoeffWriterType::Pointer coeffWriter = CoeffWriterType::New();
        coeffWriter->SetInput(composer->GetOutput());
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#else
        typedef itk::TransformFileWriter TransformWriterType;
        TransformWriterType::Pointer coeffWriter = TransformWriterType::New();
        coeffWriter->SetInput(bsplineTransformMedium);
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#endif



        sprintf(outputFileName, "%s%s", outputName.c_str(), "mediumDefField_initial.nii.gz");

        FieldType::Pointer field = FieldType::New();
        field->SetRegions( fixedRegion );
        field->SetOrigin( fixedImage->GetOrigin() );
        field->SetSpacing( fixedImage->GetSpacing() );
        field->SetDirection( fixedImage->GetDirection() );
        field->Allocate();

        typedef itk::ImageRegionIterator< FieldType > FieldIterator;
        FieldIterator fi( field, fixedRegion );
        fi.GoToBegin();

        DeformableTransformType::InputPointType  fixedPoint;
        DeformableTransformType::OutputPointType movingPoint;
        FieldType::IndexType index;
        VectorType displacement;

        while( ! fi.IsAtEnd() ) {
            index = fi.GetIndex();
            field->TransformIndexToPhysicalPoint( index, fixedPoint );
            movingPoint = bsplineTransformMedium->TransformPoint( fixedPoint );
            displacement = movingPoint - fixedPoint;
            fi.Set( displacement );
            ++fi;
        }

        std::cout << "Writing deformation field ..." << std::flush;
        //writeFile4D( convertFromVectorDFImage( field ), argv[6] );
        typedef itk::ImageFileWriter<FieldType> DFWriterType;
        DFWriterType::Pointer DFWriter = DFWriterType::New();
        DFWriter->SetInput(field);
        DFWriter->SetFileName( outputFileName );
        DFWriter->Update();
        std::cout << " Done!" << std::endl << std::endl;

        //ResampleFilterType::Pointer resampleNRR = ResampleFilterType::New();
        sprintf(outputFileName,"%s%s", outputName.c_str(),"mediumNRR_initial.nii.gz");

        resampleNRR->SetTransform( bsplineTransformMedium );
        resampleNRR->SetInput(  movingImage );
        resampleNRR->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleNRR->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleNRR->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleNRR->SetOutputDirection( fixedImage->GetDirection() );
        resampleNRR->SetDefaultPixelValue(0.0);
        std::cout << "Computing NRR output" << std::endl;
        resampleNRR->Update();

        std::cout << "Writing NRR resampled moving image..." << std::flush;
        writeFile3D( resampleNRR->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl;

    }


    //  We now pass the parameters of the medium resolution transform as the initial
    //  parameters to be used in a second stage of the registration process.

    std::cout << "Starting Registration with medium resolution transform" << std::endl;

    image_registration->SetInitialTransformParameters( bsplineTransformMedium->GetParameters() );
    image_registration->SetTransform( bsplineTransformMedium );

    optimizer->SetMaximumStepLength( NRR_MAX_STEP_TWO );
    optimizer->SetNumberOfIterations( NRR_ITERATIONS_TWO );

    if(useFixedMask)
    {
        metric->SetFixedImageMask( spatialObjectMask );
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (16)) );
    }
    else
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / (32)) );
    //#ifdef MATTES
    //  metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / 3) );
    //#endif

    try {
        //itkProbesStart( "Deformable Registration Medium" );
#if ITK_VERSION_MAJOR < 4
        image_registration->StartRegistration();
#else
        image_registration->Update();
#endif

        // itkProbesStop( "Deformable Registration Medium" );
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Deformable Registration Medium Grid completed" << std::endl << std::endl;

    finalParameters = image_registration->GetLastTransformParameters();
    bsplineTransformMedium->SetParameters( finalParameters );

    if (writeFiles == true) {
        sprintf(outputFileName, "%s%s", outputName.c_str(), "medium_nrr.nii.gz");

        ResampleFilterType::Pointer resampleNRR = ResampleFilterType::New();
        resampleNRR->SetTransform(bsplineTransformMedium);
        resampleNRR->SetInput(  movingImage );
        resampleNRR->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleNRR->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleNRR->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleNRR->SetOutputDirection( fixedImage->GetDirection() );
        resampleNRR->SetDefaultPixelValue(0.0);
        std::cout << "Computing NRR output" << std::endl;
        resampleNRR->Update();

        std::cout << "Writing medium NRR resampled moving image..." << std::flush;
        writeFile3D( resampleNRR->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl;

        sprintf(outputFileName, "%s%s", outputName.c_str(), "medium_Bspline.txt");

        itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
        transformWriter->AddTransform(bsplineTransformMedium);
        transformWriter->SetFileName(outputFileName);
        try {
            transformWriter->Update();
        } catch (itk::ExceptionObject & excp) {
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }

        sprintf(outputFileName,"%s%s", outputName.c_str(),"medium_Bspline.nii.gz");

#if ITK_VERSION_MAJOR < 4
        typedef itk::Compose3DVectorImageFilter<DeformableTransformType::ImageType, FieldType> VectorComposerType;
        VectorComposerType::Pointer composer = VectorComposerType::New();
        composer->SetInput1(bsplineTransformMedium->GetCoefficientImage()[0]);
        composer->SetInput2(bsplineTransformMedium->GetCoefficientImage()[1]);
        composer->SetInput3(bsplineTransformMedium->GetCoefficientImage()[2]);
        composer->Update();
        typedef itk::ImageFileWriter<VectorComposerType::OutputImageType> CoeffWriterType;
        CoeffWriterType::Pointer coeffWriter = CoeffWriterType::New();
        coeffWriter->SetInput(composer->GetOutput());
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#else
        typedef itk::TransformFileWriter TransformWriterType;
        TransformWriterType::Pointer coeffWriter = TransformWriterType::New();
        coeffWriter->SetInput(bsplineTransformMedium);
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#endif

        sprintf(outputFileName, "%s%s", outputName.c_str(), "medium_DefField.nii.gz");

        FieldType::Pointer field = FieldType::New();
        field->SetRegions( fixedRegion );
        field->SetOrigin( fixedImage->GetOrigin() );
        field->SetSpacing( fixedImage->GetSpacing() );
        field->SetDirection( fixedImage->GetDirection() );
        field->Allocate();

        typedef itk::ImageRegionIterator< FieldType > FieldIterator;
        FieldIterator fi( field, fixedRegion );
        fi.GoToBegin();

        DeformableTransformType::InputPointType  fixedPoint;
        DeformableTransformType::OutputPointType movingPoint;
        FieldType::IndexType index;
        VectorType displacement;

        while( ! fi.IsAtEnd() ) {
            index = fi.GetIndex();
            field->TransformIndexToPhysicalPoint( index, fixedPoint );
            movingPoint = bsplineTransformMedium->TransformPoint( fixedPoint );
            displacement = movingPoint - fixedPoint;
            fi.Set( displacement );
            ++fi;
        }

        std::cout << "Writing deformation field ..." << std::flush;
        //writeFile4D( convertFromVectorDFImage( field ), argv[6] );
        typedef itk::ImageFileWriter<FieldType> DFWriterType;
        DFWriterType::Pointer DFWriter = DFWriterType::New();
        DFWriter->SetInput(field);
        DFWriter->SetFileName( outputFileName );
        DFWriter->Update();
        std::cout << " Done!" << std::endl << std::endl;

    }


    /*********************************************************************************
     *                                                                               *
     *                     Perform Fine Deformable Registration                      *
     *                                                                               *
     *********************************************************************************/

    DeformableTransformType::Pointer  bsplineTransformFine = DeformableTransformType::New();

    //fixedImage  = fixedPyramid->GetOutput(PYRAMID_RESOLUTION_THREE) ;
    //movingImage = movingPyramid->GetOutput(PYRAMID_RESOLUTION_THREE) ;

    movingImage = readFile3D( atlasFile );
    fixedImage   = readFile3D( dataFile );

    fixedRegion = fixedImage->GetBufferedRegion();
    image_registration->SetFixedImage(  fixedImage   );
    image_registration->SetMovingImage(   movingImage  );
    //registration->SetNumberOfLevels(1);

    std::cout << "------ Multi level resolution three ------" << std::endl;
    fixedImageSize       = fixedRegion.GetSize();
    fixedImageSpacing = fixedImage->GetSpacing();
    std::cout << "ImageSize: " << fixedImageSize[0] << ", " << fixedImageSize[1] << ", " << fixedImageSize[2] << std::endl;
    std::cout << "ImageSpacing: " << fixedImageSpacing[0] << ", " <<  fixedImageSpacing[1] << ", " <<  fixedImageSpacing[2] << std::endl;

    unsigned int fineSpacing = GRID_POINT_SPACING_THREE; // Set the initial gridpoint spacing to 10mm
    std::cout << "Fine Grid Point Spacing (mm): " << fineSpacing << std::endl;

#if ITK_VERSION_MAJOR < 4

    RegionType::SizeType   gridHighSizeOnImage;
    gridHighSizeOnImage[0] = (int)floor( fixedImageSize[0] * fixedImageSpacing[0] / fineSpacing );
    gridHighSizeOnImage[1] = (int)floor( fixedImageSize[1] * fixedImageSpacing[1] / fineSpacing );
    gridHighSizeOnImage[2] = (int)floor( fixedImageSize[2] * fixedImageSpacing[2] / fineSpacing );
    std::cout << "Number of grid points: " << gridHighSizeOnImage[0] << ", " << gridHighSizeOnImage[1] << ", " << gridHighSizeOnImage[2] << ", " << std::endl;
    totalGridSize = gridSizeOnImage + gridBorderSize;
    totalGridSize = gridHighSizeOnImage + gridBorderSize;
    bsplineRegion.SetSize( totalGridSize );

    SpacingType spacingHigh = fixedImage->GetSpacing();
    OriginType  originHigh  = fixedImage->GetOrigin();

    for(unsigned int rh=0; rh<ImageDimension; rh++)
    {
        spacingHigh[rh] *= static_cast<double>(fixedImageSize[rh] - 1)  / static_cast<double>(gridHighSizeOnImage[rh] - 1);
        originHigh[rh] -= spacingHigh[rh];
    }

    gridDirection = fixedImage->GetDirection();
    SpacingType gridOriginOffsetHigh = gridDirection * spacingHigh;
    OriginType gridOriginHigh = origin - gridOriginOffsetHigh;

    bsplineTransformFine->SetGridSpacing( spacingHigh );
    bsplineTransformFine->SetGridOrigin( gridOriginHigh );
    bsplineTransformFine->SetGridRegion( bsplineRegion );
    bsplineTransformFine->SetGridDirection( gridDirection );
    bsplineTransformFine->SetBulkTransform( affineTransform );

#else

    for( unsigned int i=0; i < ImageDimension; i++ )
    {
        fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
                                     static_cast<double>(
                                         fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1 );
        meshSize[i] = (int)floor( fixedImageSize[0] * fixedImageSpacing[0] / fineSpacing ) - SplineOrder;
    }
    /*  meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );*/
    bsplineTransformFine->SetTransformDomainOrigin( fixedImage->GetOrigin() );
    bsplineTransformFine->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
    bsplineTransformFine->SetTransformDomainMeshSize( meshSize );
    bsplineTransformFine->SetTransformDomainDirection( fixedImage->GetDirection() );
    //TODO Set Affine transform to initialise NR registration

#endif

    numberOfBSplineParameters = bsplineTransformFine->GetNumberOfParameters();
    ParametersType parametersHigh( numberOfBSplineParameters );
    parametersHigh.Fill( 0.0 );

#if ITK_VERSION_MAJOR < 4
    // Now we need to initialize the BSpline coefficients of the higher resolution transform.
    counter = 0;

    for ( unsigned int k = 0; k < SpaceDimension; k++ )
    {

        typedef DeformableTransformType::ImageType ParametersImageType;
        typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
        ResamplerType::Pointer upsampler = ResamplerType::New();

        typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
        FunctionType::Pointer function = FunctionType::New();

        upsampler->SetInput( bsplineTransformMedium->GetCoefficientImage()[k] );
        upsampler->SetInterpolator( function );
        upsampler->SetTransform( identityTransform );
        upsampler->SetSize( bsplineTransformFine->GetGridRegion().GetSize() );
        upsampler->SetOutputSpacing( bsplineTransformFine->GetGridSpacing() );
        upsampler->SetOutputOrigin( bsplineTransformFine->GetGridOrigin() );

        typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType> DecompositionType;
        DecompositionType::Pointer decomposition = DecompositionType::New();

        decomposition->SetSplineOrder( SplineOrder );
        decomposition->SetInput( upsampler->GetOutput() );
        decomposition->Update();

        ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

        // copy the coefficients into the parameter array
        typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
        Iterator it( newCoefficients, bsplineTransformFine->GetGridRegion() );
        while ( !it.IsAtEnd() ) {
            parametersHigh[ counter++ ] = it.Get();
            ++it;
        }

    }
    bsplineTransformFine->SetParameters( parametersHigh );
#else

#endif

    optimizerScales = OptimizerScalesType( numberOfBSplineParameters );
    optimizerScales.Fill( 1.0 );
    optimizer->SetScales( optimizerScales );


    //  We now pass the parameters of the high resolution transform as the initial
    //  parameters to be used in a second stage of the registration process.

    std::cout << "Starting Registration with high resolution transform" << std::endl;

    image_registration->SetInitialTransformParameters( bsplineTransformFine->GetParameters() );
    image_registration->SetTransform( bsplineTransformFine );

    optimizer->SetMaximumStepLength( NRR_MAX_STEP_THREE  );
    optimizer->SetNumberOfIterations( NRR_ITERATIONS_THREE );

    if(useFixedMask)
    {
        metric->SetFixedImageMask( spatialObjectMask );
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfVoxelsInMask / (4)) );
    }
    else
        metric->SetNumberOfSpatialSamples( static_cast<int>(numberOfPixels / (4)) );

    try
    {
        //itkProbesStart( "Deformable Registration Fine" );
#if ITK_VERSION_MAJOR < 4
        image_registration->StartRegistration();
#else
        image_registration->Update();
#endif
        //itkProbesStop( "Deformable Registration Fine" );
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Deformable Registration Fine Grid completed" << std::endl;
    std::cout << std::endl << std::endl;

    finalParameters = image_registration->GetLastTransformParameters();
    bsplineTransformFine->SetParameters( finalParameters );

    //itkProbesReport( std::cout );
    /*********************************************************************************
     *                                                                               *
     *         Output the B-spline transform parameters & Coeff Image                *
     *                                                                               *
     *********************************************************************************/
    std::cout << "Write transform" << std::endl;
    std::cout << std::endl << std::endl;

    sprintf(outputFileName, "%s%s", outputName.c_str(), "fine_Bsplines.txt");

    itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
    transformWriter->AddTransform(bsplineTransformFine);
    transformWriter->SetFileName(outputFileName);
    try {
        transformWriter->Update();
    } catch (itk::ExceptionObject & excp) {
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }

    if (writeFiles == true) {
        sprintf(outputFileName,"%s%s", outputName.c_str(),"fine_Bsplines.nii.gz");

#if ITK_VERSION_MAJOR < 4
        typedef itk::Compose3DVectorImageFilter<DeformableTransformType::ImageType, FieldType> VectorComposerType;
        VectorComposerType::Pointer composer = VectorComposerType::New();
        composer->SetInput1(bsplineTransformFine->GetCoefficientImage()[0]);
        composer->SetInput2(bsplineTransformFine->GetCoefficientImage()[1]);
        composer->SetInput3(bsplineTransformFine->GetCoefficientImage()[2]);
        composer->Update();
        typedef itk::ImageFileWriter<VectorComposerType::OutputImageType> CoeffWriterType;
        CoeffWriterType::Pointer coeffWriter = CoeffWriterType::New();
        coeffWriter->SetInput(composer->GetOutput());
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#else
        typedef itk::TransformFileWriter TransformWriterType;
        TransformWriterType::Pointer coeffWriter = TransformWriterType::New();
        coeffWriter->SetInput(bsplineTransformFine);
        coeffWriter->SetFileName( outputFileName );
        coeffWriter->Update();
#endif

    }


    /*********************************************************************************
     *                                                                               *
     *                       Output the Deformation Field                            *
     *                                                                               *
     *********************************************************************************/

    if (writeFiles == true) {
        sprintf(outputFileName, "%s%s", outputName.c_str(), "fine_DefField.nii.gz");

        FieldType::Pointer field = FieldType::New();
        field->SetRegions( fixedRegion );
        field->SetOrigin( fixedImage->GetOrigin() );
        field->SetSpacing( fixedImage->GetSpacing() );
        field->SetDirection( fixedImage->GetDirection() );
        field->Allocate();

        typedef itk::ImageRegionIterator< FieldType > FieldIterator;
        FieldIterator fi( field, fixedRegion );
        fi.GoToBegin();

        DeformableTransformType::InputPointType  fixedPoint;
        DeformableTransformType::OutputPointType movingPoint;
        FieldType::IndexType index;
        VectorType displacement;

        while( ! fi.IsAtEnd() ) {
            index = fi.GetIndex();
            field->TransformIndexToPhysicalPoint( index, fixedPoint );
            movingPoint = bsplineTransformFine->TransformPoint( fixedPoint );
            displacement = movingPoint - fixedPoint;
            fi.Set( displacement );
            ++fi;
        }

        std::cout << "Writing deformation field ..." << std::flush;
        //writeFile4D( convertFromVectorDFImage( field ), argv[6] );
        typedef itk::ImageFileWriter<FieldType> DFWriterType;
        DFWriterType::Pointer DFWriter = DFWriterType::New();
        DFWriter->SetInput(field);
        DFWriter->SetFileName( outputFileName );
        DFWriter->Update();
        std::cout << " Done!" << std::endl << std::endl;
    }

    /*********************************************************************************
     *                                                                               *
     *                  Output the Non-rigid Transformed Image                       *
     *                                                                               *
     *********************************************************************************/

    if (writeFiles == true) {
        sprintf(outputFileName, "%s%s", outputName.c_str(), "fine_nrr.nii.gz");

        ResampleFilterType::Pointer resampleNRR = ResampleFilterType::New();
        resampleNRR->SetTransform( bsplineTransformFine );
        resampleNRR->SetInput(  movingImage );
        resampleNRR->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
        resampleNRR->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampleNRR->SetOutputSpacing( fixedImage->GetSpacing() );
        resampleNRR->SetOutputDirection( fixedImage->GetDirection() );
        resampleNRR->SetDefaultPixelValue(0.0);
        std::cout << "Computing NRR output" << std::endl;
        resampleNRR->Update();

        std::cout << "Writing NRR resampled moving image..." << std::flush;
        this->writeFile3D( resampleNRR->GetOutput() , outputFileName );
        std::cout << " Done!" << std::endl;

    }
    std::cout << "Finished ITK NRR" << std::endl;
    return true;
}

/**
 * writeFile outputs a 3D image
 * @param image an ITK smart pointer to the image
 * @param outputFileName the path of the file
 */
void
Registration
::writeFile3D(FixedImageType::Pointer image, std::string outputFileName) {
    typedef itk::ImageFileWriter< FixedImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(image);
    writer->SetFileName(outputFileName);
    try {
        writer->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "ExceptionObject caught:" << std::endl;
        std::cerr << err << std::endl;
        exit(1);
    }
    std::cout << "File written: " << outputFileName <<  std::endl;
}

/**
 * readFile3D reads in a 3D file of the specified name
 * @param fileName the path of the file
 * @return an ITK smart pointer to the image
 */
Registration::FixedImageType::Pointer
Registration
::readFile3D(std::string fileName) {
    typedef itk::ImageFileReader< FixedImageType >  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName);
    try {
        reader->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "ExceptionObject caught:" << std::endl;
        std::cerr << err << std::endl;
        exit(1);
    }
    return reader->GetOutput();
}


void
Registration
::Voting(std::vector<std::string> inputFilenames, std::string outputFilename, int modulo, int N, ConsensusMode consensusMode)
{
    std::string cmd_staple(STAPLE_EXEC);

    cmd_staple.append(" -m ");

    if(consensusMode == VOTING)
        cmd_staple.append(STAPLE_TYPE);
    else
        cmd_staple.append(STAPLE_TYPE2);

    cmd_staple.append(" -in ");

    /*1 create the source_files*/
    for (unsigned int i = 0; i < inputFilenames.size(); i++)
    {
        if((i+1)%N == static_cast<unsigned int>(modulo) )
        {
            cmd_staple.append(inputFilenames[i]);
            cmd_staple.append("  ");
        }
    }

    cmd_staple.append("-n 117 -outh ");
    //cmd_gm.append(" -iv 0 254 -ov 0 1 -outh ");

    cmd_staple.append(outputFilename);

    if (system(cmd_staple.c_str()) != 0 ) {
        std::cerr << "Error: External call to STAPLE was unsuccesful." << std::endl;
    }

}

bool
Registration
::PropagateITKRigid(std::string inputImage, std::string outputImageFileName, std::string targetFile, std::string inputTransform, InterpolateMode mode)
{

    std::cout << "Read rigid transform " << inputTransform << std::endl;
    // Read in ITK transform
    typedef itk::VersorRigid3DTransform< double >VersorRigid3DTransformType;
    VersorRigid3DTransformType::Pointer rigidVersorTransform;

    typedef itk::Euler3DTransform< double >Euler3DTransformType;
    Euler3DTransformType::Pointer rigidEulerTransform;

    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >  ResampleFilterType;
    ResampleFilterType::Pointer resampleAffine = ResampleFilterType::New();


    itk::TransformFileReader::Pointer rigidReader = itk::TransformFileReader::New();
    std::string ra = inputTransform;
    rigidReader->SetFileName(ra);
    try {
        rigidReader->Update();
    }
    catch(itk::ExceptionObject &ex) {
        std::cout << "Failed reading rigid " << std::endl;
        throw ex;
        return false;
    }
    if(!strcmp(rigidReader->GetTransformList()->front()->GetNameOfClass(), "VersorRigid3DTransform"))
    {
        rigidVersorTransform = static_cast<VersorRigid3DTransformType*>(rigidReader->GetTransformList()->front().GetPointer());
        resampleAffine->SetTransform( rigidVersorTransform );
    }
    else if(!strcmp(rigidReader->GetTransformList()->front()->GetNameOfClass(), "Euler3DTransform"))
    {
        rigidEulerTransform = static_cast<Euler3DTransformType*>(rigidReader->GetTransformList()->front().GetPointer());
        resampleAffine->SetTransform( rigidEulerTransform );
    }
    else
    {
        std::cerr << "Transform Input is not of VersorRigid3DTransform" << std::endl;
        return false;
    }

    FixedImageType::Pointer movingImage = this->readFile3D( targetFile );
    MovingImageType::Pointer fixedImage = readFile3D( inputImage );

    std::cout << fixedImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout <<fixedImage->GetSpacing()<< std::endl;


    if( mode == INTERP_NEAREST)
        resampleAffine->SetInterpolator( m_NearestNeighborOrientatedInterpolator );
    else if(mode == INTERP_LINEAR)
        resampleAffine->SetInterpolator( m_LinearOrientatedInterpolator );
    else if(mode == INTERP_BSPLINE)
    {
        m_BsplineOrientatedInterpolator->SetSplineOrder(3);
        resampleAffine->SetInterpolator( m_BsplineOrientatedInterpolator);
    }
    else
    {
        m_BsplineOrientatedInterpolator->SetSplineOrder(5);
        resampleAffine->SetInterpolator( m_BsplineOrientatedInterpolator);
    }


    resampleAffine->SetInput( movingImage );
    resampleAffine->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
    resampleAffine->SetOutputOrigin(  fixedImage->GetOrigin() );
    resampleAffine->SetOutputSpacing( fixedImage->GetSpacing() );
    resampleAffine->SetOutputDirection( fixedImage->GetDirection() );
    resampleAffine->SetDefaultPixelValue(0.0);
    try {
        resampleAffine->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "ExceptionObject caught:" << std::endl;
        std::cerr << err << std::endl;
        return false;
    }

    writeFile3D( resampleAffine->GetOutput() , outputImageFileName );
    return true;
}

bool
Registration
::PropagateITKAffine(std::string inputImage,
                     std::string outputImageFileName,
                     std::string targetFile,
                     std::string inputTransform,
                     InterpolateMode mode)
{

    std::cout << "Read transform " << inputTransform << std::endl;
    // Read in ITK transform
    typedef itk::AffineTransform< double, ImageDimension > AffineTransformType;
    AffineTransformType::Pointer affineTransform;

    itk::TransformFileReader::Pointer affineReader = itk::TransformFileReader::New();
    try {
        affineReader->SetFileName(inputTransform);
        affineReader->Update();
    } catch(itk::ExceptionObject &ex) {
        throw ex;
        return 0;
    }

    if(!strcmp(affineReader->GetTransformList()->front()->GetNameOfClass(), "AffineTransform")) {
        affineTransform = static_cast<AffineTransformType*>(affineReader->GetTransformList()->front().GetPointer());
    } else {
        std::cerr << "Affine Transform Input is not of AffineTransformType" << std::endl;
        return EXIT_FAILURE;
    }

    FixedImageType::Pointer fixedImage = this->readFile3D( inputImage );
    MovingImageType::Pointer movingImage = readFile3D( targetFile );

    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >  ResampleFilterType;
    ResampleFilterType::Pointer resampleAffine = ResampleFilterType::New();

    if( mode == INTERP_NEAREST)
        resampleAffine->SetInterpolator( m_NearestNeighborOrientatedInterpolator );
    else if(mode == INTERP_LINEAR)
        resampleAffine->SetInterpolator( m_LinearOrientatedInterpolator );
    else if(mode == INTERP_BSPLINE)
    {
        m_BsplineOrientatedInterpolator->SetSplineOrder(3);
        resampleAffine->SetInterpolator( m_BsplineOrientatedInterpolator);
    }
    else
    {
        m_BsplineOrientatedInterpolator->SetSplineOrder(5);
        resampleAffine->SetInterpolator( m_BsplineOrientatedInterpolator);
    }


    resampleAffine->SetTransform( affineTransform );
    resampleAffine->SetInput( movingImage );
    resampleAffine->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
    resampleAffine->SetOutputOrigin(  fixedImage->GetOrigin() );
    resampleAffine->SetOutputSpacing( fixedImage->GetSpacing() );
    resampleAffine->SetOutputDirection( fixedImage->GetDirection() );
    resampleAffine->SetDefaultPixelValue(0.0);
    try {
        resampleAffine->Update();
    } catch (itk::ExceptionObject & err) {
        std::cerr << "ExceptionObject caught:" << std::endl;
        std::cerr << err << std::endl;
        exit(1);
    }

    writeFile3D( resampleAffine->GetOutput() , outputImageFileName );
    return true;
}

bool
Registration
::PropagateITKAffineBspline(std::string inputImage, std::string targetFile, std::string outputFilename, std::string inputTransform, std::string inputBsplineTransform, InterpolateMode mode)
{

    std::cout << "About to propogate ITK NRR" << std::endl;
    // Read in ITK transform
    typedef itk::AffineTransform< double, ImageDimension > AffineTransformType;
    AffineTransformType::Pointer affineTransform;

    itk::TransformFileReader::Pointer affineReader = itk::TransformFileReader::New();
    try {
        affineReader->SetFileName(inputTransform);
        affineReader->Update();
    } catch(itk::ExceptionObject &ex) {
        throw ex;
        return 0;
    }

    if(!strcmp(affineReader->GetTransformList()->front()->GetNameOfClass(), "AffineTransform")) {
        affineTransform = static_cast<AffineTransformType*>(affineReader->GetTransformList()->front().GetPointer());
    } else {
        std::cerr << "Affine Transform Input is not of AffineTransformType" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "About to read NRR transform" << std::endl;
    /* read the bspline transform from a file */
    itk::TransformFileReader::Pointer nrrReader = itk::TransformFileReader::New();
    try {
        nrrReader->SetFileName(inputBsplineTransform);
        nrrReader->Update();
    } catch(itk::ExceptionObject &ex) {
        throw ex;
        return 0;
    }
    std::cout << "Check" << std::endl;
    /* copy the transform into a new B-Spline Transform object and set the Affine bulk transform */
#if ITK_VERSION_MAJOR < 4
    typedef itk::BSplineDeformableTransform< double, ImageDimension, SplineOrder > DeformableTransformType;
#else
    typedef itk::BSplineTransform< double, ImageDimension, SplineOrder > DeformableTransformType;
#endif

    DeformableTransformType::Pointer bSplineTransform;

    if (!strcmp(nrrReader->GetTransformList()->front()->GetNameOfClass(), "BSplineDeformableTransform")) {
        bSplineTransform  =  static_cast<DeformableTransformType*>(nrrReader->GetTransformList()->front().GetPointer());
    } else {
        std::cerr << "NRR BSpline Transform Input is not of BSplineDeformableTransformType" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Read transforms" << std::endl;

#if ITK_VERSION_MAJOR < 4
    bSplineTransform->SetBulkTransform( affineTransform );
#else
//TODO Set Affine transform to initialise NR registration
#endif

    std::cout << "Read Fixed " <<inputImage << std::endl;
    FixedImageType::Pointer fixedImage = this->readFile3D( inputImage );
    std::cout << "Read Moving " <<targetFile<< std::endl;
    MovingImageType::Pointer movingImage = readFile3D( targetFile );

    std::cout << "Resample " << std::endl;
    typedef itk::ResampleImageFilter< MovingImageType, FixedImageType >  ResampleFilterType;
    ResampleFilterType::Pointer resampleNRR = ResampleFilterType::New();

    if( mode == INTERP_NEAREST)
        resampleNRR->SetInterpolator( m_NearestNeighborOrientatedInterpolator );
    else if(mode == INTERP_LINEAR)
        resampleNRR->SetInterpolator( m_LinearOrientatedInterpolator );
    else if(mode == INTERP_BSPLINE)
    {
        m_BsplineOrientatedInterpolator->SetSplineOrder(3);
        resampleNRR->SetInterpolator( m_BsplineOrientatedInterpolator);
    }
    else
    {
        m_BsplineOrientatedInterpolator->SetSplineOrder(5);
        resampleNRR->SetInterpolator( m_BsplineOrientatedInterpolator);
    }


    resampleNRR->SetTransform( bSplineTransform );
    resampleNRR->SetInput(  movingImage );
    resampleNRR->SetSize(  fixedImage->GetLargestPossibleRegion().GetSize() );
    resampleNRR->SetOutputOrigin(  fixedImage->GetOrigin() );
    resampleNRR->SetOutputSpacing( fixedImage->GetSpacing() );
    resampleNRR->SetOutputDirection( fixedImage->GetDirection() );
    resampleNRR->SetDefaultPixelValue(0.0);
    resampleNRR->Update();
    std::cout << "Done " << std::endl;
    this->writeFile3D( resampleNRR->GetOutput() , outputFilename );

    return true;
}

}
