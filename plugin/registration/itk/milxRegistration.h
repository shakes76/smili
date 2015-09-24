/*=========================================================================
Program: milxAliBaba
Module: milxRegistration.h
Author: ?
Modified by: Nicholas Dowson
Language: C++
Created: ?

Copyright: (c) 2009 CSIRO, Australia.

This software is protected by international copyright laws.
Any unauthorised copying, distribution or reverse engineering is prohibited.

Licence:
All rights in this Software are reserved to CSIRO. You are only permitted
to have this Software in your possession and to make use of it if you have
agreed to a Software License with CSIRO.

BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef __milxRegistration_h
#define __milxRegistration_h

#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#if ITK_VERSION_MAJOR < 4
#include <itkOrientedImage.h>
#include <itkOrientationAdapter.h>
#else
#include <itkImage.h>
#include <itkOrientationAdapterBase.h>
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <zlib.h>
#include <itkOrientImageFilter.h>


namespace milx
{

/** \class milxRegistration
 * \brief Library of Registration algorithms used in milx-view
 * - Currently supports the following algorithms.
 * ITK (rigid-affine)
 * ITK (Bspline deformable NRR)
 * ITK (Diffeomorphic Demons)
 *  - Several preprocessing steps available.
 * Histogram matching.
 * Zero-mean unit variance.
 *  - Output transforms saved.
 * Use milxPropagate to apply transforms.
 *  - Parameters can be set manually or via a config file.
 *
 *  - In the end Factories should be used to perform the below... too lazy now to do this.
 *
 *  - All data read in and converted to RAI...
 *  - Will need to add code to ensure orientation is taken into account and preserved...
 *  - Hence, data saved out in original space...
 *
 */

class Registration
{

public:
    Registration();
    ~Registration();

    //ND
    //typedef unsigned int RegistrationMode;

    typedef float                           InputPixelType;
    typedef float                           OutputPixelType;
    typedef itk::Vector<InputPixelType, 3 > VectorPixelType;
    typedef itk::Image<VectorPixelType, 3 > DemonsDeformationFieldType;

    typedef itk::Image<InputPixelType,3>    InputImageType;
    typedef itk::Image<OutputPixelType,3>   OutputImageType;


    typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, double >  NearestNeighborInterpolatorType;
    typedef itk::LinearInterpolateImageFunction<InputImageType, double >           LinearInterpolatorType;
    typedef itk::BSplineInterpolateImageFunction<InputImageType, double, double >  BsplineInterpolatorType;

    // Used by ITK NMI
    typedef itk::Vector< float, 3 >                  VectorType;
#if ITK_VERSION_MAJOR < 4
    typedef itk::OrientedImage< InputPixelType, 3>   FixedImageType;
    typedef itk::OrientedImage< InputPixelType, 3 >  MovingImageType;
    typedef itk::OrientedImage< VectorType, 3 >      FieldType;
#else
    typedef itk::Image< InputPixelType, 3>   FixedImageType;
    typedef itk::Image< InputPixelType, 3 >  MovingImageType;
    typedef itk::Image< VectorType, 3 >      FieldType;
#endif

    typedef itk::NearestNeighborInterpolateImageFunction<FixedImageType, double >  NearestNeighborOrientatedInterpolatorType;
    typedef itk::LinearInterpolateImageFunction<FixedImageType, double >           LinearOrientatedInterpolatorType;
    typedef itk::BSplineInterpolateImageFunction<FixedImageType, double, double >  BsplineOrientatedInterpolatorType;

    typedef enum {
        ALADIN_RIGID=0,           // aladin2000  Rigid
        ITK_RIGID=1,              // ITK NMI     Rigid
        ALADIN_AFFINE=2,          // aladin2000  Affine
        ITK_AFFINE=3,             // ITK NMI     Affine
        ALADIN_DIFFDEMONS=4,      // aladin2000  Affine + Demons NRR
        ITK_AFFINEBSPLINE=5,      // ITK NMI     Affine + BSpline NRR
        ALIBABA_RIGID=6,          // milxAliBaba Rigid
        ALIBABA_AFFINE=7,         // milxAliBaba Affine
        ALADIN_MILXBSPLINE=8,     // aladin2000  Affine + BSpline NRR
        NIFTI_REG=9,              // NiftyReg    Affine + NiftyReg NRR
        NIFTI_ALADIN_REG=10,      // NiftyReg    Affine
        ITK_AFFINE_DIFFDEMONS=11, // ITK NMI     Affine + Demons NRR
        ITKALADIN_DIFFDEMONS=12,  // ITK aladin  Affine + Demons NRR
        DIFFDEMONS=13             // ITK aladin  Affine + Demons NRR
    }                                      RegistrationMode;

    typedef enum { INTERP_NEAREST=0, INTERP_LINEAR=1, INTERP_BSPLINE=2, INTERP_BSPLINE5=3 }   InterpolateMode;
    typedef enum { VOTING=0, MULTISTAPLE=1 }                               ConsensusMode;

    void writeFile3D(FixedImageType::Pointer image, std::string outputFileName);
    FixedImageType::Pointer readFile3D(std::string fileName);
    /**
     * Configure registration class
     *  - Allows everything to be configured via one txt file.
     *  - ie load IDs, images, atlases, labellings, RegistrationMode, Registration Setting filename
     */
    bool LoadSettings(std::string /*filename*/)
    {
        return false;
    }

    /**
     * Various Get Calls
     */
    int GetNumberOfCases()
    {
        return m_CaseFilenames.size();
    }

    int GetNumberOfCases2()
    {
        return m_CaseFilenames2.size();
    }

    int GetNumberOfAtlases()
    {
        return m_AtlasFilenames.size();
    }

    int GetNumberOfLabellings()
    {
        return m_LabellingFilenames.size();
    }

    std::vector<int> GetCaseIDs()
    {
        return m_CaseIDs;
    }

    std::vector<std::string> GetCaseFileNames()
    {
        return m_CaseFilenames;
    }

    std::vector<std::string> GetCaseFileNames2()
    {
        return m_CaseFilenames2;
    }

    std::string GetCaseFileName(int caseID)
    {
        for(int i = 0; i < this->GetNumberOfCases(); i++)
        {
            if(m_CaseIDs[i] == caseID)
                return m_CaseFilenames[i];
        }
        std::cerr << "Warning caseID " << caseID << " not found ?" <<std::endl;
        return NULL;
    }

    std::string GetCaseFileName2(int caseID)
    {
        if(m_CaseFilenames2.size() == m_CaseFilenames.size())
        {
            std::cerr << "Bailing as File2size != File1size" <<std::endl;
            return NULL;
        }
        for(int i = 0; i < this->GetNumberOfCases(); i++)
        {
            if(m_CaseIDs[i] == caseID)
                return m_CaseFilenames2[i];
        }
        std::cerr << "Warning caseID " << caseID << " not found ?" <<std::endl;
        return NULL;
    }

    std::vector<int> GetAtlasIDs()
    {
        return m_AtlasIDs;
    }
    std::vector<std::string> GetAtlasFileNames()
    {
        return m_AtlasFilenames;
    }
    std::string GetAtlasFileName(int caseID)
    {
        for(int i = 0; i < this->GetNumberOfAtlases(); i++)
        {
            if(m_AtlasIDs[i] == caseID)
                return m_AtlasFilenames[i];
        }
        std::cerr << "Warning caseID " << caseID << " not found ?" <<std::endl;
        return NULL;
    }

    std::vector<int> GetLabellingIDs()
    {
        return m_LabellingIDs;
    }
    std::vector<std::string> GetLabellingFileNames()
    {
        std::cout << "Get Labelled filenames " << m_LabellingFilenames.size() << std::endl;
        return m_LabellingFilenames;
    }
    std::vector<std::string> GetLabellingFileNames(int caseID)
    {
        std::vector<std::string> labellingNames;
        for(int i = 0; i < this->GetNumberOfLabellings(); i++)
        {
            if(m_LabellingIDs[i] == caseID)
                labellingNames.push_back(m_LabellingFilenames[i]);
        }
        return labellingNames;
    }


    /**
     * Return the list of final output filenames for the propagated atlases and labels
     *  - As multiple
     */
    std::vector<std::string> GetPropagatedOutputAtlases()
    {

        return m_PropagatedAtlases;
    }
    std::vector<std::string> GetPropagatedOutputTransforms()
    {
        return m_PropagatedTransforms;
    }
    std::string GetPropagatedOutputAtlases(int caseID)
    {
        for(int i = 0; i < this->GetNumberOfAtlases(); i++)
        {
            if(m_AtlasIDs[i] == caseID)
                return m_PropagatedAtlases[i];
        }
        std::cerr << "Warning caseID " << caseID << " not found when returning propagated atlases ?" <<std::endl;
        return NULL;
    }

    std::vector<std::string> GetPropagatedOutputLabels()
    {
        return m_PropagatedLabelling;
    }
    std::vector<std::string> GetPropagatedOutputLabels(int caseID)
    {
        std::vector<std::string> labellingNames;
        for(int i = 0; i < this->GetNumberOfLabellings(); i++)
        {
            if(m_LabellingIDs[i] == caseID)
                labellingNames.push_back(m_PropagatedLabelling[i]);
        }
        return labellingNames;
    }
    bool ClearPropagatedOutputAtlasesLabels()
    {
        m_PropagatedAtlases.clear();
        m_PropagatedLabelling.clear();
        return true;
    }


    /**
     * Set/Clear list of images to register or co-register
     */
    bool SetCase(int caseID, std::string filenames)
    {
        std::vector<int> ids;
        ids.push_back(caseID);
        std::vector<std::string> caseFilenames;
        caseFilenames.push_back(filenames);
        return this->SetCases(ids, caseFilenames);
    }
    bool SetCases(std::vector<int> caseIDs, std::vector<std::string> caseFilenames)
    {
        m_CaseIDs = caseIDs;
        m_CaseFilenames = caseFilenames;
        return true;
    }
    bool SetCasePair(int caseID, std::string filenames, std::string filenames2)
    {
        std::vector<int> ids;
        ids.push_back(caseID);
        std::vector<std::string> caseFilenames;
        caseFilenames.push_back(filenames);
        std::vector<std::string> caseFilenames2;
        caseFilenames2.push_back(filenames2);
        return this->SetCasesPair(ids, caseFilenames, caseFilenames2);
    }
    bool SetCasesPair(std::vector<int> caseIDs, std::vector<std::string> caseFilenames, std::vector<std::string> caseFilenames2)
    {
        m_CaseIDs = caseIDs;
        m_CaseFilenames = caseFilenames;
        m_CaseFilenames2 = caseFilenames2;
        return true;
    }
    bool ClearCases()
    {
        m_CaseIDs.clear();
        m_CaseFilenames.clear();
        return true;
    }

    /**
     * Set/Clear list of atlases, labellings and IDs
     *  - Case IDs should be unique.
     *  - Note: Labelling IDs are assumed to correspond to caseIDs
     *  - Previous Atlases will be removed.
     */
    bool SetAtlases(std::vector<int> atlasIDs, std::vector<std::string> atlasFilenames)
    {
        m_AtlasIDs = atlasIDs;
        m_AtlasFilenames = atlasFilenames;
        return true;
    }
    bool SetAtlas(int atlasIDs, std::string atlasFilenames)
    {
        m_AtlasIDs.push_back(atlasIDs);
        m_AtlasFilenames.push_back(atlasFilenames);
        return true;
    }

    bool SetAtlases(std::vector<int> atlasIDs, std::vector<std::string> atlasFilenames, std::vector<int> labellingIDs, std::vector<std::string> labellingFilenames)
    {
        m_AtlasIDs = atlasIDs;
        m_AtlasFilenames = atlasFilenames;
        m_LabellingIDs = labellingIDs;
        m_LabellingFilenames = labellingFilenames;
        std::cout << "Labelling filename " << m_LabellingFilenames.size() << std::endl;

        return true;
    }

    bool SetAtlases(std::vector<int> atlasIDs,
                    std::vector<std::string> atlasFilenames,
                    std::vector<int> labellingIDs,
                    std::vector<std::string> labellingFilenames,
                    std::vector<int> maskIDs,
                    std::vector<std::string> maskFilenames)
    {
        m_AtlasIDs = atlasIDs;
        m_AtlasFilenames = atlasFilenames;
        m_LabellingIDs = labellingIDs;
        m_LabellingFilenames = labellingFilenames;
        std::cout << "Labelling filename " << m_LabellingFilenames.size() << std::endl;

        m_MaskIDs = maskIDs;
        m_MaskFilenames= maskFilenames;
        std::cout << "Number of image masks used " << m_MaskFilenames.size() << std::endl;

        return true;
    }

    bool SetAtlases(std::vector<int> atlasIDs,
                    std::vector<std::string> atlasFilenames,
                    std::vector<int> labellingIDs,
                    std::vector<std::string> labellingFilenames,
                    std::vector<int> maskIDs,
                    std::vector<std::string> maskFilenames,
                    std::vector<std::string> maskFilenames2)
    {
        m_AtlasIDs = atlasIDs;
        m_AtlasFilenames = atlasFilenames;
        m_LabellingIDs = labellingIDs;
        m_LabellingFilenames = labellingFilenames;
        std::cout << "Labelling filename " << m_LabellingFilenames.size() << std::endl;

        m_MaskIDs = maskIDs;
        m_MaskFilenames= maskFilenames;
        m_MaskFilenames2= maskFilenames2;
        std::cout << "Number of image masks used " << m_MaskFilenames.size() << std::endl;
        std::cout << "Number of image masks2 used " << m_MaskFilenames2.size() << std::endl;

        return true;
    }

    bool SetAtlas(int atlasIDs, std::string atlasFilenames, std::vector<std::string> labellingFilenames)
    {
        m_AtlasIDs.push_back(atlasIDs);
        m_AtlasFilenames.push_back(atlasFilenames);
        for (unsigned int i=0; i<labellingFilenames.size(); i++)
            m_LabellingIDs.push_back(atlasIDs);
        m_LabellingFilenames = labellingFilenames;
        std::cout << "Labelling filename " << m_LabellingFilenames.size() << " " << m_LabellingIDs.size() << std::endl;
        return true;
    }

    bool ClearAtlases()
    {
        m_AtlasIDs.clear();
        m_AtlasFilenames.clear();
        m_LabellingIDs.clear();
        m_LabellingFilenames.clear();
        return true;
    }

    /**
     * Set Output prefix
     *  - Note: this supports the use of /tmp/CASEID/case_CASEID_with_ATLASID_
     */
    void SetOutputPrefix(std::string prefix)
    {
        m_OutputPrefix = prefix;
    }
    std::string GetOutputPrefix()
    {
        return m_OutputPrefix;
    }

    /**
     * SetInterpolation used for each step
     */
    // Used inside registration algorithm (only used in ITK methods)
    bool SetRegistrationInterpolation(InterpolateMode mode)
    {
        m_RegistrationInterpolation = mode;
        return true;
    }
    // Used to propogate atlas
    bool SetAtlasInterpolation(InterpolateMode mode)
    {
        m_AtlasInterpolation = mode;
        return true;
    }
    // Used to propogate labels
    bool SetLabelledImageInterpolation(InterpolateMode mode)
    {
        m_LabelledImageInterpolation = mode;
        return true;
    }
    bool SetLabelledImageInterpolation(std::vector<InterpolateMode> modes)
    {
        m_MultipleLabelledImageInterpolation = modes;
        return true;
    }

    /**
     * Enable/Disable the use of Voting
     *  - By default disabled
     */
    void SetUseConsensus(bool flag)
    {
        m_UseConsensus = flag;
    }
    bool GetUseConsensus()
    {
        return m_UseConsensus;
    }

    /**
     * Set Registration scheme to use.
     */
    bool SetRegistrationMode(int registrationScheme, std::string configFilename = "NULL")
    {
        std::cout << "milxRegistration::SetRegistrationModeSet " << registrationScheme << std::endl;
        if(registrationScheme == 0)
        {
            RegistrationMode mode = ALADIN_RIGID;
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 1)
        {
            RegistrationMode mode = ITK_RIGID;
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 2)
        {
            RegistrationMode mode = ALADIN_AFFINE;
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 3)
        {
            RegistrationMode mode = ITK_AFFINE;
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 4)
        {
            RegistrationMode mode = ALADIN_DIFFDEMONS;
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 5)
        {
            RegistrationMode mode = ITK_AFFINEBSPLINE;
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 6)
        {
            RegistrationMode mode = ALIBABA_RIGID; //ND
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 7)
        {
            RegistrationMode mode = ALIBABA_AFFINE; //ND
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 8)
        {
            RegistrationMode mode = ALADIN_MILXBSPLINE; //JF
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 9)
        {
            RegistrationMode mode = NIFTI_REG; //PB
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 10)
        {
            RegistrationMode mode = NIFTI_ALADIN_REG; //PB
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 11)
        {
            RegistrationMode mode = ITK_AFFINE_DIFFDEMONS; //PB
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 12)
        {
            RegistrationMode mode = ITKALADIN_DIFFDEMONS; //PB
            this->SetRegistrationMode(mode, configFilename);
        }
        else if(registrationScheme == 13)
        {
            RegistrationMode mode = DIFFDEMONS; //DRH
            this->SetRegistrationMode(mode, configFilename);
        }
        else
        {
            std::cerr << "Unknown registration mode: " << registrationScheme << std::endl;
            return false;
        }

        return true;
    }
    bool SetRegistrationMode(RegistrationMode registrationScheme,
                             std::string /*configFilename = "NULL"*/ )
    {
        m_RegistrationMode = registrationScheme;
        return true;
    }
    RegistrationMode GetRegistrationMode()
    {
        return this->m_RegistrationMode;
    }

    /**
     * Set Consensus Mode
     */
    bool SetConsensusMode(int consensusScheme, std::string configFilename = "NULL")
    {
        std::cout << "Set Mode " << consensusScheme << std::endl;
        if(consensusScheme == 0)
        {
            ConsensusMode mode = VOTING;
            this->SetConsensusMode(mode, configFilename);
        }
        else if(consensusScheme == 1)
        {
            ConsensusMode mode = MULTISTAPLE;
            this->SetConsensusMode(mode, configFilename);
        }
        return true;
    }
    bool SetConsensusMode(ConsensusMode consensusScheme,
                          std::string /*configFilename = "NULL" */ )
    {
        m_ConsensusMode = consensusScheme;
        return true;
    }
    ConsensusMode GetConsensusMode()
    {
        return this->m_ConsensusMode;
    }

    /**
     * Options to configure individual schemes (ie if not using config file).
     */
    bool ConfigureRegistration(RegistrationMode /*mode*/, std::string /*filename*/)
    {
        return true;
    }

    /**
     * Options to enable the amount of files to remain
     *  - 5 or more saves all files (default)
     *  - 4 do not save deformation fields
     *  - 3 do not save propagated atlases
     *  - 2 do not save transforms
     *  - 1 do not save propagated labels
     *  - 0 should save nothing
     */
    void SetFilesToSave(int value)
    {
        m_SaveFiles = value;
    }

    /**
     * Not used in this code
     */
    void UpdatePipeline() {};

    /**
     * Variables used by this class;
     */
    std::vector<int>         m_AtlasIDs;
    std::vector<std::string> m_AtlasFilenames;
    std::vector<int>         m_LabellingIDs;
    std::vector<std::string> m_LabellingFilenames;

    std::vector<int>         m_MaskIDs;
    std::vector<std::string> m_MaskFilenames;
    std::vector<std::string> m_MaskFilenames2;

    std::vector<int>         m_CaseIDs;
    std::vector<std::string> m_CaseFilenames;
    std::vector<std::string> m_CaseFilenames2;

    std::string              m_OutputPrefix;

    std::vector<std::string> m_PropagatedAtlases;
    std::vector<std::string> m_PropagatedLabelling;
    std::vector<std::string> m_PropagatedTransforms;

    InterpolateMode              m_RegistrationInterpolation;
    InterpolateMode              m_AtlasInterpolation;
    InterpolateMode              m_LabelledImageInterpolation;
    std::vector<InterpolateMode> m_MultipleLabelledImageInterpolation;

    RegistrationMode             m_RegistrationMode;

    bool                         m_UseConsensus;
    ConsensusMode                m_ConsensusMode;

    int m_SaveFiles;
    /**
     * Declare interpolators
     */
    NearestNeighborInterpolatorType::Pointer m_NearestNeighborInterpolator;
    LinearInterpolatorType::Pointer          m_LinearInterpolator;
    BsplineInterpolatorType::Pointer         m_BsplineInterpolator;

    NearestNeighborOrientatedInterpolatorType::Pointer m_NearestNeighborOrientatedInterpolator;
    LinearOrientatedInterpolatorType::Pointer          m_LinearOrientatedInterpolator;
    BsplineOrientatedInterpolatorType::Pointer         m_BsplineOrientatedInterpolator;

    void Voting(std::vector<std::string> inputFilenames,
                std::string outputFilename,
                int modulo,
                int N,
                ConsensusMode consensusMode);

    /**
     * Various code that requires refactoring
     *  - Note: In the future, the below should probably be implemented in a factory !
     */

    // Warning: This code generally assumes that
    // datafile is fixed and atlas is moving


    bool RegisterRigidITK(std::string dataFile, std::string atlasFile, std::string outputPrefix, bool Invert, std::string dataMaskFile, std::string atlasMaskFile );
    bool RegisterMilxNonRigid(std::string dataFile, std::string atlasFile, std::string outputName, bool Invert);
    bool RegisterAffineITK(std::string dataFile, std::string atlasFile, std::string outputPrefix, bool Invert, std::string dataMaskFile, std::string atlasMaskFile );
    bool RegisterNonRigidBsplineITK(std::string dataFile, std::string atlasFile, std::string outputPrefix, bool Invert, std::string dataMaskFile, std::string atlasMaskFile );


    void LoadTRSF(std::string filename1, std::string filename2, double * array, bool invert);
    void LoadTRSF(std::string filename, double * array, bool invert);

    bool PropagateAffine(std::string dataFile,
                         std::string outFile,
                         std::string targetFile,
                         double * array,
                         InterpolateMode mode,
                         itk::SpatialOrientation::ValidCoordinateOrientationFlags &flag,
                         itk::SpatialOrientationAdapter::DirectionType &dir);

    bool PropagateITKRigid(std::string dataFile, std::string outFile, std::string targetFile, std::string transformFile, InterpolateMode mode);
    bool PropagateITKAffine(std::string dataFile, std::string outFile, std::string targetFile, std::string transformFile, InterpolateMode mode);
    bool PropagateITKAffineBspline(std::string inputImage, std::string targetFile, std::string outputFilename, std::string inputTransform, std::string inputBsplineTransform, InterpolateMode mode);


    // Simple function to convert a string into an array of char
    // which can be passed as an input for the nifty-reg executables
    // which have been converted in libraries but still take (int argc, char *argv[]) as input
    unsigned int StringToCharArray(std::string & cmdline, char ** argument );

};

} // end namespace milx

#endif
