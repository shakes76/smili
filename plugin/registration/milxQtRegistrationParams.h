// This file contains all the structures and constants shared between most files of the plugin
#ifndef MILXQTRegistrationParams_H
#define MILXQTRegistrationParams_H
#include <QString>

// Type of the registration F3DNifti or AladinNifti
typedef enum
{
    AffineItk,
    DemonItk,
    F3DNifti,
    AladinNifti,
    ElastixAffine,
    ElastixBSpline,
    None
} RegType;

class milxQtRegistrationParams
{

public:
    milxQtRegistrationParams();
    ~milxQtRegistrationParams();
    void reset();


    RegType type;
    QString referenceName;
    QString floatingName;
    QString outputName;

    // Parameters for nifti F3D
    QString defOutputName;
    QString cppOutputName;
    int maxit;
    float spacing[3];
    unsigned int ln;
    unsigned int lp;
    bool nopy;
    bool useSym;
    bool cpp2Def;

    // Aladin
    bool rigOnly;
    bool affineDirect;
    float percentBlock;

    // Elastix
    QString outputFolder;
    QString parametersTxt;
    QString parameterFile;
    bool customParameterFile;
};

#endif // MILXQTRegistrationParams_H