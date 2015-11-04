#include "milxQtRegistrationParams.h"

milxQtRegistrationParams::milxQtRegistrationParams()
{
    reset();
}


void milxQtRegistrationParams::reset()
{
    referenceName = "";
    floatingName = "";
    outputName = "";
    defOutputName = "";
    cppOutputName = "";
    outputFolder = "";
    parametersTxt = "";
    parameterFile = "";
    customParameterFile = false;
}

milxQtRegistrationParams::~milxQtRegistrationParams()
{

}