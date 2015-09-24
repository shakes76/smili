#ifndef MILXQTRegistrationAdvancedOptions_H
#define MILXQTRegistrationAdvancedOptions_H

#include "milxQtRegistrationParams.h"
#include "ui_registrationAdvancedOptions.h"
#include "milxQtMain.h"


/*!
    \class milxQtRegistrationAdvancedOptions
    \brief This class is the advanced options form for the registration
*/
class milxQtRegistrationAdvancedOptions : public QDialog
{
    Q_OBJECT

public:

    /*!
        \fn milxQtRegistrationAdvancedOptions::milxQtRegistrationAdvancedOptions(QDialog *theParent = 0)
        \brief The standard constructor
    */
    milxQtRegistrationAdvancedOptions(QDialog *theParent = 0);

    /*!
        \fn milxQtRegistrationAdvancedOptions::milxQtRegistrationAdvancedOptions()
        \brief The standard destructor
    */
    virtual ~milxQtRegistrationAdvancedOptions();

    /*!
        \fn milxQtRegistrationAdvancedOptions::getParamsF3DNifti()
        \brief Get the advanced parameter of a Nifti F3D registration
    */
    milxQtRegistrationParams getParamsF3DNifti();

    /*!
    \fn milxQtRegistrationAdvancedOptions::getParamsAladinNifti()
    \brief Get the advanced parameter of a Nifti Aladin registration
    */
    milxQtRegistrationParams getParamsAladinNifti();

    /*!
        \fn milxQtRegistrationAdvancedOptions::getParamsElastixAffine()
        \brief Get the advanced parameters of an Elastix Affine registration
    */
    milxQtRegistrationParams getParamsElastixAffine();

    /*!
        \fn milxQtRegistrationAdvancedOptions::getParamsElastixBSpline()
        \brief Get the advanced parameters of an Elastix BSpline registration
    */
    milxQtRegistrationParams getParamsElastixBSpline();



    /*!
        \fn milxQtRegistrationAdvancedOptions::reset(RegType algo)
        \brief Reset the advanced form to the algorithm in parameter
    */
    void reset(RegType algo);

public slots:

    /*!
        \fn milxQtRegistrationAdvancedOptions::accept()
        \brief Click on button Ok
    */
    void accept();

    /*!
        \fn milxQtRegistrationAdvancedOptions::reject()
        \brief Click on button Cancel
    */
    void reject();

    /*!
        \fn milxQtRegistrationAdvancedOptions::reject()
        \brief Browse button clicked
    */
    void browseBtnClicked(); //

    /*!
        \fn milxQtRegistrationAdvancedOptions::reject()
        \brief Clear button clicked
    */
    void clearFileBtnClicked();

protected:
    /*!
        \fn milxQtRegistrationAdvancedOptions::createConnections()
        \brief Create connections with UI
    */
    void createConnections();

    Ui::dlgRegistrationAdvancedOptions ui; //!< User interface
    RegType currentAlgo; //!< Current algo
};

#endif // MILXQTRegistrationAdvancedOptions_H
