#ifndef MILXQTRegistrationAdvancedOptions_H
#define MILXQTRegistrationAdvancedOptions_H

#include "milxQtRegistrationStructures.h"
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
        \fn milxQtRegistrationAdvancedOptions::getParamsF3D()
        \brief Get the advanced parameter of a F3D registration
    */
	ParamsF3D getParamsF3D();
    
	/*!
        \fn milxQtRegistrationAdvancedOptions::getParamsAladin()
        \brief Get the advanced parameters of an Aladin registration
    */
	ParamsAladin getParamsAladin();

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
