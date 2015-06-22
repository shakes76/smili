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
	milxQtRegistrationAdvancedOptions(milxQtMain *theParent = 0);
	virtual ~milxQtRegistrationAdvancedOptions();
	PARAMSF3D getParamsF3D();
	PARAMSALADIN getParamsAladin();

	void reset(RegType algo);

public slots:
	void accept(); // Click on button Ok
	void reject(); // Click on button Cancel


protected:
	// User interface
    Ui::dlgRegistrationAdvancedOptions ui;

	// Main window of SMILIX
    milxQtMain *MainWindow;

	// Current algo
	RegType currentAlgo;

	// Create connections with UI
    void createConnections();
};

#endif // MILXQTRegistrationAdvancedOptions_H
