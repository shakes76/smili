#ifndef MILXQTRegistrationRegType
#define MILXQTRegistrationRegType

// Type of the registration FFD or Affine
typedef enum
{
	FFD,
	Affine
} RegType;

#endif

#ifndef MILXQTRegistrationAdvancedOptions_H
#define MILXQTRegistrationAdvancedOptions_H

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

	double getSx();
	double getSy();
	double getSz();
	int getMaxItLevel();
	int getNbLevel();
	int getFirstLevels();
	bool getUsePyramidalApproach();
	bool getUseSymmetricApproach();

	void reset(RegType algo);

public slots:
	void accept(); // Click on button Ok
	void reject(); // Click on button Cancel


protected:
	// User interface
    Ui::dlgRegistrationAdvancedOptions ui;

	// Main window of SMILIX
    milxQtMain *MainWindow;

	// Create connections with UI
    void createConnections();
};

#endif // MILXQTRegistrationAdvancedOptions_H
