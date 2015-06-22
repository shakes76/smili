#ifndef MILXQTRegistrationWindow_H
#define MILXQTRegistrationWindow_H

#include "ui_registrationWindow.h"
#include "milxQtRegistrationAdvancedOptions.h"
#include "milxQtMain.h"
#include "milxQtImage.h"
#include "milxQtFile.h"
#include "milxQtRegistrationNiftiReg.h"
#include "milxQtRegistrationStructures.h"

// Structure containing all the informations about a registration
struct RegistrationParams {
  RegType algo;
  PARAMSF3D F3D;
  PARAMSALADIN Aladin;
  PARAMSCPP2DEF cpp2Def;
  bool useCpp2Def;
};

typedef QList<milxQtImage *> QImageList;

/*!
    \class milxQtRegistrationWindow
    \brief This class is the registration form and starts the registration
*/
class milxQtRegistrationWindow : public QDialog
{
    Q_OBJECT

public:
	milxQtRegistrationWindow(milxQtMain *theParent = 0);
	virtual ~milxQtRegistrationWindow();
	void setup(RegType regType);
	void removeFiles(RegistrationParams reg); // Remove files created during a registration

public slots:
	void accept(); // Click on button Ok
	void reject(); // Click on button Cancel
	void referenceComboChange(int newIndex); // Reference combo box changed
    void algoComboChange(int newIndex); // Algo combo box changed
	void cpp2defFinished(); // Cpp2Def is done
	void registrationFinished(); // Registration is done
	void advancedOptionsClicked(); // Open the advanced option windows

protected:
	// User interface
    Ui::dlgRegistrationWindow ui;

	// Main window of SMILIX
    milxQtMain *MainWindow;

	// NiftiReg implementation
	milxQtRegistrationNifti *niftiReg;

	// Advanced options window
	milxQtRegistrationAdvancedOptions *advancedOptionsWindow;

	// Create connections with UI
    void createConnections();

	// Get the list of handled images for the registration
	void getListOfHandledImages();

	// Create a temporary file and copy its path in the buffer
	bool createTmpFile(char buffer[FILENAME_MAX + 1]);

    // Return the id of the selected images (id in the QImageList imageList)
    QList<int> getSelectedImages();

	// List of opened images
	QImageList imageList;
	
    // Main registration function
    void registration();

    // Registration queue
	QList<RegistrationParams> regQueue;

    // Current registration in progress
	RegistrationParams currentReg;

	// Is there any work in process
	bool workInProcess;
};

#endif // MILXQTRegistrationWindow_H
