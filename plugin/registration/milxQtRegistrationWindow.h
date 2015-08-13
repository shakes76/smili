#ifndef MILXQTRegistrationWindow_H
#define MILXQTRegistrationWindow_H

#include "ui_registrationWindow.h"
#include "milxQtRegistrationAdvancedOptions.h"
#include "milxQtMain.h"
#include "milxQtImage.h"
#include "milxQtFile.h"
#include "milxQtRegistrationNiftiReg.h"
#include "milxQtRegistrationImage.h"

typedef QList<milxQtImage *> QImageList;

/*!
    \class milxQtRegistrationWindow
    \brief This class is the registration form and starts the registration
*/
class milxQtRegistrationWindow : public QDialog
{
    Q_OBJECT

public:
	milxQtRegistrationWindow(QWidget * theParent);
	virtual ~milxQtRegistrationWindow();
	void initUI();
	void setAlgo(RegType regType);
	void updateImageListCombo();
	void updateOpenImages();
	bool isImageInList(QString path);
	QList<milxQtRegistrationImage *> images;
	void updateParameters(); // update images parameters
	ParamsF3D getParamsF3D();
	ParamsAladin getParamsAladin();
	void addImage(milxQtRegistrationImage *); // Add an image to the list
	void disableUI(); // Disable the user interface
	void enableUI(); // Enable the user interface
	void performRegistrations(); // perform the next registration
	void workCompleted(); // everything has been completed
	void computeAtlas(); // compute the average (Atlas) of all the registrations

public slots:
	void accept(); // Click on button Ok
	void reject(); // Click on button Cancel
	void referenceComboChange(int newIndex); // Reference combo box changed
    void algoComboChange(int newIndex); // Algo combo box changed
	void advancedOptionsClicked(); // Open the advanced option windows
	void addImageClicked(); // Button add image clicked
	void selectAllClicked(); // Button select all clicked
	void unselectAllClicked(); // Button unselect all clicked
	void browseBtnClicked(); // Button browse clicked
	void regComplete(); // A registration has been completed
	void clearList(); // Clear the list of images
	void averageComputed(); // The average (Atlas) has been computed

protected:
	// User interface
    Ui::dlgRegistrationWindow ui;

	// Main window of SMILIX
    milxQtMain *MainWindow;

	// Advanced options window
	milxQtRegistrationAdvancedOptions *advancedOptionsWindow;

	// Create connections with UI
    void createConnections();

	// Is there any work in progress
	bool workInProgress;

	// Do we need to compute the average of the registrations
	bool computeAverage;
	milxQtRegistrationNifti *niftiReg; // nifti reg to compute average
	QString atlasPath; // Path to the outputed atlas
	bool openResults; // do we need to open the results
};

#endif // MILXQTRegistrationWindow_H
