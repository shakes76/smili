#ifndef MILXQTRegistrationImage_H
#define MILXQTRegistrationImage_H

#include "milxQtRegistrationStructures.h"
#include "milxQtRegistrationImage.h"
#include "milxQtRegistrationNiftiReg.h"
#include "milxQtMain.h"
#include <QObject>
#include <QFile>


/*!
    \class milxQtRegistrationImage
    \brief Contain all the information required for a registration
*/
class milxQtRegistrationImage : public QObject
{
	Q_OBJECT

public:
	milxQtRegistrationImage(QObject * parent, milxQtImage * imageWindow, milxQtMain * mainW);
	milxQtRegistrationImage(QObject * parent, QString, milxQtMain * mainW);
	~milxQtRegistrationImage();

	QString getPath(); // Path of the file
	void setChecked(bool); // Do we have to perform a registration on this image
	void copyPath(char dest[FILENAME_MAX], QString src); // Copy a path from a QString to a char *
	bool isChecked(); // is the image checked (if we have to perform a registration)
	void setParams(ParamsF3D); // set F3D registration parameters
	void setParams(ParamsAladin); // set Aladin registration parameters
	void setRegType(RegType); // set the registration type	
	void setReference(milxQtRegistrationImage *); // set the reference image
	void setOutputFolder(QString); // Set the output folder
	bool isOpened(); // Is the image opened
	void startRegistration(); // Start the registration
	void createFiles(); // Create the files required for the registration
	QString createFile(QString pathtemplate); // Create a file and return the path
	void deleteTmpFiles(); // Delete the temporary files
	void copyAndReplace(QString src, QString dst); // Copy and replace
	void setOpenResults(bool open); // Open the results after the registration
	QString getOutputPath(); // Path of the output file
	QString getOutputFolder(); // Return the output folder
	QString createAtlasFile(); // Create the file for the Atlas
	bool isWorkDone(); // Is the registration/deformation field calculation completed

signals:
	void done(); // The registration and transformations have been done


public slots:
	void registrationCompleted(); // the registration is completed
	void cpp2defCompleted(); // the transformation to deformation field is completed
	void setIsRef(bool); // Set the image as the reference
	bool isRef(); // Is the image the reference image

protected:
	void init(QString);
	milxQtImage * window;
	bool openedImage;
	bool openResults;
	QString path;
	QString outputFolder;
	bool checked;
	bool workDone;
	bool isRefImg;
	ParamsF3D paramsF3D;
	ParamsAladin paramsAladin;
	RegType type;
	milxQtRegistrationImage * reference;
	milxQtRegistrationNifti *niftiReg;
	milxQtMain *MainWindow;

	// Filepath used for the registration
	QString tmp_img;
	QString tmp_ref;
	QString tmp_cpp;
	QString output_def;
	QString output_path;
};

#endif // MILXQTRegistrationImage_H

