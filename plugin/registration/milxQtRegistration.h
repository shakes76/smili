#ifndef MILXQTRegistration_H
#define MILXQTRegistration_H

#include "milxQtRegistrationStructures.h"
#include "milxQtRegistration.h"
#include "milxQtRegistrationNiftiReg.h"
#include "milxQtMain.h"
#include <QObject>
#include <QFile>


/*!
    \class milxQtRegistration
    \brief Contain all the informations and functions required to register an image
*/
class milxQtRegistration : public QObject
{
    Q_OBJECT

public:

    /*!
        \fn milxQtRegistration::milxQtRegistration(QObject * parent, milxQtImage * imageWindow, milxQtMain * mainW)
        \brief The standard constructor
    */
    milxQtRegistration(QObject * parent, milxQtImage * imageWindow, milxQtMain * mainW);

    /*!
        \fn	milxQtRegistration(QObject * parent, QString, milxQtMain * mainW)
        \brief The standard constructor
    */
    milxQtRegistration(QObject * parent, QString, milxQtMain * mainW);

    /*!
        \fn milxQtRegistration::~milxQtRegistration();
        \brief The standard destructor
    */
    ~milxQtRegistration();

    /*!
        \fn milxQtRegistration::getPath()
        \brief Get the path of the file
    */
    QString getPath();

    /*!
        \fn milxQtRegistration::setChecked(bool)
        \brief Do we have to perform a registration on this image (is the image check in the form)
    */
    void setChecked(bool);

    /*!
        \fn milxQtRegistration::copyPath(char dest[FILENAME_MAX], QString src)
        \brief Copy a path from a QString to a char *
    */
    void copyPath(char dest[FILENAME_MAX], QString src);

    /*!
        \fn milxQtRegistration::isChecked()
        \brief Is the image checked (if we have to perform a registration)
    */
    bool isChecked();

    /*!
        \fn milxQtRegistration::setParams(ParamsF3D)
        \brief Set the F3D registration parameters
    */
    void setParams(ParamsF3D);

    /*!
        \fn milxQtRegistration::setParams(ParamsAladin)
        \brief Set the Aladin registration parameters
    */
    void setParams(ParamsAladin);

    /*!
        \fn milxQtRegistration::setRegType(RegType)
        \brief Set the registration type
    */
    void setRegType(RegType);

    /*!
        \fn milxQtRegistration::setReference(milxQtRegistration *)
        \brief Set the reference image
    */
    void setReference(milxQtRegistration *);

    /*!
        \fn milxQtRegistration::setOutputFolder(QString)
        \brief Set the output folder for the registration
    */
    void setOutputFolder(QString);

    /*!
        \fn milxQtRegistration::milxQtRegistration()
        \brief Is the image open in SMILi
    */
    bool isOpened();

    /*!
        \fn milxQtRegistration::milxQtRegistration()
        \brief Start the registration
    */
    void startRegistration();

    /*!
        \fn milxQtRegistration::createFiles()
        \brief Create the files required for the registration
    */
    void createFiles();

    /*!
        \fn milxQtRegistration::createFile(QString pathtemplate)
        \brief Create a file and return the path
    */
    QString createFile(QString pathtemplate);

    /*!
        \fn milxQtRegistration::deleteTmpFiles()
        \brief Delete the temporary files
    */
    void deleteTmpFiles();

    /*!
        \fn milxQtRegistration::copyAndReplace(QString src, QString dst)
        \brief Copy and replace
    */
    void copyAndReplace(QString src, QString dst);

    /*!
        \fn milxQtRegistration::setOpenResults(bool open)
        \brief Open the results after the registration
    */
    void setOpenResults(bool open);

    /*!
        \fn milxQtRegistration::getOutputPath()
        \brief Return the path of the output file
    */
    QString getOutputPath();

    /*!
        \fn milxQtRegistration::getOutputFolder()
        \brief Return the path of the output folder
    */
    QString getOutputFolder();

    /*!
        \fn milxQtRegistration::createAtlasFile()
        \brief Create the file for the atlas and return the filepath
    */
    QString createAtlasFile();

    /*!
        \fn milxQtRegistration::isWorkDone()
        \brief Is the registration/deformation field calculation completed
    */
    bool isWorkDone();

signals:

    /*!
        \fn milxQtRegistration::done()
        \brief The registration and transformations have been done
    */
    void done();


public slots:

    /*!
        \fn milxQtRegistration::registrationCompleted()
        \brief The registration is completed
    */
    void registrationCompleted();

    /*!
        \fn milxQtRegistration::cpp2defCompleted()
        \brief The transformation to deformation field is completed
    */
    void cpp2defCompleted();

    /*!
        \fn milxQtRegistration::setIsRef(bool)
        \brief Set the image as the reference
    */
    void setIsRef(bool);

    /*!
        \fn milxQtRegistration::isRef()
        \brief Is the image the reference image
    */
    bool isRef();

protected:
    /*!
        \fn milxQtRegistration::init()
        \brief Initialise the class/object with default parameters
    */
    void init(QString);

    milxQtImage * window; //!< MilxQt window of the image, if the image is already opened in SMILI
    bool openedImage; //!<  Is the image already opened in SMILI
    bool openResults; //!< Do we have to open the results after the registration
    QString path; //!< Path of the image
    QString outputFolder; //!< Path of the output folder
    bool checked; //!< Is the image checked: do we need to perform a registration
    bool workDone; //!< Is the registration done
    bool isRefImg; //!< Is this image the reference image
    ParamsF3D paramsF3D; //!< Parameters for a F3D registration
    ParamsAladin paramsAladin; //!< Parameters for an Aladin registration
    RegType type; //!< Type of the registration
    milxQtRegistration * reference; //!< Reference image for the registration
    milxQtRegistrationNifti *niftiReg; //!< Class containing the algorithms for the registration
    milxQtMain *MainWindow; //!< Main window of SMILI

    // Filepath use for the registration
    QString tmp_img; //!< Temp path/file for the image
    QString tmp_ref; //!< Temp path/file for the reference
    QString tmp_cpp; //!< Temp path/file for the cpp file
    QString output_def; //!< Output path for the deformation field
    QString output_path; //!< Output path/file for the image
};

#endif // MILXQTRegistration_H

