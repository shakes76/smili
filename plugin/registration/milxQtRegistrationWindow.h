#ifndef MILXQTRegistrationWindow_H
#define MILXQTRegistrationWindow_H

#include "ui_registrationWindow.h"
#include "milxQtRegistrationAdvancedOptions.h"
#include "milxQtMain.h"
#include "milxQtImage.h"
#include "milxQtFile.h"
#include "milxQtRegistrationNiftiReg.h"
#include "milxQtRegistration.h"

typedef QList<milxQtImage *> QImageList;

/*!
    \class milxQtRegistrationWindow
    \brief This class is the registration window/form
*/
class milxQtRegistrationWindow : public QDialog
{
    Q_OBJECT

public:
    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief The standard constructor
    */
    milxQtRegistrationWindow(QWidget * theParent);

    /*!
        \fn milxQtRegistrationWindow::~milxQtRegistrationWindow(QWidget *parent = 0)
        \brief The standard destructor
    */
    virtual ~milxQtRegistrationWindow();

    /*!
        \fn milxQtRegistrationWindow::initUI()
        \brief The standard constructor
    */
    void initUI();

    /*!
        \fn milxQtRegistrationWindow::setAlgo()
        \brief Set the type of algorithm and change the form accordingly (display the correct fields)
    */
    void setAlgo(RegType regType);

    /*!
        \fn milxQtRegistrationWindow::updateImageListCombo()
        \brief Update the list of images in the combo box
    */
    void updateImageListCombo();

    /*!
        \fn milxQtRegistrationWindow::updateOpenImages()
        \brief Update the list of open images in SMILI
    */
    void updateOpenImages();

    /*!
        \fn milxQtRegistrationWindow::isImageInList()
        \brief Is the image already listed in our list of images
    */
    bool isImageInList(QString path);

    /*!
        \fn milxQtRegistrationWindow::updateParameters()
        \brief Update images parameters, set the values of the form to images
    */
    void updateParameters();

    /*!
        \fn milxQtRegistrationWindow::getParamsF3D()
        \brief Return the parameters for a F3D registration
    */
    ParamsF3D getParamsF3D();

    /*!
        \fn milxQtRegistrationWindow::getParamsAladin()
        \brief Return the parameters for an Aladin registration
    */
    ParamsAladin getParamsAladin();

    /*!
        \fn milxQtRegistrationWindow::addImage(milxQtRegistration *)
        \brief Add an image to our list of images
    */
    void addImage(milxQtRegistration *);

    /*!
        \fn milxQtRegistrationWindow::disableUI()
        \brief Disable the user interface
    */
    void disableUI();

    /*!
        \fn milxQtRegistrationWindow::enableUI()
        \brief Enable the user interface
    */
    void enableUI();

    /*!
        \fn milxQtRegistrationWindow::performRegistrations()
        \brief Perform the next registration
    */
    void performRegistrations();

    /*!
        \fn milxQtRegistrationWindow::workCompleted()
        \brief Everything has been completed
    */
    void workCompleted();

    /*!
        \fn milxQtRegistrationWindow::computeAtlas()
        \brief Compute the average (Atlas) of all the registrations
    */
    void computeAtlas();

    QList<milxQtRegistration *> images; //!< List of images for the combobox and list of selectable images

public slots:
    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Click on button Ok
    */
    void accept();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Click on button Cancel
    */
    void reject();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Reference combo box changed
    */
    void referenceComboChange(int newIndex);

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Algo combo box changed
    */
    void algoComboChange(int newIndex);

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Open the advanced option windows
    */
    void advancedOptionsClicked();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Button add image clicked
    */
    void addImageClicked();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Button select all clicked
    */
    void selectAllClicked();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Button unselect all clicked
    */
    void unselectAllClicked();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Button browse clicked
    */
    void browseBtnClicked();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief A registration has been completed
    */
    void regComplete();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Clear the list of images
    */
    void clearList();

    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief The average (Atlas) has been computed
    */
    void averageComputed();

protected:
    /*!
        \fn milxQtRegistrationWindow::milxQtRegistrationWindow(QWidget *parent = 0)
        \brief Create connections with UI
    */
    void createConnections();

    Ui::dlgRegistrationWindow ui; //!< User interface
    milxQtMain *MainWindow; //!< Main window of SMILIX
    milxQtRegistrationAdvancedOptions *advancedOptionsWindow; //!< Advanced options window
    bool workInProgress; //!< Is there any work in progress
    bool computeAverage; //!< Do we need to compute the average of the registrations
    milxQtRegistrationNifti *niftiReg; //!< nifti reg to compute average
    QString atlasPath; //!< Path to the outputed atlas
    bool openResults; //!< Store if we need to open the results
};

#endif // MILXQTRegistrationWindow_H
