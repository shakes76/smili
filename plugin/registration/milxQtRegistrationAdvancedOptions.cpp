#include "milxQtRegistrationAdvancedOptions.h"

#include <QCheckBox>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>

milxQtRegistrationAdvancedOptions::milxQtRegistrationAdvancedOptions(QDialog *theParent) : QDialog(theParent)
{
    ui.setupUi(this);

    setWindowModality(Qt::ApplicationModal);
    setWindowTitle(tr("Registration Advanced Options"));

    createConnections();
    reset(F3DNifti);
}

milxQtRegistrationAdvancedOptions::~milxQtRegistrationAdvancedOptions()
{

}

// Reset the interface
void milxQtRegistrationAdvancedOptions::reset(RegType algo)
{
    QSettings settings("Smili", "Registration Plugin");
    currentAlgo = algo;


    if (algo == F3DNifti)
    {
        ui.titleSplinesOptions->setVisible(true);
        ui.descriptionSplinesOptions->setVisible(true);
        ui.sx->setVisible(true);
        ui.sy->setVisible(true);
        ui.sz->setVisible(true);
        ui.spinBoxSx->setVisible(true);
        ui.spinBoxSy->setVisible(true);
        ui.spinBoxSz->setVisible(true);
        ui.vSpacerSplines->changeSize(20, 20, QSizePolicy::Fixed, QSizePolicy::Fixed);
        ui.vSpacerOptimisations->changeSize(0, 0, QSizePolicy::Fixed, QSizePolicy::Fixed);
        ui.vSpacerSymApproach->changeSize(0, 0, QSizePolicy::Fixed, QSizePolicy::Fixed);

        ui.optimisationsForm->setVerticalSpacing(7);
        ui.lblMaxIt->setVisible(true);
        ui.spinBoxMaxIt->setVisible(true);
        ui.lblLn->setVisible(true);
        ui.spinBoxLn->setVisible(true);
        ui.lblLp->setVisible(true);
        ui.spinBoxLp->setVisible(true);
        ui.checkBoxNopy->setVisible(true);
        ui.lblPctOfBlock->setVisible(false);
        ui.spinBoxPctBlock->setVisible(false);
        ui.checkBoxRigOnly->setVisible(false);
        ui.checkBoxAffDirect->setVisible(false);

        ui.titleSymmetricApproach->setVisible(true);
        ui.checkBoxSym->setVisible(true);

        ui.lblParameterFile->setVisible(false);
        ui.lblParameterFileExplanation->setVisible(false);
        ui.lineParameterFile->setVisible(false);
        ui.btnBrowseParameterFile->setVisible(false);
        ui.lblPctOfBlock->setVisible(false);
        ui.btnClearParameterFile->setVisible(false);

        ui.spinBoxSx->setValue(settings.value("F3D/sx", -5).toFloat());
        ui.spinBoxSy->setValue(settings.value("F3D/sy", -5).toFloat());
        ui.spinBoxSz->setValue(settings.value("F3D/sz", -5).toFloat());

        ui.spinBoxMaxIt->setValue(settings.value("F3D/maxit", 300).toInt());
        ui.spinBoxLn->setValue(settings.value("F3D/ln", 3).toInt());
        ui.spinBoxLp->setValue(settings.value("F3D/lp", 3).toInt());

        ui.checkBoxNopy->setChecked(settings.value("F3D/nopy", false).toBool());
        ui.checkBoxSym->setChecked(settings.value("F3D/sym", false).toBool());
    }

    else if (algo == AladinNifti)
    {
        ui.titleSplinesOptions->setVisible(false);
        ui.descriptionSplinesOptions->setVisible(false);
        ui.sx->setVisible(false);
        ui.sy->setVisible(false);
        ui.sz->setVisible(false);
        ui.spinBoxSx->setVisible(false);
        ui.spinBoxSy->setVisible(false);
        ui.spinBoxSz->setVisible(false);
        ui.checkBoxNopy->setVisible(false);

        ui.vSpacerSplines->changeSize(0, 0, QSizePolicy::Fixed, QSizePolicy::Fixed);
        ui.vSpacerOptimisations->changeSize(20, 20, QSizePolicy::Fixed, QSizePolicy::Fixed);
        ui.vSpacerSymApproach->changeSize(0, 0, QSizePolicy::Fixed, QSizePolicy::Fixed);

        ui.optimisationsForm->setVerticalSpacing(7);
        ui.lblMaxIt->setVisible(true);
        ui.spinBoxMaxIt->setVisible(true);
        ui.lblLn->setVisible(true);
        ui.spinBoxLn->setVisible(true);
        ui.lblLp->setVisible(true);
        ui.spinBoxLp->setVisible(true);
        ui.checkBoxNopy->setVisible(false);
        ui.lblPctOfBlock->setVisible(true);
        ui.spinBoxPctBlock->setVisible(true);
        ui.checkBoxRigOnly->setVisible(true);
        ui.checkBoxAffDirect->setVisible(true);

        ui.titleSymmetricApproach->setVisible(true);
        ui.checkBoxSym->setVisible(true);

        ui.lblParameterFile->setVisible(false);
        ui.lblParameterFileExplanation->setVisible(false);
        ui.lineParameterFile->setVisible(false);
        ui.btnBrowseParameterFile->setVisible(false);
        ui.btnClearParameterFile->setVisible(false);

        ui.spinBoxMaxIt->setValue(settings.value("aladin/maxit", 5).toInt());
        ui.spinBoxLn->setValue(settings.value("aladin/ln", 3).toInt());
        ui.spinBoxLp->setValue(settings.value("aladin/lp", 3).toInt());
        ui.spinBoxPctBlock->setValue(settings.value("aladin/%v", 50).toFloat());

        ui.checkBoxAffDirect->setChecked(settings.value("aladin/aF3Direct", false).toBool());
        ui.checkBoxRigOnly->setChecked(settings.value("aladin/rigOnly", false).toBool());
        ui.checkBoxSym->setChecked(settings.value("aladin/sym", false).toBool());
    }

    else if (algo == ElastixAffine || algo == ElastixBSpline)
    {
        ui.titleSplinesOptions->setVisible(false);
        ui.descriptionSplinesOptions->setVisible(false);
        ui.sx->setVisible(false);
        ui.sy->setVisible(false);
        ui.sz->setVisible(false);
        ui.spinBoxSx->setVisible(false);
        ui.spinBoxSy->setVisible(false);
        ui.spinBoxSz->setVisible(false);
        ui.checkBoxNopy->setVisible(false);

        ui.vSpacerSplines->changeSize(0, 0, QSizePolicy::Fixed, QSizePolicy::Fixed);
        ui.vSpacerOptimisations->changeSize(20, 20, QSizePolicy::Fixed, QSizePolicy::Fixed);
        ui.vSpacerSymApproach->changeSize(0, 0, QSizePolicy::Fixed, QSizePolicy::Fixed);

        ui.optimisationsForm->setVerticalSpacing(0);
        ui.lblMaxIt->setVisible(true);
        ui.spinBoxMaxIt->setVisible(true);
        ui.lblLn->setVisible(false);
        ui.spinBoxLn->setVisible(false);
        ui.lblLp->setVisible(false);
        ui.spinBoxLp->setVisible(false);
        ui.checkBoxNopy->setVisible(false);
        ui.lblPctOfBlock->setVisible(false);
        ui.spinBoxPctBlock->setVisible(false);
        ui.checkBoxRigOnly->setVisible(false);
        ui.checkBoxAffDirect->setVisible(false);

        ui.titleSymmetricApproach->setVisible(false);
        ui.checkBoxSym->setVisible(false);

        ui.lblParameterFile->setVisible(true);
        ui.lineParameterFile->setVisible(true);
        ui.lblParameterFileExplanation->setVisible(true);
        ui.btnBrowseParameterFile->setVisible(true);
        ui.btnClearParameterFile->setVisible(true);

        if (algo == ElastixAffine)
        {
            ui.spinBoxMaxIt->setValue(settings.value("elastixAffine/maxit", 250).toInt());
            ui.lineParameterFile->setText(settings.value("elastixAffine/parameterFile", "").toString());
        }
        else if (algo == ElastixBSpline)
        {
            ui.spinBoxMaxIt->setValue(settings.value("elastixBSpline/maxit", 500).toInt());
            ui.lineParameterFile->setText(settings.value("elastixBSpline/parameterFile", "").toString());
        }
    }

    this->resize(10, 10);
}

// Create connections with the user interface
void milxQtRegistrationAdvancedOptions::createConnections()
{
    connect(this->ui.btnBrowseParameterFile, SIGNAL(clicked()), this, SLOT(browseBtnClicked()));
    connect(this->ui.btnClearParameterFile, SIGNAL(clicked()), this, SLOT(clearFileBtnClicked()));
}

void milxQtRegistrationAdvancedOptions::browseBtnClicked()
{
    QFileDialog fileOpener;
    fileOpener.setFileMode(QFileDialog::FileMode::ExistingFile);
    QString filename = fileOpener.getOpenFileName(this, tr("Select File"));
    this->ui.lineParameterFile->setText(filename);
}

void milxQtRegistrationAdvancedOptions::clearFileBtnClicked()
{
    this->ui.lineParameterFile->setText("");
}


milxQtRegistrationParams milxQtRegistrationAdvancedOptions::getParamsElastixAffine()
{
    milxQtRegistrationParams params;

    if (this->ui.lineParameterFile->text() != "" && QFile::exists(this->ui.lineParameterFile->text()))
    {
        params.parameterFile = this->ui.lineParameterFile->text();
        params.customParameterFile = true;

        return params;
    }

    // Default parameters
    params.parametersTxt =
        "(FixedInternalImagePixelType \"float\")\r\n"
        "(MovingInternalImagePixelType \"float\")\r\n"
        "(UseDirectionCosines \"true\")\r\n"
        "(Registration \"MultiResolutionRegistration\")\r\n"
        "(Interpolator \"BSplineInterpolator\")\r\n"
        "(ResampleInterpolator \"FinalBSplineInterpolator\")\r\n"
        "(Resampler \"DefaultResampler\")\r\n"
        "(FixedImagePyramid \"FixedRecursiveImagePyramid\")\r\n"
        "(MovingImagePyramid \"MovingRecursiveImagePyramid\")\r\n"
        "(Optimizer \"AdaptiveStochasticGradientDescent\")\r\n"
        "(Transform \"AffineTransform\")\r\n"
        "(Metric \"AdvancedMattesMutualInformation\")\r\n"
        "(AutomaticScalesEstimation \"true\")\r\n"
        "(AutomaticTransformInitialization \"true\")\r\n"
        "(HowToCombineTransforms \"Compose\")\r\n"
        "(NumberOfHistogramBins 32)\r\n"
        "(ErodeMask \"false\")\r\n"
        "(NumberOfResolutions 4)\r\n"
        "(MaximumNumberOfIterations " + QString::number(ui.spinBoxMaxIt->value()) + ")\r\n"
        "(NumberOfSpatialSamples 2048)\r\n"
        "(NewSamplesEveryIteration \"true\")\r\n"
        "(ImageSampler \"Random\")\r\n"
        "(BSplineInterpolationOrder 1)\r\n"
        "(FinalBSplineInterpolationOrder 3)\r\n"
        "(DefaultPixelValue 0)\r\n"
        "(WriteResultImage \"true\")\r\n"
        "(ResultImageFormat \"nii\")\r\n";

    return params;
}


milxQtRegistrationParams milxQtRegistrationAdvancedOptions::getParamsElastixBSpline()
{
    milxQtRegistrationParams params;

    if (this->ui.lineParameterFile->text() != "" && QFile::exists(this->ui.lineParameterFile->text()))
    {
        params.parameterFile = this->ui.lineParameterFile->text();
        params.customParameterFile = true;

        return params;
    }

    // Default parameters
    params.parametersTxt =
        "(FixedInternalImagePixelType \"float\")\r\n"
        "(MovingInternalImagePixelType \"float\")\r\n"
        "(UseDirectionCosines \"true\")\r\n"
        "(Registration \"MultiResolutionRegistration\")\r\n"
        "(Interpolator \"BSplineInterpolator\")\r\n"
        "(ResampleInterpolator \"FinalBSplineInterpolator\")\r\n"
        "(Resampler \"DefaultResampler\")\r\n"
        "(FixedImagePyramid \"FixedRecursiveImagePyramid\")\r\n"
        "(MovingImagePyramid \"MovingRecursiveImagePyramid\")\r\n"
        "(Optimizer \"AdaptiveStochasticGradientDescent\")\r\n"
        "(Transform \"BSplineTransform\")\r\n"
        "(Metric \"AdvancedMattesMutualInformation\")\r\n"
        "(FinalGridSpacingInPhysicalUnits 16)\r\n"
        "(HowToCombineTransforms \"Compose\")\r\n"
        "(NumberOfHistogramBins 32)\r\n"
        "(ErodeMask \"false\")\r\n"
        "(NumberOfResolutions 4)\r\n"
        "(MaximumNumberOfIterations " + QString::number(ui.spinBoxMaxIt->value()) + ")\r\n"
        "(NumberOfSpatialSamples 2048)\r\n"
        "(NewSamplesEveryIteration \"true\")\r\n"
        "(ImageSampler \"Random\")\r\n"
        "(BSplineInterpolationOrder 1)\r\n"
        "(FinalBSplineInterpolationOrder 3)\r\n"
        "(DefaultPixelValue 0)\r\n"
        "(WriteResultImage \"true\")\r\n"
        "(ResultImageFormat \"nii\")\r\n";


    return params;
}


milxQtRegistrationParams milxQtRegistrationAdvancedOptions::getParamsF3DNifti()
{
    milxQtRegistrationParams params;

    if (currentAlgo != F3DNifti)
    {
        this->reset(F3DNifti);
    }

    // Get input Sx
    params.spacing[0] = ui.spinBoxSx->value();

    // Get input Sy
    params.spacing[1] = ui.spinBoxSy->value();

    // Get input Sz
    params.spacing[2] = ui.spinBoxSz->value();

    // Get input max it level
    params.maxit = ui.spinBoxMaxIt->value();

    // Get input nb level
    params.ln = ui.spinBoxLn->value();

    // Get input first levels
    params.lp = ui.spinBoxLp->value();

    // Get input use pyramidal approach
    params.nopy = ui.checkBoxNopy->isChecked();

    // Get input ise symmetric approach
    params.useSym = ui.checkBoxSym->isChecked();

    return params;
}


milxQtRegistrationParams milxQtRegistrationAdvancedOptions::getParamsAladinNifti()
{
    milxQtRegistrationParams params;

    if (currentAlgo != AladinNifti)
    {
        this->reset(AladinNifti);
    }

    // Get input max it level
    params.maxit = ui.spinBoxMaxIt->value();

    // Get input nb level
    params.ln = ui.spinBoxLn->value();

    // Get input first levels
    params.lp = ui.spinBoxLp->value();

    // Get percentage of block
    params.percentBlock = ui.spinBoxPctBlock->value();

    // Aladin direct
    params.affineDirect = ui.checkBoxAffDirect->isChecked();

    // Rigid only
    params.rigOnly = ui.checkBoxRigOnly->isChecked();

    // Get input ise symmetric approach
    params.useSym = ui.checkBoxSym->isChecked();

    return params;
}


// Btn Ok clicked
void milxQtRegistrationAdvancedOptions::accept()
{
    // Save parameters
    QSettings settings("Smili", "Registration Plugin");

    if (currentAlgo == F3DNifti)
    {
        settings.setValue("F3D/sx", ui.spinBoxSx->value());
        settings.setValue("F3D/sy", ui.spinBoxSy->value());
        settings.setValue("F3D/sz", ui.spinBoxSz->value());

        settings.setValue("F3D/maxit", ui.spinBoxMaxIt->value());
        settings.setValue("F3D/ln", ui.spinBoxLn->value());
        settings.setValue("F3D/lp", ui.spinBoxLp->value());

        settings.setValue("F3D/nopy", ui.checkBoxNopy->isChecked());
        settings.setValue("F3D/sym", ui.checkBoxSym->isChecked());
    }
    else if (currentAlgo == AladinNifti)
    {
        settings.setValue("aladin/maxit", ui.spinBoxMaxIt->value());
        settings.setValue("aladin/ln", ui.spinBoxLn->value());
        settings.setValue("aladin/lp", ui.spinBoxLp->value());
        settings.setValue("aladin/%v", ui.spinBoxPctBlock->value());

        settings.setValue("aladin/aF3Direct", ui.checkBoxAffDirect->isChecked());
        settings.setValue("aladin/rigOnly", ui.checkBoxRigOnly->isChecked());
        settings.setValue("aladin/sym", ui.checkBoxSym->isChecked());
    }
    else if (currentAlgo == ElastixAffine)
    {
        settings.setValue("elastixAffine/maxit", ui.spinBoxMaxIt->value());
        settings.setValue("elastixAffine/parameterFile", ui.lineParameterFile->text());
    }
    else if (currentAlgo == ElastixBSpline)
    {
        settings.setValue("elastixBSpline/maxit", 500);
        settings.setValue("elastixBSpline/parameterFile", ui.lineParameterFile->text());
    }

    // Close the window
    this->close();
    this->hide();
}

// Btn Cancel clicked
void milxQtRegistrationAdvancedOptions::reject()
{
    // Close the window
    this->close();
    this->hide();
}
