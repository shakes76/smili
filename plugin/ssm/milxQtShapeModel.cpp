/*=========================================================================
  The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO)
  ABN 41 687 119 230.
  All rights reserved.

  Licensed under the CSIRO BSD 3-Clause License
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License in the file LICENSE.md or at

  https://stash.csiro.au/projects/SMILI/repos/smili/browse/license.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
=========================================================================*/
#include <ctime>

#include "milxQtShapeModel.h"

//Used for generateModesModel
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkTensor.h>
#include <vtkTensorGlyph.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkMath.h>
#include <vtkSphericalTransform.h>
#include <vnl/vnl_double_3.h>
#include <vnl/algo/vnl_scatter_3x3.h>
//Improve Rendering
#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>

#include <milxQtFile.h>
#include <milxQtPlot.h>

milxQtShapeModel::milxQtShapeModel(QWidget *theParent, bool contextSystem) : milxQtRenderWindow(theParent, contextSystem)
{
    reset();

    m_StandardSSM = ShapeModelType::New();
    m_SSM = m_StandardSSM;

    m_meanModel = NULL;
    m_modesVectorModel = NULL;
    m_modesTensorModel = NULL;
    m_correspondences = NULL;

    ///Set strings
    milxQtWindow::prefix = "SSM: ";

    createActions();
    setupTooltips();
    createConnections();
}

milxQtShapeModel::~milxQtShapeModel()
{

}

void milxQtShapeModel::LoadModel(const QString filename)
{
    reset();

    if(m_SSM->LoadModel(filename.toStdString().c_str()))
        m_loaded = true;
    else
    {
        printError("Could not open Robust SSM file.");
        return;
    }

    generateSSM();
}

void milxQtShapeModel::SetInputCollection(vtkPolyDataCollection* meshes)
{
    size_t n = meshes->GetNumberOfItems();

    working(-1);
    reset();

    meshes->InitTraversal();
    for(size_t j = 0; j < n; j ++)
    {
        m_SSM->AddShape(meshes->GetNextItem()); //!< Set model to shape filter

        qApp->processEvents(); ///Keep UI responsive
    }
    done(-1);

    m_loaded = true;
}

void milxQtShapeModel::SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames)
{
    const int n = meshes->GetNumberOfItems();

    working(-1);
    reset();

    ///Setup temporary combo widget for selecting case ids when there are
    ///Multiple integer matches in the filename
    int index = 0;
    QDialog *casePossibilities = new QDialog(this);
    QComboBox *integerValues = new QComboBox(this);
    QPushButton *okButton = new QPushButton(this);
        okButton->setText("Ok");
    QLabel *lblMessage = new QLabel(this);
        lblMessage->setText("Choose Case ID from possibilities.");
    QFormLayout *formLayout = new QFormLayout(this);
        formLayout->addRow(lblMessage);
        formLayout->addRow(tr("&Select Case ID from first file: "), integerValues);
        formLayout->addRow(okButton);
        casePossibilities->setLayout(formLayout);

    connect(okButton, SIGNAL(clicked(bool)), casePossibilities, SLOT(accept()));

    meshes->InitTraversal();
    for(int j = 0; j < n; j ++)
    {
        QFileInfo fi(filenames[j]);
        QRegExp rx("(\\d+)", Qt::CaseSensitive, QRegExp::RegExp2); ///Capture integer expression

        ///Read all integers from filename and save in list
        int pos = 0;
        QStringList list;
        while ((pos = rx.indexIn(fi.baseName(), pos)) != -1)
        {
            list << rx.cap(1);
            pos += rx.matchedLength(); //Move
        }

        ///If more than one integer in filename, ask user to select correct integer
        if(list.size() > 1 && j == 0)
        {
            printInfo("Please choose Case ID from possibilities for first file.");
            for(int k = 0; k < list.size(); k ++)
                integerValues->addItem(list[k]);

            casePossibilities->exec();

            index = integerValues->currentIndex();
        }

        int caseID = 0;
        if(!list.isEmpty())
            caseID = list[index].toInt(); //!< Extract case ID from filename by using integer user provided (or only one present)

        m_caseIDs.append(caseID); //!< Add case ID to list

        m_SSM->AddShape(meshes->GetNextItem()); //!< Set model to shape filter

        qApp->processEvents(); ///Keep UI responsive
    }
    done(-1);

    qDebug() << "Case IDs: " << m_caseIDs << endl;

    m_loaded = true;
}

QPointer<milxQtModel> milxQtShapeModel::getMeanModel()
{
    if(!m_loaded)
        return NULL;
    if(!m_modelled)
        generateSSM();
    if(!m_meaned)
    {
        generateMeanModel(); ///Displays Mean mesh so
        RemoveActor(m_meanModel->GetActor()); ///Dont display mean mesh
    }

    return m_meanModel;
}

vtkDataSet* milxQtShapeModel::GetDataSet()
{
    if(m_modesVectorModel)
        return m_modesVectorModel->GetOutput();
    else if(m_modesTensorModel)
        return m_modesTensorModel->GetOutput();
    else if(m_correspondences)
        return m_correspondences->GetOutput();
    else
        return m_meanModel->GetOutput();
}

void milxQtShapeModel::createMenu(QMenu *menu)
{
    if(!menu)
        return;

    menu->clear();

    menu->addSeparator()->setText(tr("View"));
    menu->addAction(actionMean);
    menu->addAction(actionAligned);
    menu->addAction(actionOriginal);
    menu->addAction(actionModesAsVectors);
    menu->addAction(actionModesAsTensors);
    menu->addAction(actionModesAsCollection);
    menu->addAction(actionCorrespond);
    plotMenu = menu->addMenu("Plot");
    plotMenu->addAction(actionCompact);
    plotMenu->addAction(actionSpecificity);
    plotMenu->addAction(actionGeneralise);
    plotMenu->addAction(actionValues);
    plotMenu->addAction(actionModes);
    plotMenu->addAction(actionParameters);
    menu->addSeparator()->setText(tr("Validation/Output"));
    menu->addAction(actionOriginalMeshes);
    menu->addAction(actionAlignedMeshes);
    menu->addAction(actionAlignment);
    menu->addAction(actionPointIDs);
    menu->addAction(actionCoordinates);
    menu->addSeparator();
    menu->addAction(actionReplaceOriginal);
    //Align
    alignMenu = menu->addMenu("Procrustes Mode");
    alignMenu->addAction(actionSimilarity);
    alignMenu->addAction(actionRigid);
    alignMenu->addAction(actionAffine);
    ///Change View of Volume
    menu->addMenu(milxQtRenderWindow::viewMenu);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewXY);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZX);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZY);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewFileAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewFileAct);
    milxQtRenderWindow::enableActionBasedOnView();
    menu->addMenu(milxQtRenderWindow::colourMapMenu);
    menu->addSeparator();
    menu->addMenu(milxQtRenderWindow::windowPropertiesMenu);
    menu->addAction(milxQtRenderWindow::refreshAct);
    menu->addAction(milxQtRenderWindow::resetAct);
}

void milxQtShapeModel::generateSSM()
{
    if(!m_loaded)
        return;

    bool ok;
    float precision = QInputDialog::getDouble(this, tr("Please Provide the precision of the model"),
                              tr("Precision:"), 0.9, 0.0, 1.0, 5, &ok);

    QMessageBox msgBox;
    msgBox.setText("The mean surface/shape will be used in the model");
    msgBox.setInformativeText("Do you want to use the first surface instead of the mean?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();

    if(ret == QMessageBox::Yes)
    {
        m_StandardSSM->TimepointModeOn();
    }

    if(!ok) //cancelled
        return;

    ///Compute SSM
    try
    {
        m_StandardSSM->SetPrecision(precision);
        m_StandardSSM->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        printError("Exception caught while generating an SSM!");
        cout << err << endl;
        return;
    }

    const size_t n = m_StandardSSM->GetNumberOfShapes();

    ///Compute the sum of eigenvalues, which is a measure of the variation of the modelling
    double eigenSum = 0.0;
    vtkSmartPointer<vtkFloatArray> eigenVals = m_StandardSSM->GetPCA()->GetEvals();
        printInfo("Eigenvalues Info - #tuples: " + QString::number(eigenVals->GetNumberOfTuples()) + ", #components: " + QString::number(eigenVals->GetNumberOfComponents()));
        for(vtkIdType j = 0; j < eigenVals->GetNumberOfTuples(); j ++)
            eigenSum += eigenVals->GetValue(j);
        printInfo("Sum of Eigenvalues is " + QString::number(eigenSum));
        printInfo("Mean Scale: " + QString::number(m_StandardSSM->GetMeanShapeSize()));
        printInfo("Number of Modes: " + QString::number(m_StandardSSM->GetNumberOfModes()));

    vtkSmartPointer<vtkLookupTable> tmpLookupTable = vtkLookupTable::SafeDownCast(milxQtRenderWindow::lookupTable);
        tmpLookupTable->SetTableRange(0.0, n+1); ///Build colour lookup table
        tmpLookupTable->Build();
    milxQtRenderWindow::lookupTable = tmpLookupTable;

    milxQtRenderWindow::setView(milxQtRenderWindow::defaultView); //Default view

    m_modelled = true;
}

void milxQtShapeModel::generateMeanModel(vtkSmartPointer<vtkPolyData> shape)
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();

    printInfo("Generating Mean Model.");
    if(!m_meanModel)
    {
        m_meanModel = new milxQtModel;
        m_meanModel->SetInput(m_SSM->GetMeanShape());
    }
    if(!shape)
    {
        m_meanModel->setName("Mean");
        m_meanModel->SetPoints(m_SSM->GetMeanShape()->GetPoints());
    }
    else
    {
        m_meanModel->setName("Custom Shape");
        printDebug("Adding Custom Model.");
        m_meanModel->SetInput(shape);
    }
    m_meanModel->generateModel();

    AddActor(m_meanModel->GetActor()); //!< Add data to general display

    m_meaned = true;
    actionMean->setChecked(true);
}

void milxQtShapeModel::generateModels(const bool display)
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();

    const int n = m_SSM->GetNumberOfShapes();

    working(-1);
    printInfo("Generating Original Models.");
    for(int j = 0; j < n; j ++)
    {
        double colour[3];

        milxQtRenderWindow::lookupTable->GetColor(j, colour); //!< Pull colour for data

        ///Build model
        QPointer<milxQtModel> shape = new milxQtModel;
            shape->setName( QString::number(j) );
            shape->setLargeDataSetMode(true);
            shape->SetInput(m_SSM->GetShape(j));
            shape->generatePoints(colour[0], colour[1], colour[2]);
            shape->SetOpacity(0.1);

        m_models.append(shape);

        if(display)
            AddActor(m_models.last()->GetActor()); //!< Add data to general display

        qApp->processEvents(); ///Keep UI responsive
        milxQtRenderWindow::Render();
    }
    printInfo("Done. Displaying Original Models.");

    if(display)
    {
        milxQtRenderWindow::generateRender();
        actionOriginal->setChecked(true);
    }
    m_shapesModelled = true;
    done(-1);
}

void milxQtShapeModel::generateAlignedModels(const bool display)
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();

    emit working(-1);
    const int n = m_SSM->GetNumberOfShapes();

    printInfo("Generating Aligned Models.");
    for(int j = 0; j < n; j ++)
    {
        double colour[3];

        milxQtRenderWindow::lookupTable->GetColor(j, colour); //!< Pull colour for data

        ///Build model
        QPointer<milxQtModel> shape = new milxQtModel;
            shape->setName( QString::number(j) );
            shape->setLargeDataSetMode(true);
            //printDebug("Get Aligned Shape " + QString::number(j));
            shape->SetInput(m_SSM->GetProcrustesAlignedSurface(j)); //Crashes sometimes?
//            shape->SetInput(m_SSM->GetMeanShape());
//            shape->SetPoints(m_SSM->GetAlignedPoints(j));
            shape->generatePoints(colour[0], colour[1], colour[2]);
            shape->SetOpacity(0.1);
            shape->Update();

        //printDebug("Add to Aligned Models Set");
        m_alignedModels.append(shape);

        if(display)
            AddActor(m_alignedModels.last()->GetActor()); //!< Add data to general display

        qApp->processEvents(); ///Keep UI responsive
        milxQtRenderWindow::Render();
    }
    printInfo("Done. Displaying Aligned Models.");

    if(display)
    {
        milxQtRenderWindow::generateRender();
        actionAligned->setChecked(true);
    }
    m_alignedShapesModelled = true;
    done(-1);
}

void milxQtShapeModel::generateModes()
{
    bool flgNormal = true;

    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();
    if(!m_meaned)
        generateMeanModel();

    bool ok1, ok2;
    int mode = QInputDialog::getInt(this, tr("Please Provide the mode to view (1 to n)"),
                              tr("Mode:"), 1, 1, m_StandardSSM->GetNumberOfModes(), 1, &ok1);
    float modeWeight = QInputDialog::getDouble(this, tr("Please Provide the weight of the mode"),
                              tr("Weight:"), 3.0, -50.0, 50.0, 2, &ok2);

    m_mode = mode;
    m_vector = false;
    m_tensor = false;
    m_modes = false;

    if(!ok1 || !ok2) //cancelled
        return;

    QMessageBox msgBox;
    msgBox.setText("A parameterised surface will be generated");
    msgBox.setInformativeText("Do you want to display this instead of the mean surface?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();

    /**
    Code contributed by Kaikai Shen, 2010.
    Smart Pointered and Modified by Shekhar Chandra, 2010.
    */
    working(-1);
    vtkSmartPointer<vtkPolyData> meanShape = m_meanModel->GetOutput();
    //Ensure mean is actually mean
    meanShape->SetPoints(m_StandardSSM->GetMeanShape()->GetPoints());

    vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
        normals->SetInput(meanShape);
        normals->ComputeCellNormalsOff();
        normals->ComputePointNormalsOn();
        normals->SplittingOff();
        // normals->AutoOrientNormalsOn();
        // normals->ConsistencyOn();
        // normals->FlipNormalsOn();
        normals->Update();

    vtkSmartPointer<vtkFloatArray> normArray =
        vtkFloatArray::SafeDownCast(normals->GetOutput()->GetPointData()->GetNormals());

    vtkSmartPointer<vtkFloatArray> b = vtkFloatArray::New();
        b->SetNumberOfValues(mode);
        b->FillComponent(0, 0.0);
        b->SetValue(mode-1, modeWeight); //!< View mode
        printInfo("Generating display for mode " + QString::number(mode-1));

    vtkSmartPointer<vtkPolyData> varShape = vtkSmartPointer<vtkPolyData>::New();
        varShape->DeepCopy(m_StandardSSM->GetParameterisedShape(b));

    vtkSmartPointer<vtkFloatArray> varScalars = vtkSmartPointer<vtkFloatArray>::New();
        varScalars->SetName("Variation");
    vtkSmartPointer<vtkFloatArray> varVectors = vtkSmartPointer<vtkFloatArray>::New();
        varVectors->SetNumberOfComponents(3);
    vtkSmartPointer<vtkFloatArray> varShapeVectors = vtkSmartPointer<vtkFloatArray>::New();
        varShapeVectors->SetNumberOfComponents(3);
    vtkSmartPointer<vtkFloatArray> varTensors = vtkSmartPointer<vtkFloatArray>::New();
        varTensors->SetNumberOfComponents(9);

    printInfo("Computing Variation");
    for(int i = 0; i < meanShape->GetNumberOfPoints(); i++)
    {
        vtkFloatingPointType* meanPoint = meanShape->GetPoint(i);
        vtkFloatingPointType* varPoint  = varShape->GetPoint(i);

        vtkFloatingPointType xVal = varPoint[0] - meanPoint[0];
        vtkFloatingPointType yVal = varPoint[1] - meanPoint[1];
        vtkFloatingPointType zVal = varPoint[2] - meanPoint[2];

        float normalPoint[3];
        normArray->GetTupleValue(i, normalPoint);
        vtkFloatingPointType var;
        if(flgNormal)
        {
            var = xVal*normalPoint[0] + yVal*normalPoint[1] + zVal*normalPoint[2];
            var = var*var;
        }
        else
        {
            var = xVal*xVal + yVal*yVal + zVal*zVal;
        }
        // std::cout << var << std::endl;
        varVectors->InsertNextTuple3(xVal, yVal, zVal);
        varShapeVectors->InsertNextTuple3(-xVal, -yVal, -zVal);
        varScalars->InsertNextValue(var);

        vnl_scatter_3x3<vtkFloatingPointType> C;
        for(int j = 0; j < m_StandardSSM->GetNumberOfShapes(); j++)
        {
            vtkFloatingPointType point[3];
            m_StandardSSM->GetPCA()->GetInput(j)->GetPoint(i, point);
            // std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
            vnl_double_3 point_vector;
            point_vector[0] = point[0] - meanPoint[0];
            point_vector[1] = point[1] - meanPoint[1];
            point_vector[2] = point[2] - meanPoint[2];
            // std::cout << point_vector << std::endl;
            C.add_outer_product(point_vector);
        }
        // C /= m_StandardSSM->GetNumberOfShapes();
        // std::cout << meanPoint[0] << " " << meanPoint[1] << " " << meanPoint[2] << std::endl << std::endl;
        // std::cout << C << std::endl;
        vtkSmartPointer<vtkTensor> tens = vtkSmartPointer<vtkTensor>::New();
        tens->SetComponent(0,0, C(0,0));
        tens->SetComponent(0,1, C(0,1));
        tens->SetComponent(0,2, C(0,2));
        tens->SetComponent(1,0, C(1,0));
        tens->SetComponent(1,1, C(1,1));
        tens->SetComponent(1,2, C(1,2));
        tens->SetComponent(2,0, C(2,0));
        tens->SetComponent(2,1, C(2,1));
        tens->SetComponent(2,2, C(2,2));

        varTensors->InsertNextTuple(tens->T);
        // std::cout << tens->T[0] << " " << tens->T[1] << " " << tens->T[2] << std::endl;

        qApp->processEvents(); ///Keep UI responsive
    }

    vtkFloatingPointType maxScalar = std::numeric_limits<vtkFloatingPointType>::min();
    for(int i = 0; i < varScalars->GetNumberOfTuples(); i++)
    {
        if(varScalars->GetValue(i) > maxScalar)
            maxScalar = varScalars->GetValue(i);
    }
    printInfo("Max sigma: " + QString::number(sqrt(maxScalar)) + "mm.");
    for(int i = 0; i < varScalars->GetNumberOfTuples(); i++)
    {
        vtkFloatingPointType var = varScalars->GetValue(i);
        varScalars->SetValue(i, var/maxScalar);
    }

    meanShape->GetPointData()->SetVectors(varVectors);
    varShape->GetPointData()->SetVectors(varShapeVectors);
    meanShape->GetPointData()->SetScalars(varScalars);
    //scalars dont make sense for varShape
    meanShape->GetPointData()->SetTensors(varTensors);
    varShape->GetPointData()->SetTensors(varTensors);

    if(ret == QMessageBox::Yes)
        generateMeanModel(varShape);
    else
        generateMeanModel();
    done(-1);

    m_modes = true;
}

void milxQtShapeModel::generateCollectionBasedOnMode()
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();
    if(!m_meaned)
        generateMeanModel();

    const int n = m_SSM->GetNumberOfShapes();
    bool ok1, ok2, ok3;
    int mode = QInputDialog::getInt(this, tr("Please Provide the mode to view (1 to n)"),
                              tr("Mode:"), 1, 1, n, 1, &ok1);
    float modeWeight = QInputDialog::getDouble(this, tr("Please Provide the range weight c of the mode (so that -c < x < c)"),
                              tr("Weight:"), 3.0, -5.0, 5.0, 2, &ok2);
    int numberOfSurfaces = QInputDialog::getInt(this, tr("Please Provide the number of surfaces to generate"),
                              tr("Number:"), 16, 1, 1e9, 1, &ok3);

    m_mode = mode;

    if(!ok1 || !ok2 || !ok3) //cancelled
        return;

    //Init counters etc.
    float weightValue = -modeWeight;
    const float weightStep = 2.0*modeWeight/(numberOfSurfaces-1);

    vtkPolyDataCollection* meshes = vtkPolyDataCollection::New();
    QStringList surfaceNames;
    for(int j = 0; j < numberOfSurfaces; j ++)
    {
        vtkSmartPointer<vtkFloatArray> b = vtkFloatArray::New();
            b->SetNumberOfValues(mode);
            b->FillComponent(0, 0.0);
            b->SetValue(mode-1, weightValue); //!< View mode
            printInfo("Generating surface for mode " + QString::number(mode-1));

        vtkPolyData *shape = vtkPolyData::New();
            shape->DeepCopy(m_SSM->GetParameterisedShape(b));

        meshes->AddItem(shape);
        surfaceNames.append(QString::number(weightStep));

        weightValue += weightStep;
    }

    printInfo("Mode Surfaces Generated");
    emit collectionAvailable(meshes, surfaceNames);
}

void milxQtShapeModel::generateVectorField()
{
    if(!m_loaded)
        return;
    if(!m_modes)
        generateModes();
    if(!m_modes) //Failed or cancelled, stop
        return;
    if(!m_meaned)
        generateMeanModel();

    printInfo("Generating Vector Field.");
    vtkSmartPointer<vtkPolyData> meanShape = m_meanModel->GetOutput();

    if(!m_modesVectorModel)
        m_modesVectorModel = new milxQtModel; //smart deletion
        m_modesVectorModel->setName("Modes As Vectors");
        m_modesVectorModel->setLargeDataSetMode(true); //May Improve performance for large dataset
        m_modesVectorModel->SetInput(meanShape);
        m_modesVectorModel->generateVectorField(1.0);
        m_modesVectorModel->ImmediateModeRenderingOn(); //May Improve performance for large dataset

    AddActor(m_modesVectorModel->GetActor());
    m_vector = true;
    actionModesAsVectors->setChecked(true);
}

void milxQtShapeModel::generateTensorField()
{
    if(!m_loaded)
        return;
    if(!m_modes)
        generateModes();
    if(!m_modes) //Failed or cancelled, stop
        return;
    if(!m_meaned)
        generateMeanModel();

    printInfo("Generating Tensor Field.");
    if(!m_modesTensorModel)
        m_modesTensorModel = new milxQtModel; //smart deletion
        m_modesTensorModel->setName("Modes As Tensors");
        m_modesTensorModel->setLargeDataSetMode(true); //May Improve performance for large dataset
        m_modesTensorModel->SetInput(m_meanModel->GetOutput());
        m_modesTensorModel->generateTensorField(1.0);
        m_modesTensorModel->ImmediateModeRenderingOn(); //May Improve performance for large dataset

    AddActor(m_modesTensorModel->GetActor());
    m_tensor = true;
    actionModesAsTensors->setChecked(true);
}

void milxQtShapeModel::generateCorrespondences()
{
    if(!m_loaded)
        return;

    const int n = m_SSM->GetNumberOfShapes();
    const int noOfPoints = m_SSM->GetShape(0)->GetNumberOfPoints(); ///\todo Assume all shapes have the smae number of points

    printInfo("Generating Correspondences.");
    working(-1);
    vtkSmartPointer<vtkDoubleArray> maxDeviations = vtkSmartPointer<vtkDoubleArray>::New();
        maxDeviations->SetNumberOfComponents(3);
        maxDeviations->SetNumberOfTuples(noOfPoints);
        maxDeviations->SetName("Deviations");

    ///Initialise deviation array
    coordinate origin(0.0);
    for(int l = 0; l < noOfPoints; l ++) ///For all shapes
        maxDeviations->SetTupleValue(l, origin.data_block());

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        vtkSmartPointer<vtkPolyData> currentShape = m_SSM->GetProcrustesAlignedSurface(j);

        for(int k = 0; k < n; k ++) ///For each shape
        {
            if(j == k)
                continue;

            vtkSmartPointer<vtkPolyData> nextShape = m_SSM->GetProcrustesAlignedSurface(k);

            for(int l = 0; l < noOfPoints; l ++) ///For all shapes
            {
                coordinate currentPoint(currentShape->GetPoint(l));
                coordinate nextPoint(nextShape->GetPoint(l));
                coordinate maxVector(maxDeviations->GetTuple(l));
                coordinate deviation = currentPoint - nextPoint;

                double length = deviation.squared_magnitude();
                double currentMaxLength = maxVector.squared_magnitude();

                if(length > currentMaxLength)
                {
                    //deviation += currentPoint;
                    maxDeviations->SetTupleValue(l, deviation.data_block());
                }
            }

            qApp->processEvents(); ///Keep UI responsive
        }
    }

    printInfo("Showing Hedgehog plot.");
    if(!m_correspondences)
        m_correspondences = new milxQtModel; //smart deletion
        m_correspondences->SetInput(m_SSM->GetMeanShape());
        m_correspondences->SetVectors(maxDeviations);
        m_correspondences->generateHedgehog();

    AddActor(m_correspondences->GetActor());
    m_correspond = true;
    actionCorrespond->setChecked(true);
    done(-1);
}

//Slots
void milxQtShapeModel::mean()
{
    if(!m_loaded && !m_modelled)
        return;

    if(actionMean->isChecked())
    {
        if(!m_meaned)
            generateMeanModel();

        AddActor(m_meanModel->GetActor()); //!< Add data to general display
    }
    else
    {
        if(m_meaned)
            RemoveActor(m_meanModel->GetActor());
    }

    milxQtRenderWindow::Render();
}

void milxQtShapeModel::aligned()
{
    if(!m_loaded && !m_modelled)
        return;

    int n = m_SSM->GetNumberOfShapes();

    if(actionAligned->isChecked())
    {
        if(!m_alignedShapesModelled)
            generateAlignedModels();
        else
            for(int j = 0; j < n; j ++)
                AddActor(m_alignedModels[j]->GetActor()); //!< Add data to general display
    }
    else
    {
        if(m_alignedShapesModelled)
        {
            for(int j = 0; j < n; j ++)
                RemoveActor(m_alignedModels[j]->GetActor());
        }
    }

    milxQtRenderWindow::Render();
}

void milxQtShapeModel::original()
{
    if(!m_loaded && !m_modelled)
        return;

    const int n = m_SSM->GetNumberOfShapes();

    if(actionOriginal->isChecked())
    {
        if(!m_shapesModelled)
            generateModels();
        else
            for(int j = 0; j < n; j ++)
                AddActor(m_models[j]->GetActor()); //!< Add data to general display
    }
    else
    {
        if(m_shapesModelled)
        {
            for(int j = 0; j < n; j ++)
                RemoveActor(m_models[j]->GetActor());
        }
    }

    milxQtRenderWindow::Render();
}

void milxQtShapeModel::modesAsVectors()
{
    if(!m_loaded && !m_modelled)
        return;

    if(actionModesAsVectors->isChecked())
    {
        m_modes = false; //Force re-computation of modes
        generateVectorField();

        if(m_vector)
            AddActor(m_modesVectorModel->GetActor());
    }
    else
    {
        if(m_modes && m_vector)
            RemoveActor(m_modesVectorModel->GetActor());
    }

    milxQtRenderWindow::Render();
}

void milxQtShapeModel::modesAsTensors()
{
    if(!m_loaded && !m_modelled)
        return;

    if(actionModesAsTensors->isChecked())
    {
        m_modes = false; //Force re-computation of modes
        generateTensorField();

        if(m_tensor)
            AddActor(m_modesTensorModel->GetActor());
    }
    else
    {
        if(m_modes && m_tensor)
            RemoveActor(m_modesTensorModel->GetActor());
    }

    milxQtRenderWindow::Render();
}

void milxQtShapeModel::modesAsCollection()
{
    if(!m_loaded && !m_modelled)
        return;

    generateCollectionBasedOnMode();
}

void milxQtShapeModel::correspondences()
{
    if(!m_loaded)
        return;

    if(actionCorrespond->isChecked())
    {
        if(!m_correspond)
            generateCorrespondences();

        AddActor(m_correspondences->GetActor()); //!< Add data to general display
    }
    else
    {
        if(m_correspond)
            RemoveActor(m_correspondences->GetActor());
    }

    milxQtRenderWindow::Render();
}

void milxQtShapeModel::compactness()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> values = m_StandardSSM->GetPCA()->GetEvals();

    ///Init table
    for(int j = 0; j < 3; j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < values->GetNumberOfTuples(); k ++)
            column->InsertNextValue(0.0);
        table->AddColumn(column);
    }
    printDebug("Compactness Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    ///Set column names
    table->GetColumn(0)->SetName("Mode");
    table->GetColumn(1)->SetName("Value");
    table->GetColumn(2)->SetName("Precision");

    ///Compute total of eigenvalues
    double total = 0.0;
    for (int j = 0; j < values->GetNumberOfTuples(); j++)
        total += values->GetValue(j);
    printInfo("Total of eigenvalues is " + QString::number(total));

    ///Compute compactness
    double compactness = 0.0;
    for (int j = 0; j < values->GetNumberOfTuples(); j++)
    {
        compactness += values->GetValue(j)/total;

        table->SetValue(j, 0, j+1);
        table->SetValue(j, 1, values->GetValue(j));
        table->SetValue(j, 2, compactness);

        qApp->processEvents();
    }
    table->Update();
    table->Dump();

    QPointer<milxQtPlot> plot = new milxQtPlot;
        plot->setPlotType2D(0, 2);
        plot->setName("Compactness of the Shape Model");
        plot->setConsole(console);
        plot->legend(false);
        plot->SetSource(table);
        plot->generatePlot(); //emits result, not connected
//        plot->scatterPlot(table, 0, 2);

    emit resultAvailable( qobject_cast<milxQtRenderWindow*>(plot) );
}

void milxQtShapeModel::specificity()
{
    if(!m_modelled)
        generateSSM();

    srand(time(NULL));

    bool ok;
    int n = QInputDialog::getInt(this, tr("Please Provide the number of random shapes to use"),
                                          tr("Number of Random Shapes:"), 30, 2, 1000000, 1, &ok);

    if(!ok)
        return;

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < 3; j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < n; k ++)
            column->InsertNextValue(0.0);
        table->AddColumn(column);
    }
    printDebug("Specificity Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    ///Set column names
    table->GetColumn(0)->SetName("Shape");
    table->GetColumn(1)->SetName("Minimum RMSE");
    table->GetColumn(2)->SetName("Specificity");

    working(-1);
    double specific = 0;
    for(int i = 0; i < n; i++)
    {
      vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
          b->SetNumberOfComponents(1);
          b->SetNumberOfTuples(m_SSM->GetNumberOfModes());
          for(int j = 0; j < m_SSM->GetNumberOfModes(); j++)
              b->SetTuple1(j, rand()%35/10 - 1.75);

      vtkSmartPointer<vtkPolyData> shape = m_SSM->GetParameterisedShape(b);
      double minSumSquares = std::numeric_limits<double>::max();
      for(int k = 0; k < m_SSM->GetNumberOfShapes(); k ++)
      {
          double value = m_SSM->PolyDataDistance(shape, m_SSM->GetShape(k))/m_SSM->GetShape(0)->GetNumberOfPoints();
          if(value < minSumSquares)
              minSumSquares = value;
      }
      specific += minSumSquares;

      table->SetValue(i, 0, i);
      table->SetValue(i, 1, minSumSquares);
      table->SetValue(i, 2, specific);

      qApp->processEvents();
    }
    specific /= n;
    printInfo("Specificity Score: " + QString::number(specific));
    table->Update();
    table->Dump();
    done(-1);

    QPointer<milxQtPlot> plot = new milxQtPlot;
        plot->setPlotType2D(0, 1);
        plot->setName("Specificity of the Shape Model");
        plot->setConsole(console);
        plot->legend(false);
        plot->SetSource(table);
        plot->generatePlot(); //emits result, not connected
//        plot->scatterPlot(table, 0, 1);

    emit resultAvailable( qobject_cast<milxQtRenderWindow*>(plot) );
}

void milxQtShapeModel::generalisability()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < 2; j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < m_SSM->GetNumberOfShapes(); k ++)
            column->InsertNextValue(0.0);
        table->AddColumn(column);
    }
    printDebug("Generalisability Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    ///Set column names
    table->GetColumn(0)->SetName("Shape");
    table->GetColumn(1)->SetName("Reconstruction Error");

    working(-1);
    double generalisability = 0.0;
    for(int i = 0; i < m_SSM->GetNumberOfShapes(); i++)
    {
      /// Leave one out experiment
      vtkPolyData *shape = vtkPolyData::New(); //can't use smart ptrs here
      shape->DeepCopy(m_SSM->GetShape(0));
      m_SSM->RemoveShape(0);
      m_SSM->Update();

      vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
          b->SetNumberOfComponents(1);
          b->SetNumberOfTuples(m_SSM->GetNumberOfModes());
          b->FillComponent(0, 0.0);

      double tx, ty, tz, scale, theta, phi, psi;
      m_SSM->GetSurfaceSimilarityParameters(shape, scale, tx, ty, tz, theta, phi, psi, b);
      vtkSmartPointer<vtkPolyData> shape2 = m_SSM->GetSurface(scale, tx, ty, tz, theta, phi, psi, b);
      double epsiSqr = m_SSM->PolyDataDistance(shape, shape2)/shape->GetNumberOfPoints();

      table->SetValue(i, 0, i);
      table->SetValue(i, 1, epsiSqr);

      generalisability += epsiSqr;
      m_SSM->AddShape(shape);

      qApp->processEvents();
    }
    printInfo("Generalisability Score: " + QString::number(generalisability));
    table->Update();
    table->Dump();
    done(-1);

    QPointer<milxQtPlot> plot = new milxQtPlot;
        plot->setPlotType2D(0, 1);
        plot->setName("Generalisability of the Shape Model");
        plot->setConsole(console);
        plot->legend(false);
        plot->SetSource(table);
        plot->generatePlot(); //emits result, not connected
//        plot->scatterPlot(table, 0, 1);

    emit resultAvailable( qobject_cast<milxQtRenderWindow*>(plot) );
}

void milxQtShapeModel::eigenvalues()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> values = m_StandardSSM->GetPCA()->GetEvals();

    ///Init table
    for(int j = 0; j < 3; j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < values->GetNumberOfTuples(); k ++)
            column->InsertNextValue(0.0);
        table->AddColumn(column);
    }
    printDebug("Compactness Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    ///Set column names
    table->GetColumn(0)->SetName("Mode");
    table->GetColumn(1)->SetName("Value");
    table->GetColumn(2)->SetName("Precision");

    ///Compute total of eigenvalues
    double total = 0.0;
    for (int j = 0; j < values->GetNumberOfTuples(); j++)
        total += values->GetValue(j);
    printInfo("Total of eigenvalues is " + QString::number(total));

    ///Compute compactness
    double compactness = 0.0;
    for (int j = 0; j < values->GetNumberOfTuples(); j++)
    {
        compactness += values->GetValue(j)/total;

        table->SetValue(j, 0, j+1);
        table->SetValue(j, 1, values->GetValue(j));
        table->SetValue(j, 2, compactness);

        qApp->processEvents();
    }
    table->Update();
    table->Dump();

    QPointer<milxQtPlot> plot = new milxQtPlot;
        plot->setPlotType2D(0, 1);
        plot->setName("Eigenvalues of the Shape Model");
        plot->setConsole(console);
        plot->legend(false);
        plot->SetSource(table);
        plot->generatePlot(); //emits result, not connected
//        plot->scatterPlot(table, 0, 2);

    emit resultAvailable( qobject_cast<milxQtRenderWindow*>(plot) );
}

void milxQtShapeModel::eigenmodes()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < m_SSM->GetNumberOfModes(); j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < m_SSM->GetNumberOfShapes(); k ++)
            column->InsertNextValue(0.0);
        const QString nameOfColumn = "Mode " + QString::number(j+1);
        column->SetName(nameOfColumn.toStdString().c_str());
        table->AddColumn(column);
    }
    printDebug("Modes Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    for(int j = 0; j < m_SSM->GetNumberOfShapes(); j ++)
    {
        double scale, tx, ty, tz, theta, phi, psi; //Unused
        vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
            b->SetNumberOfComponents(1);
            b->SetNumberOfTuples(m_SSM->GetNumberOfModes());

        m_SSM->GetSurfaceSimilarityParameters(m_SSM->GetShape(j), scale, tx, ty, tz, theta, phi, psi, b);

        for(int k = 0; k < m_SSM->GetNumberOfModes(); k ++)
            table->SetValue(j, k, b->GetValue(k));

        qApp->processEvents();
    }

    QPointer<milxQtPlot> plot = new milxQtPlot;
        plot->setPlotType3D(0, 1, 2);
        plot->setName("Specificity of the Shape Model");
        plot->setConsole(console);
        plot->legend(false);
        plot->SetSource(table);
        plot->generatePlot(); //emits result, not connected
//        plot->scatterPlot(table, 0, 1, 2);

    emit resultAvailable(plot);
}

void milxQtShapeModel::parameters()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < m_SSM->GetNumberOfModes(); j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < m_SSM->GetNumberOfShapes(); k ++)
            column->InsertNextValue(0.0);
        const QString nameOfColumn = "Mode " + QString::number(j+1);
        column->SetName(nameOfColumn.toStdString().c_str());
        table->AddColumn(column);
    }
    printDebug("Parameters Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    for(int j = 0; j < m_SSM->GetNumberOfShapes(); j ++)
    {
        double scale, tx, ty, tz, theta, phi, psi; //Unused
        vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
            b->SetNumberOfComponents(1);
            b->SetNumberOfTuples(m_SSM->GetNumberOfModes());

        m_SSM->GetSurfaceSimilarityParameters(m_SSM->GetShape(j), scale, tx, ty, tz, theta, phi, psi, b);

        for(int k = 0; k < m_SSM->GetNumberOfModes(); k ++)
            table->SetValue(j, k, b->GetValue(k));

        qApp->processEvents();
    }

    QPointer<milxQtPlot> plot = new milxQtPlot;
        plot->setPlotTypeSurface();
        plot->setName("Training Shape Parameters of the Shape Model");
        plot->setConsole(console);
        plot->legend(false);
        plot->SetSource(table);
        plot->generatePlot(); //emits result, not connected
//        plot->surfacePlot(table);

    emit resultAvailable(plot);
}

void milxQtShapeModel::procrustes()
{
    //Change the Procrustes mode of the SSM based on option set
    if(actionRigid->isChecked())
        m_SSM->SetPoseType(1); //rigid
    else if(actionSimilarity->isChecked())
        m_SSM->SetPoseType(2); //similarity
    else
        m_SSM->SetPoseType(3); //affine

    m_SSM->RemoveProcrustesAlignedPoints();
    //same as m_SSM->SetValid(false); //force recomputation
    m_modelled = false;
    m_meaned = false;
    m_modes = false;
    m_alignedShapesModelled = false;
    generateMeanModel();
    aligned();
    modesAsVectors();
    modesAsTensors();
}

void milxQtShapeModel::alignment()
{
    if(!m_loaded && !m_modelled)
        return;

    working(-1);
    milxQtRenderWindow::renderer->InteractiveOff(); //!<Pause interaction so camera position is fixed

    if(!m_alignedShapesModelled)
    {
        generateAlignedModels(false);

        QMessageBox msgBox;
            msgBox.setText("Models now aligned. Set the window to desired view and try again.");
            msgBox.exec();

        milxQtRenderWindow::renderer->InteractiveOn(); //!<Resume interaction
        done(-1);
        return;
    }

    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString(), filenamePrefix;

    QFileDialog *fileSaver = new QFileDialog(this);

    filenamePrefix = fileSaver->getSaveFileName(this,
                                          tr("Select File Name Prefix to Save"),
                                          path,
                                          tr("All Files (*.*)"));

    if(filenamePrefix.isEmpty())
    {
        done(-1);
        return;
    }

    QFileInfo fi(filenamePrefix);
    QString filename = fi.path()+ "/" + fi.baseName(); //!< Strip the extension

    const int n = m_SSM->GetNumberOfShapes();

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        m_alignedModels[j]->SetOpacity(1.0);
            m_alignedModels[j]->generateSurface(); ///Show as surface

        qApp->processEvents(); ///Keep UI responsive
    }

    outputSnapshots(filename);

    milxQtRenderWindow::renderer->InteractiveOn(); //!<Resume interaction
    done(-1);
}

void milxQtShapeModel::alignedMeshes()
{
    if(!m_loaded && !m_modelled)
        return;

    working(-1);

    if(!m_alignedShapesModelled)
        generateAlignedModels(false);

    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString(), filenamePrefix;

    QFileDialog *fileSaver = new QFileDialog(this);

    filenamePrefix = fileSaver->getSaveFileName(this,
                                          tr("Select File Name Prefix to Save"),
                                          path,
                                          tr("All Files (*.*)"));

    if(filenamePrefix.isEmpty())
    {
        done(-1);
        return;
    }

    QFileInfo fi(filenamePrefix);
    QString filename = fi.path()+ "/" + fi.baseName(); //!< Strip the extension
    QString ext = fi.suffix();
    if(ext.isEmpty())
        ext = "vtk";

    const int n = m_SSM->GetNumberOfShapes();

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        m_alignedModels[j]->SetOpacity(1.0);
            m_alignedModels[j]->generateSurface(); ///Show as surface

        QPointer<milxQtFile> writer = new milxQtFile;
            QString paddedNumber = QString("%1").arg(j, 3, 10, QChar('0')); //Zero pad the number
            QString newFilename = filename + paddedNumber + "." + ext;
            writer->saveModel(newFilename, m_alignedModels[j]);

        qApp->processEvents(); ///Keep UI responsive
    }

    done(-1);
}

void milxQtShapeModel::originalMeshes()
{
    if(!m_loaded && !m_modelled)
        return;

    working(-1);

    if(!m_shapesModelled)
        generateModels(false);

    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString(), filenamePrefix;

    QFileDialog *fileSaver = new QFileDialog(this);

    filenamePrefix = fileSaver->getSaveFileName(this,
                                          tr("Select File Name Prefix to Save"),
                                          path,
                                          tr("All Files (*.*)"));

    if(filenamePrefix.isEmpty())
    {
        done(-1);
        return;
    }

    QFileInfo fi(filenamePrefix);
    QString filename = fi.path()+ "/" + fi.baseName(); //!< Strip the extension
    QString ext = fi.suffix();
    if(ext.isEmpty())
        ext = "vtk";

    const int n = m_SSM->GetNumberOfShapes();

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        m_models[j]->SetOpacity(1.0);
            m_models[j]->generateSurface(); ///Show as surface

        QPointer<milxQtFile> writer = new milxQtFile;
            QString paddedNumber = QString("%1").arg(j, 3, 10, QChar('0')); //Zero pad the number
            QString newFilename = filename + paddedNumber + "." + ext;
            writer->saveModel(newFilename, m_models[j]);

        qApp->processEvents(); ///Keep UI responsive
    }

    done(-1);
}

void milxQtShapeModel::pointIds()
{
    if(!m_loaded && !m_modelled)
        return;

    working(-1);
    milxQtRenderWindow::renderer->InteractiveOff(); //!<Pause interaction so camera position is fixed

    if(!m_alignedShapesModelled)
    {
        generateAlignedModels(false);

        QMessageBox msgBox;
            msgBox.setText("Models now aligned. Set the window to desired view and try again.");
            msgBox.exec();

        milxQtRenderWindow::renderer->InteractiveOn(); //!<Resume interaction
        done(-1);
        return;
    }

    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString(), filenamePrefix;

    QFileDialog *fileSaver = new QFileDialog(this);

    filenamePrefix = fileSaver->getSaveFileName(this,
                                          tr("Select File Name Prefix to Save"),
                                          path,
                                          tr("All Files (*.*)"));

    if(filenamePrefix.isEmpty())
    {
        done(-1);
        return;
    }

    QFileInfo fi(filenamePrefix);
    QString filename = fi.path()+ "/" + fi.baseName(); //!< Strip the extension

    const int noOfPoints = m_SSM->GetShape(0)->GetNumberOfPoints();
    const int n = m_SSM->GetNumberOfShapes();

    ///Create scalars based on the point ids
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
        colorLookupTable->SetTableRange(0, noOfPoints);
        colorLookupTable->Build();

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        colors->SetNumberOfComponents(3);
        colors->SetName("Colors");

    for(int j = 0; j < noOfPoints; j ++)
    {
        double dcolor[3];

        colorLookupTable->GetColor(j, dcolor);

        unsigned char color[3];

        color[0] = static_cast<unsigned char>(255.0 * dcolor[0]);
        color[1] = static_cast<unsigned char>(255.0 * dcolor[1]);
        color[2] = static_cast<unsigned char>(255.0 * dcolor[2]);

        colors->InsertNextTupleValue(color);
    }

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        m_alignedModels[j]->SetScalars(colors); ///Set the point ID colours
            m_alignedModels[j]->SetOpacity(1.0);
            m_alignedModels[j]->generateSurface(); ///Show as surface

        qApp->processEvents(); ///Keep UI responsive
    }

    outputSnapshots(filename);

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        m_alignedModels[j]->SetOpacity(0.1);
            m_alignedModels[j]->generatePoints(); ///Revert to points

        qApp->processEvents(); ///Keep UI responsive
    }

    milxQtRenderWindow::renderer->InteractiveOn(); //!<Resume interaction
    done(-1);
}

void milxQtShapeModel::coordinates()
{
    if(!m_loaded && !m_modelled)
        return;

    working(-1);
    milxQtRenderWindow::renderer->InteractiveOff(); //!<Pause interaction so camera position is fixed

    if(!m_alignedShapesModelled)
    {
        generateAlignedModels(false);

        QMessageBox msgBox;
            msgBox.setText("Models now aligned. Set the window to desired view and try again.");
            msgBox.exec();

        milxQtRenderWindow::renderer->InteractiveOn(); //!<Resume interaction
        done(-1);
        return;
    }

    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString(), filenamePrefix;

    QFileDialog *fileSaver = new QFileDialog(this);

    filenamePrefix = fileSaver->getSaveFileName(this,
                                          tr("Select File Name Prefix to Save"),
                                          path,
                                          tr("All Files (*.*)"));

    if(filenamePrefix.isEmpty())
    {
        done(-1);
        return;
    }

    QFileInfo fi(filenamePrefix);
    QString filename = fi.path()+ "/" + fi.baseName(); //!< Strip the extension

    const int noOfPoints = m_SSM->GetShape(0)->GetNumberOfPoints();
    const int n = m_SSM->GetNumberOfShapes();

    ///Create scalars based on the point ids
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
        colorLookupTable->SetTableRange(0, 2.0*vtkMath::Pi() * 2.0*vtkMath::Pi()); //2 pi by 2 pi
        colorLookupTable->Build();

    for(int j = 0; j < n; j ++) ///For all shapes
    {
        vtkSmartPointer<vtkPolyData> shape = m_SSM->GetProcrustesAlignedSurface(j);

        vtkSmartPointer<vtkSphericalTransform> transform = vtkSmartPointer<vtkSphericalTransform>::New();

        vtkSmartPointer<vtkTransformPolyDataFilter> transformShape = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
            transformShape->SetInput(shape);
            transformShape->SetTransform(transform->GetInverse());
            transformShape->Update();
            vtkSmartPointer<vtkPolyData> transformedShape = transformShape->GetOutput();

        vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
            colors->SetNumberOfComponents(3);
            colors->SetNumberOfTuples(noOfPoints);
            colors->SetName("Colors");

        for(int k = 0; k < noOfPoints; k ++)
        {
            double dcolor[3], point[3];
            transformedShape->GetPoint(k, point);
            //cerr << "Point: (" << point[0] << ", " << point[1] << ", " << point[2] << ")" << endl;

            colorLookupTable->GetColor(point[1]*point[2], dcolor);

            unsigned char color[3];

            color[0] = static_cast<unsigned char>(255.0 * dcolor[0]);
            color[1] = static_cast<unsigned char>(255.0 * dcolor[1]);
            color[2] = static_cast<unsigned char>(255.0 * dcolor[2]);

            colors->SetTupleValue(k, color);
        }

        m_alignedModels[j]->SetScalars(colors); ///Set the point ID colours
            m_alignedModels[j]->SetOpacity(1.0);
            m_alignedModels[j]->generateSurface(); ///Show as surface

        qApp->processEvents(); ///Keep UI responsive
    }

    outputSnapshots(filename);

    milxQtRenderWindow::renderer->InteractiveOn(); //!<Resume interaction
    done(-1);
}

void milxQtShapeModel::replaceOriginal()
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateAlignedModels(false);
    if(!m_shapesModelled)
        generateModels(false);

    emit working(-1);
    const int n = m_SSM->GetNumberOfShapes();

    printInfo("Replacing Original Models.");
    for(int j = 0; j < n; j ++)
    {
//        m_SSM->SetShape(m_SSM->GetProcrustesAlignedSurface(j), j); //!< Set model to shape filter
        m_SSM->GetShape(j)->SetPoints(m_SSM->GetAlignedPoints(j)); //!< Set model to shape filter

//        coordinate otherCentroid = m_models[j]->centroid();
        double otherCentroidSize = m_models[j]->centroidSize();

        m_models[j]->SetPoints(m_SSM->GetAlignedPoints(j));

        vtkSmartPointer<vtkTransform> transformer = vtkSmartPointer<vtkTransform>::New();
//            transformer->Translate(otherCentroid[0], otherCentroid[1], otherCentroid[2]);
            transformer->Scale(otherCentroidSize, otherCentroidSize, otherCentroidSize);

        m_models[j]->SetTransform(transformer);
        m_models[j]->refresh();

        qApp->processEvents(); ///Keep UI responsive
    }

    done(-1);
}

void milxQtShapeModel::updateLookupTable()
{
    if(m_modesVectorModel)
    {
        m_modesVectorModel->SetLookupTable(lookupTable);
        m_modesVectorModel->updateLookupTable();
    }
    else if(m_modesTensorModel)
    {
        m_modesTensorModel->SetLookupTable(lookupTable);
        m_modesTensorModel->updateLookupTable();
    }
    else if(m_correspondences)
    {
        m_correspondences->SetLookupTable(lookupTable);
        m_correspondences->updateLookupTable();
    }
    else
    {
        m_meanModel->SetLookupTable(lookupTable);
        m_meanModel->updateLookupTable();
    }
}

void milxQtShapeModel::reset()
{
    m_loaded = false;
    m_modelled = false;
    m_meaned = false;
    m_shapesModelled = false;
    m_alignedShapesModelled = false;
    m_modes = false;
    m_vector = false;
    m_tensor = false;
    m_correspond = false;
    m_mode = 1;
}

void milxQtShapeModel::outputSnapshots(const QString filename)
{
    const int n = m_SSM->GetNumberOfShapes();

    QPointer<milxQtRenderWindow> offScreenWin = new milxQtRenderWindow;
        offScreenWin->generateRender(); //Call this, otherwise blank output always
        offScreenWin->SetBackground(1, 1, 1);
        offScreenWin->OffScreenRenderingOn();
        vtkRenderer *ren = offScreenWin->GetRenderer();
        ren->SetActiveCamera(milxQtRenderWindow::renderer->GetActiveCamera()); //!<transfer fixed camera setup

    printInfo("Found " + QString::number(n) + " shapes to check.");
    for(int j = 0; j < n; j ++) ///For all shapes
    {
        QString tmpName = filename + QString::number(j) + ".png";
        printInfo("Saving " + tmpName);

        ren->AddActor(m_alignedModels[j]->GetActor()); //!< Add data to general display
        offScreenWin->Render();

        qApp->processEvents(); ///Keep UI responsive

        vtkSmartPointer<vtkWindowToImageFilter> windowToImage = vtkSmartPointer<vtkWindowToImageFilter>::New();
            windowToImage->SetInput(offScreenWin->GetRenderWindow());
            windowToImage->Update();

        vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New(); //!< Write PNG image
            writer->SetFileName(tmpName.toStdString().c_str());
            writer->SetInput(windowToImage->GetOutput());
            writer->Write();

        ren->RemoveActor(m_alignedModels[j]->GetActor()); //!< Remove data from general display

        qApp->processEvents(); ///Keep UI responsive
    }
}

void milxQtShapeModel::createActions()
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection
    contextMenu->setTitle(QApplication::translate("MainWindow", "Shape Modelling", 0, QApplication::UnicodeUTF8));

    actionMean = new QAction(this);
        actionMean->setText(QApplication::translate("SSM", "&Mean Model", 0, QApplication::UnicodeUTF8));
        actionMean->setShortcut(tr("Alt+m"));
        actionMean->setCheckable(true);
        actionMean->setChecked(false);
    actionAligned = new QAction(this);
        actionAligned->setText(QApplication::translate("SSM", "&Aligned Points", 0, QApplication::UnicodeUTF8));
        actionAligned->setShortcut(tr("Alt+a"));
        actionAligned->setCheckable(true);
        actionAligned->setChecked(false);
    actionOriginal = new QAction(this);
        actionOriginal->setText(QApplication::translate("SSM", "&Original Points", 0, QApplication::UnicodeUTF8));
        actionOriginal->setShortcut(tr("Alt+o"));
        actionOriginal->setCheckable(true);
        actionOriginal->setChecked(false);
    actionModesAsVectors = new QAction(this);
        actionModesAsVectors->setText(QApplication::translate("SSM", "Modes as &Vectors", 0, QApplication::UnicodeUTF8));
        actionModesAsVectors->setShortcut(tr("Alt+v"));
        actionModesAsVectors->setCheckable(true);
        actionModesAsVectors->setChecked(false);
    actionModesAsTensors = new QAction(this);
        actionModesAsTensors->setText(QApplication::translate("SSM", "Modes as &Tensors", 0, QApplication::UnicodeUTF8));
        actionModesAsTensors->setShortcut(tr("Alt+t"));
        actionModesAsTensors->setCheckable(true);
        actionModesAsTensors->setChecked(false);
    actionModesAsCollection = new QAction(this);
        actionModesAsCollection->setText(QApplication::translate("SSM", "Modes as Surface &Collection", 0, QApplication::UnicodeUTF8));
        actionModesAsCollection->setShortcut(tr("Alt+c"));
    actionCorrespond = new QAction(this);
        actionCorrespond->setText(QApplication::translate("SSM", "Correspondences as a &Hedgehog", 0, QApplication::UnicodeUTF8));
        actionCorrespond->setShortcut(tr("Alt+h"));
        actionCorrespond->setCheckable(true);
        actionCorrespond->setChecked(false);

    //plots menu
    actionCompact = new QAction(this);
        actionCompact->setText(QApplication::translate("SSM", "Compactness", 0, QApplication::UnicodeUTF8));
        actionCompact->setShortcut(tr("Shift+Alt+c"));
    actionSpecificity = new QAction(this);
        actionSpecificity->setText(QApplication::translate("SSM", "Specificity", 0, QApplication::UnicodeUTF8));
        actionSpecificity->setShortcut(tr("Shift+Alt+s"));
    actionGeneralise = new QAction(this);
        actionGeneralise->setText(QApplication::translate("SSM", "Generalisability", 0, QApplication::UnicodeUTF8));
        actionGeneralise->setShortcut(tr("Shift+Alt+g"));
    actionValues = new QAction(this);
        actionValues->setText(QApplication::translate("SSM", "Eigenvalues", 0, QApplication::UnicodeUTF8));
        actionValues->setShortcut(tr("Shift+Alt+v"));
    actionModes = new QAction(this);
        actionModes->setText(QApplication::translate("SSM", "Primary Eigenmodes", 0, QApplication::UnicodeUTF8));
        actionModes->setShortcut(tr("Shift+Alt+r"));
    actionParameters = new QAction(this);
        actionParameters->setText(QApplication::translate("SSM", "Training Shape Parameters", 0, QApplication::UnicodeUTF8));
        actionParameters->setShortcut(tr("Shift+Alt+p"));

    //Procrustes menu
    actionRigid = new QAction(this);
        actionRigid->setText(QApplication::translate("SSM", "Rigid", 0, QApplication::UnicodeUTF8));
        actionRigid->setShortcut(tr("Ctrl+Alt+r"));
        actionRigid->setCheckable(true);
    actionSimilarity = new QAction(this);
        actionSimilarity->setText(QApplication::translate("SSM", "Similarity", 0, QApplication::UnicodeUTF8));
        actionSimilarity->setShortcut(tr("Ctrl+Alt+s"));
        actionSimilarity->setCheckable(true);
        actionSimilarity->setChecked(true);
    actionAffine = new QAction(this);
        actionAffine->setText(QApplication::translate("SSM", "Affine", 0, QApplication::UnicodeUTF8));
        actionAffine->setShortcut(tr("Ctrl+Alt+a"));
        actionAffine->setCheckable(true);
    alignGroup = new QActionGroup(this);
        alignGroup->addAction(actionRigid);
        alignGroup->addAction(actionSimilarity);
        alignGroup->addAction(actionAffine);

    actionAlignment = new QAction(this);
        actionAlignment->setText(QApplication::translate("SSM", "Output &Alignment as Images", 0, QApplication::UnicodeUTF8));
        actionAlignment->setShortcut(tr("Alt+r"));
    actionOriginalMeshes = new QAction(this);
        actionOriginalMeshes->setText(QApplication::translate("SSM", "Output &Original Meshes", 0, QApplication::UnicodeUTF8));
        actionOriginalMeshes->setShortcut(tr("Shift+Alt+o"));
    actionAlignedMeshes = new QAction(this);
        actionAlignedMeshes->setText(QApplication::translate("SSM", "Output &Aligned Meshes", 0, QApplication::UnicodeUTF8));
        actionAlignedMeshes->setShortcut(tr("Shift+Alt+a"));
    actionPointIDs = new QAction(this);
        actionPointIDs->setText(QApplication::translate("SSM", "Output &Point IDs as Images", 0, QApplication::UnicodeUTF8));
        actionPointIDs->setShortcut(tr("Alt+p"));
    actionCoordinates = new QAction(this);
        actionCoordinates->setText(QApplication::translate("SSM", "Output &Coordinates as Images", 0, QApplication::UnicodeUTF8));
        actionCoordinates->setShortcut(tr("Alt+c"));

    actionReplaceOriginal = new QAction(this);
        actionReplaceOriginal->setText(QApplication::translate("SSM", "Replace Originals with Aligned", 0, QApplication::UnicodeUTF8));
        actionReplaceOriginal->setShortcut(tr("Shift+Alt+r"));
}

void milxQtShapeModel::createConnections()
{
    //Operations
    connect(actionMean, SIGNAL(triggered()), this, SLOT(mean()));
    connect(actionAligned, SIGNAL(triggered()), this, SLOT(aligned()));
    connect(actionOriginal, SIGNAL(triggered()), this, SLOT(original()));
    connect(actionModesAsVectors, SIGNAL(triggered()), this, SLOT(modesAsVectors()));
    connect(actionModesAsTensors, SIGNAL(triggered()), this, SLOT(modesAsTensors()));
    connect(actionModesAsCollection, SIGNAL(triggered()), this, SLOT(modesAsCollection()));
    connect(actionCorrespond, SIGNAL(triggered()), this, SLOT(correspondences()));

    //Plot
    connect(actionCompact, SIGNAL(triggered()), this, SLOT(compactness()));
    connect(actionSpecificity, SIGNAL(triggered()), this, SLOT(specificity()));
    connect(actionGeneralise, SIGNAL(triggered()), this, SLOT(generalisability()));
    connect(actionValues, SIGNAL(triggered()), this, SLOT(eigenvalues()));
    connect(actionModes, SIGNAL(triggered()), this, SLOT(eigenmodes()));
    connect(actionParameters, SIGNAL(triggered()), this, SLOT(parameters()));

    //Procrustes
    connect(actionRigid, SIGNAL(triggered()), this, SLOT(procrustes()));
    connect(actionSimilarity, SIGNAL(triggered()), this, SLOT(procrustes()));
    connect(actionAffine, SIGNAL(triggered()), this, SLOT(procrustes()));

    connect(actionAlignment, SIGNAL(triggered()), this, SLOT(alignment()));
    connect(actionOriginalMeshes, SIGNAL(triggered()), this, SLOT(originalMeshes()));
    connect(actionAlignedMeshes, SIGNAL(triggered()), this, SLOT(alignedMeshes()));
    connect(actionPointIDs, SIGNAL(triggered()), this, SLOT(pointIds()));
    connect(actionCoordinates, SIGNAL(triggered()), this, SLOT(coordinates()));
    connect(actionReplaceOriginal, SIGNAL(triggered()), this, SLOT(replaceOriginal()));
}

void milxQtShapeModel::setupTooltips()
{
    actionMean->setToolTip("Show the mean of the Statistical Shape Model");
    actionMean->setStatusTip("Show the mean of the Statistical Shape Model");
    actionAligned->setToolTip("Show the aligned shapes of the Statistical Shape Model");
    actionAligned->setStatusTip("Show the aligned shapes of the Statistical Shape Model");
    actionOriginal->setToolTip("Show the original (unmodified) shapes of the Statistical Shape Model");
    actionOriginal->setStatusTip("Show the original (unmodified) shapes of the Statistical Shape Model");
    actionModesAsVectors->setToolTip("Show the first mode of the Statistical Shape Model as a vector field");
    actionModesAsVectors->setStatusTip("Show the first mode of the Statistical Shape Model as a vector field");
    actionModesAsTensors->setToolTip("Show the first mode of the Statistical Shape Model as an ellipsoidal tensor field");
    actionModesAsTensors->setStatusTip("Show the first mode of the Statistical Shape Model as an ellipsoidal tensor field");
    actionModesAsCollection->setToolTip("Show the modes of the Statistical Shape Model as a collection of surfaces");
    actionModesAsCollection->setStatusTip("Show the modes of the Statistical Shape Model as a collection of surfaces");
    actionCorrespond->setToolTip("Show a measure of the maximum deviation of the correspondences between shapes as a hedgehog plot.");
    actionCorrespond->setStatusTip("Show a measure of the maximum deviation of the correspondences between shapes as a hedgehog plot.");

    //Plot
    actionCompact->setToolTip("Plot the compactness of the shape model.");
    actionCompact->setStatusTip("Plot the compactness of the shape model.");
    actionSpecificity->setToolTip("Plot the specificity of the shape model.");
    actionSpecificity->setStatusTip("Plot the specificity of the shape model.");
    actionGeneralise->setToolTip("Plot the generalisability of the shape model.");
    actionGeneralise->setStatusTip("Plot the generalisability of the shape model.");
    actionValues->setToolTip("Plot the eigenvalues within the shape model.");
    actionValues->setStatusTip("Plot the eigenvalues within the shape model.");
    actionModes->setToolTip("Plot the three primary modes of the training shapes within the shape model.");
    actionModes->setStatusTip("Plot the three primary modes of the training shapes within the shape model.");
    actionParameters->setToolTip("Plot the shape parameters of the training shapes within the shape model.");
    actionParameters->setStatusTip("Plot the shape parameters of the training shapes within the shape model.");

    //Procrustes
    actionRigid->setToolTip("Procrustes alignment of surfaces to rigid landmark transform.");
    actionRigid->setStatusTip("Procrustes alignment of surfaces to rigid landmark transform.");
    actionSimilarity->setToolTip("Procrustes alignment of surfaces to rigid+scale landmark transform.");
    actionSimilarity->setStatusTip("Procrustes alignment of surfaces to rigid+scale landmark transform.");
    actionAffine->setToolTip("Procrustes alignment of surfaces to affine landmark transform.");
    actionAffine->setStatusTip("Procrustes alignment of surfaces to affine landmark transform.");

    actionAlignment->setToolTip("Output each aligned shape as a PNG image via Off-screen rendering.");
    actionAlignment->setStatusTip("Output each aligned shape as a PNG image via Off-screen rendering.");
    actionOriginalMeshes->setToolTip("Output each original model as a poly data file.");
    actionOriginalMeshes->setStatusTip("Output each original model as a poly data file.");
    actionAlignedMeshes->setToolTip("Output each aligned shape as a poly data file.");
    actionAlignedMeshes->setStatusTip("Output each aligned shape as a poly data file.");
    actionPointIDs->setToolTip("Output each aligned shape with point IDs coloured as a PNG image via Off-screen rendering.");
    actionPointIDs->setStatusTip("Output each aligned shape with point IDs coloured as a PNG image via Off-screen rendering.");
    actionCoordinates->setToolTip("Output the coordinates (as scalars) of each shape as a PNG image via Off-screen rendering.");
    actionCoordinates->setStatusTip("Output the coordinates (as scalars) of each shape as a PNG image via Off-screen rendering.");
    actionReplaceOriginal->setToolTip("Replace the original meshes with the aligned ones within the SSM. Aligned meshes are rescaled to original mesh size.");
    actionReplaceOriginal->setStatusTip("Replace the original meshes with the aligned ones within the SSM. Aligned meshes are rescaled to original mesh size.");
}

void milxQtShapeModel::contextMenuEvent(QContextMenuEvent *contextEvent)
{
    createMenu(contextMenu);

    contextMenu->exec(contextEvent->globalPos());
}
