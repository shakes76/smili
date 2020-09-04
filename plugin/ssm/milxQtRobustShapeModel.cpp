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

#include "milxQtRobustShapeModel.h"

//Qt
#include <QMenu>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <QFormLayout>
#include <QComboBox>
#include <QPushButton>

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

milxQtRobustShapeModel::milxQtRobustShapeModel(QWidget *theParent, bool contextSystem) : milxQtShapeModel(theParent, contextSystem)
{
    reset();

    m_RobustSSM = RobustShapeModelType::New();
        m_SSM = m_RobustSSM;
        m_RobustSSM->DebugOn();

    m_meanModel = NULL;
    m_modesVectorModel = NULL;
    m_modesTensorModel = NULL;
    m_correspondences = NULL;

    ///Set strings
    milxQtWindow::prefix = "RSSM: ";

    createActions();
    setupTooltips();
    createConnections();
}

milxQtRobustShapeModel::~milxQtRobustShapeModel()
{

}

void milxQtRobustShapeModel::LoadModel(const QString filename)
{
    reset();

    if(m_RobustSSM->LoadCompactModel(filename.toStdString().c_str()))
        m_loaded = true;
    else
    {
        printError("Could not open Robust SSM file.");
        return;
    }

    generateSSM();
}

void milxQtRobustShapeModel::SetInputCollection(vtkPolyDataCollection* meshes)
{
    size_t n = meshes->GetNumberOfItems();

    working(-1);
    reset();

    meshes->InitTraversal();
    for(size_t j = 0; j < n; j ++)
    {
	vtkPolyData * shape = meshes->GetNextItem();

	vtkFloatArray *weights = vtkFloatArray::New(); //owned by RSSM class
	    weights->SetNumberOfTuples(shape->GetNumberOfPoints());
	    weights->SetNumberOfComponents(1);
	    weights->FillComponent(0, 1.0);

	m_RobustSSM->AddShape(shape, weights); //!< Set model to shape filter

        qApp->processEvents(); ///Keep UI responsive
    }
    done(-1);

    m_loaded = true;
}

void milxQtRobustShapeModel::SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames)
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

	vtkPolyData * shape = meshes->GetNextItem();

	vtkFloatArray *weights = vtkFloatArray::New(); //owned by RSSM class
	    weights->SetNumberOfTuples(shape->GetNumberOfPoints());
	    weights->SetNumberOfComponents(1);
	    weights->FillComponent(0, 1.0);

	m_RobustSSM->AddShape(shape, weights); //!< Set model to shape filter

        qApp->processEvents(); ///Keep UI responsive
    }
    done(-1);

    qDebug() << "Case IDs: " << m_caseIDs << endl;

    m_loaded = true;
}

void milxQtRobustShapeModel::SetInputCollection(vtkPolyDataCollection* meshes, vtkPolyData *atlasSurface, QStringList &filenames)
{
  const int n = meshes->GetNumberOfItems();

  working(-1);
  reset();

  if(!atlasSurface->GetPointData()->GetScalars())
  {
      printError("No Scalars present on atlas surface. Aborting.");
      return;
  }

  vtkFloatArray *atlasWeights = vtkFloatArray::SafeDownCast(atlasSurface->GetPointData()->GetScalars());
  if(!atlasWeights)
      printWarning("Array not floating point type?");

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

  //Set Focused mode
  printInfo("Using Focused Model Mode");
  m_RobustSSM->SetTolerance(0);
  m_RobustSSM->RandomInitialisationOff();
  m_RobustSSM->PreAlignedOff();
  m_RobustSSM->MissingModeOff();
  m_RobustSSM->SetMaxIterations(3);

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

      vtkPolyData * shape = meshes->GetNextItem();

      vtkFloatArray *weights = vtkFloatArray::New(); //owned by RSSM class
      weights->DeepCopy(atlasWeights);

      m_RobustSSM->AddShape(shape, weights); //!< Set model to shape filter

      qApp->processEvents(); ///Keep UI responsive
    }
  done(-1);

  qDebug() << "Case IDs: " << m_caseIDs << endl;

  m_loaded = true;
}

void milxQtRobustShapeModel::createMenu(QMenu *menu)
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

void milxQtRobustShapeModel::generateSSM()
{
    if(!m_loaded)
        return;

    bool ok;
    float precision = QInputDialog::getDouble(this, tr("Please Provide the precision of the model"),
                              tr("Precision:"), 0.9, 0.0, 1.0, 5, &ok);

    if(!ok) //cancelled
        return;

    ///Compute SSM
    try
    {
        m_RobustSSM->SetPrecision(precision);
        m_RobustSSM->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        printError("Exception caught while generating an SSM!");
        cout << err << endl;
        return;
    }

    const size_t n = m_RobustSSM->GetNumberOfShapes();

    ///Compute the sum of eigenvalues, which is a measure of the variation of the modelling
    double eigenSum = 0.0;
    vtkSmartPointer<vtkFloatArray> eigenVals = m_RobustSSM->GetPCA()->GetEigenValues();
        printInfo("Eigenvalues Info - #tuples: " + QString::number(eigenVals->GetNumberOfTuples()) + ", #components: " + QString::number(eigenVals->GetNumberOfComponents()));
        for(vtkIdType j = 0; j < eigenVals->GetNumberOfTuples(); j ++)
            eigenSum += eigenVals->GetValue(j);
        printInfo("Sum of Eigenvalues is " + QString::number(eigenSum));
        printInfo("Mean Scale: " + QString::number(m_RobustSSM->GetMeanShapeSize()));
        printInfo("Number of Modes: " + QString::number(m_RobustSSM->GetNumberOfModes()));

    vtkSmartPointer<vtkLookupTable> tmpLookupTable = vtkLookupTable::SafeDownCast(milxQtRenderWindow::lookupTable);
        tmpLookupTable->SetTableRange(0.0, n+1); ///Build colour lookup table
        tmpLookupTable->Build();
    milxQtRenderWindow::lookupTable = tmpLookupTable;

    milxQtRenderWindow::setView(milxQtRenderWindow::defaultView); //Default view

    m_modelled = true;
}

void milxQtRobustShapeModel::generateMeanModel(vtkSmartPointer<vtkPolyData> shape)
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();

    printInfo("Generating Mean Model.");
    if(!m_meanModel)
    {
        m_meanModel = new milxQtModel;
        m_meanModel->SetInput(m_RobustSSM->GetMeanShape());
        m_meanModel->SetScalars(m_RobustSSM->GetWeights(0));
    }
    if(!shape)
    {
        m_meanModel->setName("Mean");
        m_meanModel->SetPoints(m_RobustSSM->GetMeanShape()->GetPoints());
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

void milxQtRobustShapeModel::generateModes()
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
                              tr("Mode:"), 1, 1, m_RobustSSM->GetNumberOfModes(), 1, &ok1);
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
    meanShape->SetPoints(m_RobustSSM->GetMeanShape()->GetPoints());

    vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
#if VTK_MAJOR_VERSION <= 5
        normals->SetInput(meanShape);
#else
        normals->SetInputData(meanShape);
#endif
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
        varShape->DeepCopy(m_RobustSSM->GetParameterisedShape(b));

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
        double* meanPoint = meanShape->GetPoint(i);
        double* varPoint  = varShape->GetPoint(i);

        double xVal = varPoint[0] - meanPoint[0];
        double yVal = varPoint[1] - meanPoint[1];
        double zVal = varPoint[2] - meanPoint[2];

        float normalPoint[3];
        normArray->GetTupleValue(i, normalPoint);
        double var;
        if(flgNormal)
        {
            var = xVal*normalPoint[0] + yVal*normalPoint[1] + zVal*normalPoint[2];
            var = var*var;
        }
        else
        {
            var = xVal*xVal + yVal*yVal + zVal*zVal;
        }
//        std::cout << var << std::endl;
        varVectors->InsertNextTuple3(xVal, yVal, zVal);
        varShapeVectors->InsertNextTuple3(-xVal, -yVal, -zVal);
        varScalars->InsertNextValue(var);

        vnl_scatter_3x3<float> C;
        for(int j = 0; j < m_RobustSSM->GetNumberOfShapes(); j++)
        {
            double point[3];
            m_RobustSSM->GetShape(j)->GetPoint(i, point);
            // std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
            vnl_vector_fixed<float, 3> point_vector;
            point_vector[0] = point[0] - meanPoint[0];
            point_vector[1] = point[1] - meanPoint[1];
            point_vector[2] = point[2] - meanPoint[2];
            // std::cout << point_vector << std::endl;
            C.add_outer_product(point_vector);
        }
        // C /= m_RobustSSM->GetNumberOfShapes();
        // std::cout << meanPoint[0] << " " << meanPoint[1] << " " << meanPoint[2] << std::endl << std::endl;
//        std::cout << C << std::endl;
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
//        std::cout << tens->T[0] << " " << tens->T[1] << " " << tens->T[2] << std::endl;

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

void milxQtRobustShapeModel::generateCollectionBasedOnMode()
{
    if(!m_loaded)
        return;
    if(!m_modelled)
        generateSSM();
    if(!m_meaned)
        generateMeanModel();

    const int n = m_RobustSSM->GetNumberOfShapes();
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
            shape->DeepCopy(m_RobustSSM->GetParameterisedShape(b));

        meshes->AddItem(shape);
        surfaceNames.append(QString::number(weightStep));

        weightValue += weightStep;
    }

    printInfo("Mode Surfaces Generated");
    emit collectionAvailable(meshes, surfaceNames);
}

void milxQtRobustShapeModel::generateCorrespondences()
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
void milxQtRobustShapeModel::compactness()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> values = m_RobustSSM->GetPCA()->GetEvals();

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
    //table->Update();
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

void milxQtRobustShapeModel::specificity()
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
          b->SetNumberOfTuples(m_RobustSSM->GetNumberOfModes());
          for(int j = 0; j < m_RobustSSM->GetNumberOfModes(); j++)
              b->SetTuple1(j, rand()%35/10 - 1.75);

      vtkSmartPointer<vtkPolyData> shape = m_RobustSSM->GetParameterisedShape(b);
      double minSumSquares = std::numeric_limits<double>::max();
      for(int k = 0; k < m_RobustSSM->GetNumberOfShapes(); k ++)
      {
          double value = m_RobustSSM->PolyDataDistance(shape, m_RobustSSM->GetShape(k))/m_RobustSSM->GetShape(0)->GetNumberOfPoints();
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
    //table->Update();
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

void milxQtRobustShapeModel::generalisability()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < 2; j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < m_RobustSSM->GetNumberOfShapes(); k ++)
            column->InsertNextValue(0.0);
        table->AddColumn(column);
    }
    printDebug("Generalisability Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    ///Set column names
    table->GetColumn(0)->SetName("Shape");
    table->GetColumn(1)->SetName("Reconstruction Error");

    working(-1);
    double generalisability = 0.0;
    for(int i = 0; i < m_RobustSSM->GetNumberOfShapes(); i++)
    {
      /// Leave one out experiment
      vtkPolyData *shape = vtkPolyData::New(); //can't use smart ptrs here
      shape->DeepCopy(m_RobustSSM->GetShape(0));
      vtkFloatArray *weights = vtkFloatArray::New(); //can't use smart ptrs here
      weights->DeepCopy(m_RobustSSM->GetWeights(0));
      m_RobustSSM->RemoveShape(0); //removes weights as well
      m_RobustSSM->Update();

      vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
          b->SetNumberOfComponents(1);
          b->SetNumberOfTuples(m_RobustSSM->GetNumberOfModes());
          b->FillComponent(0, 0.0);

      vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
      shape->GetPointData()->SetScalars(weights);
      m_RobustSSM->GetSurfaceSimilarityParameters(shape, weights, matrix, b);
      vtkSmartPointer<vtkPolyData> shape2 = m_RobustSSM->GetSurface(matrix, b, false);
      double epsiSqr = m_RobustSSM->PolyDataDistance(shape, shape2)/shape->GetNumberOfPoints();

      table->SetValue(i, 0, i);
      table->SetValue(i, 1, epsiSqr);

      generalisability += epsiSqr;
      m_RobustSSM->AddShape(shape, weights);

      qApp->processEvents();
    }
    printInfo("Generalisability Score: " + QString::number(generalisability));
    //table->Update();
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

void milxQtRobustShapeModel::eigenvalues()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> values = m_RobustSSM->GetPCA()->GetEvals();

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
    //table->Update();
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

void milxQtRobustShapeModel::eigenmodes()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < m_RobustSSM->GetNumberOfModes(); j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < m_RobustSSM->GetNumberOfShapes(); k ++)
            column->InsertNextValue(0.0);
        const QString nameOfColumn = "Mode " + QString::number(j+1);
        column->SetName(nameOfColumn.toStdString().c_str());
        table->AddColumn(column);
    }
    printDebug("Modes Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    for(int j = 0; j < m_RobustSSM->GetNumberOfShapes(); j ++)
    {
        vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
            b->SetNumberOfComponents(1);
            b->SetNumberOfTuples(m_RobustSSM->GetNumberOfModes());

        vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
        m_RobustSSM->GetShape(j)->GetPointData()->SetScalars(m_RobustSSM->GetWeights(j));
        m_RobustSSM->GetSurfaceSimilarityParameters(m_RobustSSM->GetShape(j), m_RobustSSM->GetWeights(j), matrix, b);

        for(int k = 0; k < m_RobustSSM->GetNumberOfModes(); k ++)
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

void milxQtRobustShapeModel::parameters()
{
    if(!m_modelled)
        generateSSM();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    ///Init table
    for(int j = 0; j < m_RobustSSM->GetNumberOfModes(); j ++)
    {
        vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
        for(int k = 0; k < m_RobustSSM->GetNumberOfShapes(); k ++)
            column->InsertNextValue(0.0);
        const QString nameOfColumn = "Mode " + QString::number(j+1);
        column->SetName(nameOfColumn.toStdString().c_str());
        table->AddColumn(column);
    }
    printDebug("Parameters Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));

    for(int j = 0; j < m_RobustSSM->GetNumberOfShapes(); j ++)
    {
        vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
            b->SetNumberOfComponents(1);
            b->SetNumberOfTuples(m_RobustSSM->GetNumberOfModes());

        vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
        m_RobustSSM->GetShape(j)->GetPointData()->SetScalars(m_RobustSSM->GetWeights(j));
        m_RobustSSM->GetSurfaceSimilarityParameters(m_RobustSSM->GetShape(j), m_RobustSSM->GetWeights(j), matrix, b);

        for(int k = 0; k < m_RobustSSM->GetNumberOfModes(); k ++)
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

void milxQtRobustShapeModel::reset()
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

void milxQtRobustShapeModel::createActions()
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection
    contextMenu->setTitle("Robust Shape Modelling");

    actionMean = new QAction(this);
        actionMean->setText("&Mean Model");
        actionMean->setShortcut(tr("Alt+m"));
        actionMean->setCheckable(true);
        actionMean->setChecked(false);
    actionAligned = new QAction(this);
        actionAligned->setText("&Aligned Points");
        actionAligned->setShortcut(tr("Alt+a"));
        actionAligned->setCheckable(true);
        actionAligned->setChecked(false);
    actionOriginal = new QAction(this);
        actionOriginal->setText("&Original Points");
        actionOriginal->setShortcut(tr("Alt+o"));
        actionOriginal->setCheckable(true);
        actionOriginal->setChecked(false);
    actionModesAsVectors = new QAction(this);
        actionModesAsVectors->setText("Modes as &Vectors");
        actionModesAsVectors->setShortcut(tr("Alt+v"));
        actionModesAsVectors->setCheckable(true);
        actionModesAsVectors->setChecked(false);
    actionModesAsTensors = new QAction(this);
        actionModesAsTensors->setText("Modes as &Tensors");
        actionModesAsTensors->setShortcut(tr("Alt+t"));
        actionModesAsTensors->setCheckable(true);
        actionModesAsTensors->setChecked(false);
    actionModesAsCollection = new QAction(this);
        actionModesAsCollection->setText("Modes as Surface &Collection");
        actionModesAsCollection->setShortcut(tr("Alt+c"));
    actionCorrespond = new QAction(this);
        actionCorrespond->setText("Correspondences as a &Hedgehog");
        actionCorrespond->setShortcut(tr("Alt+h"));
        actionCorrespond->setCheckable(true);
        actionCorrespond->setChecked(false);

    //plots menu
    actionCompact = new QAction(this);
        actionCompact->setText("Compactness");
        actionCompact->setShortcut(tr("Shift+Alt+c"));
    actionSpecificity = new QAction(this);
        actionSpecificity->setText("Specificity");
        actionSpecificity->setShortcut(tr("Shift+Alt+s"));
    actionGeneralise = new QAction(this);
        actionGeneralise->setText("Generalisability");
        actionGeneralise->setShortcut(tr("Shift+Alt+g"));
    actionValues = new QAction(this);
        actionValues->setText("Eigenvalues");
        actionValues->setShortcut(tr("Shift+Alt+v"));
    actionModes = new QAction(this);
        actionModes->setText("Primary Eigenmodes");
        actionModes->setShortcut(tr("Shift+Alt+r"));
    actionParameters = new QAction(this);
        actionParameters->setText("Training Shape Parameters");
        actionParameters->setShortcut(tr("Shift+Alt+p"));

    //Procrustes menu
    actionRigid = new QAction(this);
        actionRigid->setText("Rigid");
        actionRigid->setShortcut(tr("Ctrl+Alt+r"));
        actionRigid->setCheckable(true);
    actionSimilarity = new QAction(this);
        actionSimilarity->setText("Similarity");
        actionSimilarity->setShortcut(tr("Ctrl+Alt+s"));
        actionSimilarity->setCheckable(true);
        actionSimilarity->setChecked(true);
    actionAffine = new QAction(this);
        actionAffine->setText("Affine");
        actionAffine->setShortcut(tr("Ctrl+Alt+a"));
        actionAffine->setCheckable(true);
    alignGroup = new QActionGroup(this);
        alignGroup->addAction(actionRigid);
        alignGroup->addAction(actionSimilarity);
        alignGroup->addAction(actionAffine);

    actionAlignment = new QAction(this);
        actionAlignment->setText("Output &Alignment as Images");
        actionAlignment->setShortcut(tr("Alt+r"));
    actionOriginalMeshes = new QAction(this);
        actionOriginalMeshes->setText("Output &Original Meshes");
        actionOriginalMeshes->setShortcut(tr("Shift+Alt+o"));
    actionAlignedMeshes = new QAction(this);
        actionAlignedMeshes->setText("Output &Aligned Meshes");
        actionAlignedMeshes->setShortcut(tr("Shift+Alt+a"));
    actionPointIDs = new QAction(this);
        actionPointIDs->setText("Output &Point IDs as Images");
        actionPointIDs->setShortcut(tr("Alt+p"));
    actionCoordinates = new QAction(this);
        actionCoordinates->setText("Output &Coordinates as Images");
        actionCoordinates->setShortcut(tr("Alt+c"));

    actionReplaceOriginal = new QAction(this);
        actionReplaceOriginal->setText("Replace Originals with Aligned");
        actionReplaceOriginal->setShortcut(tr("Shift+Alt+r"));
}

void milxQtRobustShapeModel::contextMenuEvent(QContextMenuEvent *contextEvent)
{
    createMenu(contextMenu);

    contextMenu->exec(contextEvent->globalPos());
}
