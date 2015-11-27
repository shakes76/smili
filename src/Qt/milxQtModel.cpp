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
#include "milxQtModel.h"

//Qt
#include <QInputDialog>
#include <QColorDialog>
//Data
#include <vtkMath.h>
#include <vtkTable.h>
//Mapper
#include <vtkGraphicsFactory.h>
#include <vtkCamera.h>
//Graphing
#include <vtkDoubleArray.h>
#include <vtkDataSetMapper.h>
#include <vtkIdFilter.h>
#include <vtkSelectVisiblePoints.h>
//Contouring
#include <vtkGlyph3D.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkPolygonalSurfacePointPlacer.h>
#include <vtkPolygonalSurfaceContourLineInterpolator.h>
//Filters
#include <vtkTextProperty.h>
#include <vtkCoordinate.h> //for scalarbar
#include <vtkImageQuantizeRGBToIndex.h>
#include <vtkImageToPolyDataFilter.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
//#include <vtkStreamLine.h>
#include <vtkStreamTracer.h>
#include <vtkLogLookupTable.h>
//SMILI
#include "milxMath.h"
#include "milxColourMap.h"

#include "milxQtFile.h"

///milxQtModel
milxQtModel::milxQtModel(QWidget *theParent, bool contextSystem) : milxQtRenderWindow(theParent, contextSystem)
{
    ///Initialise variables
    modelled = false;
    resetFlags();

    colourRed = defaultColour;
    colourGreen = defaultColour;
    colourBlue = defaultColour;

    ///Set internals to zero
    modelCentroid.fill(0.0);
    modelCovarianceMatrix.set_size(3, 3);
    modelCovarianceMatrix.fill(0.0);

    ///Connect Model Object for progress updates
    linkProgressEventOf(milx::VTKProgressUpdates->GetUpdateObject()); //link internal observer

    ///Set strings
    milxQtWindow::prefix = "Mdl: ";

    ///Allocate critical aspects
    milxQtWindow::setDeletableOnClose(true);

    //View defaults
    milxQtRenderWindow::setDefaultView(AXIAL);
    setView(AXIAL);
    setDefaultOrientation(RADIOLOGICAL);

    createActions();
    createConnections();

    surfaceAct->setChecked(true);
}

milxQtModel::~milxQtModel()
{
    ///Smart pointers handle deletion
}

//Sets
void milxQtModel::AddInput(vtkSmartPointer<vtkPolyData> mesh)
{
    printInfo("Appending Model data");

    model.Append(mesh);

    resetFlags();
    loaded = true;
    appended = true;
}

void milxQtModel::AddArray(vtkSmartPointer<vtkDataArray> array)
{
    model.AddArray(array);
    resetFlags();

    loaded = true;
}

void milxQtModel::SetInput(vtkSmartPointer<vtkPolyData> mesh)
{
    model.SetInput(mesh);
    resetFlags();
    loaded = true;
}
/*
void milxQtModel::SetInputPointSet(vtkSmartPointer<vtkPointSet> mesh)
{
    if(!modelled)
    {
        modelMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        modelActor = vtkSmartPointer<vtkLODActor>::New();
    }

    printDebug("Generating Model for Point Set.");
    if(largeMode)
        modelMapper->GlobalImmediateModeRenderingOn();
#if VTK_MAJOR_VERSION <= 5
    modelMapper->SetInputConnection(mesh->GetProducerPort());
#else
    modelMapper->SetInputData(mesh);
#endif
    modelMapper->ScalarVisibilityOn();
    modelMapper->SetScalarModeToUsePointData();
    modelMapper->SetColorModeToMapScalars();
    modelMapper->Update();

    ///Setup actor for rendering
    modelActor->SetMapper(modelMapper);
    modelActor->GetProperty()->SetColor(colourRed, colourGreen, colourBlue);
    modelActor->GetProperty()->SetSpecularColor(1, 1, 1);
    modelActor->GetProperty()->SetSpecular(0.25);
    modelActor->GetProperty()->SetSpecularPower(10);
    modelActor->GetProperty()->SetAmbient(0.2);
    modelActor->GetProperty()->SetDiffuse(0.8);

    modelled = true;
    loaded = true;
    milxQtRenderWindow::AddActor(modelActor);
    milxQtRenderWindow::generateRender();
}*/

//void milxQtModel::SetNumberOfInputs(const int inputs)
//{
//    if(!appendPolyData)
//        appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
//    appendPolyData->SetNumberOfInputs(inputs);
//}

//void milxQtModel::SetInput(const int index, vtkSmartPointer<vtkPolyData> mesh)
//{
//    if(!appendPolyData)
//        appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
//    if(loaded && !appended)
//        appendPolyData->AddInput(currentMesh);
//    appendPolyData->SetInput(index, mesh);
//    linkProgressEventOf(appendPolyData); //Keeps UI responsive
//    currentMesh = appendPolyData->GetOutput();
//    resetFlags();
//    loaded = true;
//    appended = true;
//}

void milxQtModel::SetPoints(vtkSmartPointer<vtkPoints> modelPoints)
{
    model.SetPoints(modelPoints);
    resetFlags();
    loaded = true;
}

void milxQtModel::SetPolys(vtkSmartPointer<vtkCellArray> modelPolys)
{
    model.SetPolys(modelPolys);
    resetFlags();
    loaded = true;
}

void milxQtModel::SetScalars(vtkSmartPointer<vtkDataArray> modelScalars)
{
    model.SetScalars(modelScalars);
    resetFlags();

    if(modelScalars)
		    scalarsSet = true;
    loaded = true;
}

void milxQtModel::RemoveScalars()
{
    model.RemoveScalars();
}

void milxQtModel::SetVectors(vtkSmartPointer<vtkDataArray> modelVectors)
{
    model.SetVectors(modelVectors);
    resetFlags();

    if(modelVectors)
		    scalarsSet = true;
    loaded = true;
}

void milxQtModel::SetTransform(vtkSmartPointer<vtkTransform> newTransform)
{
    if(loaded)
    {
        emit working(-1);
        ///Transform this model
        model.SetTransform(newTransform);

        if(normalsModel)
            normalsModel->SetTransform(newTransform);
        if(centroidModel)
            centroidModel->SetTransform(newTransform);
        if(outlineModel)
            outlineModel->SetTransform(newTransform);
        emit done(-1);

        resetFlags();
    }
}

void milxQtModel::SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph)
{
    model.SetGraph(graph);

    resetFlags();
    loaded = true;
}

//Gets
vtkSmartPointer<vtkPolyData> milxQtModel::GetPolyDataInput()
{
    return model.PreviousResult();
}

vtkSmartPointer<vtkPolyData> milxQtModel::GetOutput()
{
    return model.Result();
}

//Operators
milxQtModel& milxQtModel::operator=(const milxQtModel &operand)
{
    if(this != &operand)
    {
        ///Copy Flags' States
        loaded = operand.loaded;
        computedCentroid = operand.computedCentroid;
        computedCovariance = operand.computedCovariance;

        ///Copy internal computed variables
        modelCentroid = operand.modelCentroid;
        modelCovarianceMatrix = operand.modelCovarianceMatrix;
    }

    return *this;
}

void milxQtModel::refresh()
{
    if(!modelled)
        reset();
    else
    {
        model.Update();
    #if VTK_MAJOR_VERSION <=5
        modelMapper->SetInput(model.Result());
    #else
        modelMapper->SetInputData(model.Result());
    #endif
        modelMapper->Update();
        milxQtRenderWindow::refresh();
    }
}

coordinate& milxQtModel::centroid()
{
    if(!computedCentroid)
    {
        modelCentroid = milx::Math<coordinateType>::Centroid(model.Result()->GetPoints());
        computedCentroid = true;
    }

    if(centroidAct->isChecked())
    {
        printInfo("Displaying Centroid: (" + QString::number(modelCentroid[0]) + ", " + QString::number(modelCentroid[1]) + ", " + QString::number(modelCentroid[2]) + ")");
        if(!centroidModel)
            centroidModel = new milxQtModel;

        vtkPoints *centroidPoint = vtkPoints::New();
        centroidPoint->InsertNextPoint(modelCentroid.data_block());
        centroidModel->SetPoints(centroidPoint);
        centroidModel->generatePointModel(2.0, 1.0, 0.0, 0.0);

        milxQtRenderWindow::AddActor(centroidModel->GetActor());
        milxQtRenderWindow::generateRender();
    }
    else
    {
        if(centroidModel)
            milxQtRenderWindow::RemoveActor(centroidModel->GetActor());
    }

    return modelCentroid;
}

double milxQtModel::centroidSize(bool average)
{
    coordinate currentCentroid = centroid();
    return milx::Math<coordinateType>::CentroidSize(model.Result()->GetPoints(), currentCentroid, average);
}

vnl_matrix<double>& milxQtModel::covarianceMatrix()
{
    if(!computedCovariance)
    {
        coordinate currentCentroid = centroid();
        modelCovarianceMatrix = milx::Math<coordinateType>::CovarianceMatrix(model.Result()->GetPoints(), currentCentroid);
        computedCovariance = true;
    }

    return modelCovarianceMatrix;
}

void milxQtModel::generateDelaunayGraph(float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        model.DelaunayGraph3D();
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateDelaunay2DTriangulation()
{
    if(loaded)
    {
        emit working(-1);
        model.Delaunay2DTriangulation();

        vtkSmartPointer<vtkPolyDataMapper> delaunayMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        #if VTK_MAJOR_VERSION <=5
            delaunayMapper->SetInput(model.Result());
        #else
            delaunayMapper->SetInputData(model.Result());
        #endif

        if(!delaunayActor)
            delaunayActor = vtkSmartPointer<vtkLODActor>::New();
            delaunayActor->SetMapper(delaunayMapper);

        milxQtRenderWindow::AddActor(delaunayActor);
        milxQtRenderWindow::generateRender();
        emit done(-1);
    }
}

void milxQtModel::generateDelaunayTriangulation(float red, float green, float blue)
{
    if(loaded)
    {
        if(!genDelaunayTriAct->isChecked())
            milxQtRenderWindow::RemoveActor(delaunayActor);
        else
        {
            emit working(-1);
            ///Generate Delaunay Graph
            model.DelaunayTriangulation();

            vtkSmartPointer<vtkDataSetMapper> delaunayMapper = vtkSmartPointer<vtkDataSetMapper>::New();
            #if VTK_MAJOR_VERSION <=5
                delaunayMapper->SetInput(model.Result());
            #else
                delaunayMapper->SetInputData(model.Result());
            #endif

            if(!delaunayActor)
                delaunayActor = vtkSmartPointer<vtkLODActor>::New();
                delaunayActor->SetMapper(delaunayMapper);
    //            delaunayActor->GetProperty()->SetColor(1,0,0);

            milxQtRenderWindow::AddActor(delaunayActor);
            milxQtRenderWindow::generateRender();
            emit done(-1);
        }
    }
}

void milxQtModel::generateWireframe(float red, float green, float blue)
{
    if(!loaded)
        return;
    if(!modelled)
        generateModel(red, green, blue);

    modelActor->GetProperty()->SetRepresentationToWireframe();

    wireframeAct->setChecked(true);
}

void milxQtModel::generateVertices(float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateVertices();
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateVerticesAs(const GlyphType glyphType, float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateVerticesAs2DGlyphs(glyphType);
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateTubes(float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateTubes();
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generatePoints(float red, float green, float blue)
{
    if(!loaded)
        return;
    if(!modelled)
        generateModel(red, green, blue);

    modelActor->GetProperty()->SetRepresentationToPoints();

    pointsAct->setChecked(true);
}

void milxQtModel::generatePointModel(double newScale, float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        double bounds[6];
        model.Result()->GetBounds(bounds);
        double scaling = ( (bounds[1]-bounds[0]) + (bounds[3]-bounds[2]) + (bounds[5]-bounds[4]) )/150;

        if(newScale == 1) //default scaling, so do automatically
            scaling *= newScale;
        else
            scaling = newScale;

        model.GeneratePointModel(scaling);
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateSampledPoints(float distance, float red, float green, float blue)
{
    if(loaded)
    {
        bool ok = false;

        if(distance <= 0)
        {
            distance = QInputDialog::getDouble(this, tr("Please Provide the sample distance"),
                                             tr("Distance:"), 0.1, 0.0, 2147483647, 7, &ok);
        }
        else
            ok = true;

        if(!ok)
            return;

        emit working(-1);
        model.GenerateSampledPoints(distance);
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateVectorField(double newScale, float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        bool useNormals = false;
        vtkSmartPointer<vtkFloatArray> vectorArray;
        if(!GetVectors()) //use scalars instead
        {
            printInfo("No Vectors found.");
            if(!GetNormals())
            {
                printInfo("Generating Vector Field from Scalars (with normals).");
                model.GenerateNormals();
            }

            ///Use Glyphs (vtkGlyph3D) at each point in the set with normals
            vectorArray = vtkFloatArray::SafeDownCast(model.GetNormals());
            useNormals = true;
        }
        else
        {
            printInfo("Generating Vector Field from Vectors.");
            vectorArray = vtkFloatArray::SafeDownCast(model.GetVectors());
        }

        if(newScale == 0.0 && vectorArray) //auto rescale
        {
            coordinate meanDirection(0.0); //used to set scale
            for(vtkIdType j = 0; j < model.GetNumberOfPoints(); j ++)
            {
                coordinate vectorData(vectorArray->GetTuple3(j));

                meanDirection += vectorData;
            }
            meanDirection /= model.GetNumberOfPoints();
            printDebug("Mean Vector: " + QString::number(meanDirection[0]) + ", " + QString::number(meanDirection[1]) + ", " + QString::number(meanDirection[2]) + ", ");
            printDebug("Mean Vector L2 Norm: " + QString::number(meanDirection.two_norm()));

            bool ok1 = false;
            double newScaling = 1.0/meanDirection.two_norm();
            newScaling = QInputDialog::getDouble(this, tr("Please enter scaling for the vector field"),
                                tr("Scaling:"), newScaling, -2147483647, 2147483647, 7, &ok1);

            if(!ok1)
                newScale = 0.0;
            else
                newScale = newScaling;
            printInfo("Vectors will be scaled by " + QString::number(newScale) + " for display.");
        }

        ///Use Glyphs (vtkGlyph3D) at each point in the set
        model.GenerateVectorField(newScale, useNormals);
        emit done(-1);
        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateTensorField(double newScale, float red, float green, float blue)
{
    if(!loaded)
        return;

    if(!GetTensors())
    {
        printError("No tensors found in polydata. Exiting.");
        return;
    }

    emit working(-1);
    if(newScale == 0.0) //auto rescale
    {
        vnl_matrix<double> meanDirection(3, 3, 0.0); //used to set scale
        for(vtkIdType j = 0; j < model.GetNumberOfPoints(); j ++)
        {
            vnl_matrix<double> vectorData(GetTensors()->GetTuple9(j), 3 ,3);

            meanDirection += vectorData;
        }
        meanDirection /= model.GetNumberOfPoints();
        printDebug("Mean Vector L2 Norm: " + QString::number(meanDirection.frobenius_norm()));

        bool ok1 = false;
        double newScaling = 1.0/meanDirection.frobenius_norm();
        newScaling = QInputDialog::getDouble(this, tr("Please enter scaling for the vector field"),
                            tr("Scaling:"), newScaling, -2147483647, 2147483647, 7, &ok1);

        if(!ok1)
            newScale = 0.0;
        else
            newScale = newScaling;
        printInfo("Vectors will be scaled by " + QString::number(newScale) + " for display.");
    }

    printInfo("Generating Tensor Field.");
    model.GenerateTensorField(newScale);
    emit done(-1);

    ///Generate the model
    generateModel(red, green, blue);
}

void milxQtModel::generateStreamLines(vtkSmartPointer<vtkImageData> vectorFieldData, float timestep)
{
    if(!loaded)
        return;

    if(!GetPoints())
    {
        printError("No points found in polydata. Exiting.");
        return;
    }

    bool ok = false;
    if(timestep == 0.0)
    {
      timestep = QInputDialog::getDouble(this, tr("Please Provide the timestep"),
                                         tr("Timestep:"), 0.1, 0.0, 100.0, 5, &ok);
    }

    if(!ok)
        return;

    emit working(-1);

    double bounds[6];
    vectorFieldData->GetBounds(bounds);

    const double maxProp = bounds[1]-bounds[0];
    const double maxSteps = maxProp/timestep;

    ///Not that in order for the stream lines to work, one requires both cell and point data.
    ///For the seeds one only needs point data.
    printInfo("Generating Streamlines.");
    // Streamline itself
    vtkSmartPointer<vtkStreamTracer> streamLines = vtkSmartPointer<vtkStreamTracer>::New();
    #if VTK_MAJOR_VERSION <= 5
        streamLines->SetInput(vectorFieldData);
        streamLines->SetSource(model.Result());
    #else
        streamLines->SetInputData(vectorFieldData);
        streamLines->SetSourceData(model.Result());
    #endif
//        streamLines->SetMaximumPropagationTime(200);
        streamLines->SetMaximumPropagation(maxProp);
        streamLines->SetMaximumNumberOfSteps(maxSteps);
        streamLines->SetInitialIntegrationStep(timestep);
//        streamLines->SetIntegrationDirectionToForward();
        streamLines->SetIntegrationDirectionToBoth();
        streamLines->SetIntegratorTypeToRungeKutta4();
//        streamLines->SetInterpolatorTypeToDataSetPointLocator();
//        streamLines->SetInterpolatorTypeToCellLocator();
//        streamLines->SetIntegrationStepUnit(1); //LENGTH_UNIT
//        streamLines->VorticityOn();
        linkProgressEventOf(streamLines);
        streamLines->Update();
    emit done(-1);

    model.SetInput(streamLines->GetOutput());

    ///Generate the model
    generateModel();
}

void milxQtModel::generateHedgehog(double newScale, float red, float green, float blue)
{
    if(loaded)
    {
        emit working(-1);
        printInfo("Generated Hedgehog.");
        model.GenerateHedgehog(newScale);
        emit done(-1);

        ///Generate the model
        generateModel(red, green, blue);
    }
}

void milxQtModel::generateSurface(float red, float green, float blue)
{
    if(!loaded)
        return;
    if(!modelled)
        generateModel(red, green, blue);

    modelActor->GetProperty()->SetRepresentationToSurface();

    surfaceAct->setChecked(true);
}

void milxQtModel::generateNormals(int pointNormals)
{
    if(!loaded)
        return;

    emit working(-1);
    printInfo("Generating Normals for the model");
    model.GenerateNormals(pointNormals);
    emit done(-1);
}

void milxQtModel::generateModel(float red, float green, float blue)
{
    if(loaded)
    {
        if(!modelled)
        {
            modelMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            modelActor = vtkSmartPointer<vtkLODActor>::New();
        }

        if(model.Result()->GetPointData()->GetScalars())
            scalarsSet = true;

        printDebug("Generating Model.");
        if(largeMode)
            modelMapper->GlobalImmediateModeRenderingOn();
    #if VTK_MAJOR_VERSION <=5
        modelMapper->SetInput(model.Result());
    #else
        modelMapper->SetInputData(model.Result());
    #endif // VTK_MAJOR_VERSION
        modelMapper->ScalarVisibilityOn();
        modelMapper->SetScalarModeToUsePointData();
//        modelMapper->SetColorModeToMapScalars(); //force colour unsigned char scalars with this option
        modelMapper->SetColorModeToDefault();
        if(scalarsSet && !customScalarBar)
        {
            vtkSmartPointer<vtkLookupTable> tmpLookupTable = vtkLookupTable::SafeDownCast(lookupTable);

            tmpLookupTable->SetTableRange(model.Result()->GetScalarRange());
            tmpLookupTable->Build();

            updateLookupTable();
        }
        modelMapper->Update();

        colourRed = red;
        colourGreen = green;
        colourBlue = blue;

        ///Setup actor for rendering
        modelActor->SetMapper(modelMapper);
        modelActor->GetProperty()->SetColor(colourRed, colourGreen, colourBlue);
        modelActor->GetProperty()->SetSpecularColor(1, 1, 1);
        modelActor->GetProperty()->SetSpecular(0.25);
        modelActor->GetProperty()->SetSpecularPower(10);
        modelActor->GetProperty()->SetAmbient(0.2);
        modelActor->GetProperty()->SetDiffuse(0.8);

        modelled = true; //dont move down, ordering necessary
        milxQtRenderWindow::AddActor(modelActor);
        if(!rendered) ///keep view that may be set by user
        {
            milxQtRenderWindow::generateRender();
            if(milxQtRenderWindow::useDefaultView)
            	milxQtRenderWindow::setView(milxQtRenderWindow::defaultView); //Default view
        }
        milxQtRenderWindow::refresh();
    }
}

void milxQtModel::generateLabels()
{
    if(loaded)
    {
        emit working(-1);
        if(!labelled)
            modelLabelsActor = vtkSmartPointer<vtkActor2D>::New();

        vtkSmartPointer<vtkIdFilter> ids = vtkSmartPointer<vtkIdFilter>::New();
    #if VTK_MAJOR_VERSION <=5
        ids->SetInput(model.Result());
    #else
        ids->SetInputData(model.Result());
    #endif // VTK_MAJOR_VERSION
        ids->PointIdsOn();

        vtkSmartPointer<vtkSelectVisiblePoints> visPts = vtkSmartPointer<vtkSelectVisiblePoints>::New();
        visPts->SetInputConnection( ids->GetOutputPort() );
        visPts->SetRenderer( renderer );

        vtkSmartPointer<vtkLabeledDataMapper> labeller = vtkSmartPointer<vtkLabeledDataMapper>::New();
        labeller->SetInputConnection( visPts->GetOutputPort() );
        linkProgressEventOf(labeller);
        labeller->SetLabelModeToLabelFieldData();

        modelLabelsActor->SetMapper( labeller );
        emit done(-1);
        labelled = true;

        milxQtRenderWindow::AddActor2D(modelLabelsActor);
        milxQtRenderWindow::generateRender();
    }
}

void milxQtModel::generatePointIDsScalars()
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateVertexScalars();
        emit done(-1);

        ///Generate the model
        generateModel();
    }
}

void milxQtModel::generateIsoSurface(vtkSmartPointer<vtkImageData> img, int contourNumber, double value)
{
    bool ok1 = false, ok2 = false;

    if(contourNumber < 0)
    {
//      contourNumber = QInputDialog::getInteger(this, tr("Please enter minimum iso surface value"),
//                        tr("Contour Number:"), 0, -2147483647, 2147483647, 1, &ok1);
        contourNumber = 0;
        ok1 = true;
        value = QInputDialog::getDouble(this, tr("Please enter label/iso surface value"),
                          tr("Value:"), 0.5, -2147483647, 2147483647, 7, &ok2);
    }
    else
        ok1 = ok2 = true;

    if(ok1 && ok2)
    {
        emit working(-1);
        vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
        #if VTK_MAJOR_VERSION <=5
            imageReorient->SetInput(img);
        #else
            imageReorient->SetInputData(img);
        #endif // VTK_MAJOR_VERSION
            imageReorient->SetFilteredAxis(1);
            imageReorient->FlipAboutOriginOn();
            linkProgressEventOf(imageReorient);
            imageReorient->Update(); //ITK image would have been flipped

        model.IsoSurface(imageReorient->GetOutput(), value);
        setName("Iso Surface");
        loaded = true;
        emit done(-1);

        refresh();
    }
    else
        loaded = false;
}

void milxQtModel::generatePolyDataFromImage(vtkSmartPointer<vtkImageData> img)
{
    emit working(-1);
    vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    #if VTK_MAJOR_VERSION <=5
        imageReorient->SetInput(img);
    #else
        imageReorient->SetInputData(img);
    #endif
        imageReorient->SetFilteredAxis(1);
        imageReorient->FlipAboutOriginOn();
        imageReorient->Update(); //ITK image would have been flipped

//    vtkSmartPointer<vtkImageQuantizeRGBToIndex> quant = vtkSmartPointer<vtkImageQuantizeRGBToIndex>::New();
//        quant->SetInput(img);
//        quant->SetNumberOfColors(16);

    vtkSmartPointer<vtkImageToPolyDataFilter> imgToModel = vtkSmartPointer<vtkImageToPolyDataFilter>::New();
//        imgToModel->SetInputConnection(quant->GetOutputPort());
#if VTK_MAJOR_VERSION <=5
    imgToModel->SetInput(imageReorient->GetOutput());
#else
    imgToModel->SetInputData(imageReorient->GetOutput());
#endif // VTK_MAJOR_VERSION
//        imgToModel->SetLookupTable(quant->GetLookupTable());
//        imgToModel->SetColorModeToLUT();
    imgToModel->SetOutputStyleToPolygonalize();
//        imgToModel->SetError(0);
    imgToModel->DecimationOn();
//        imgToModel->SetDecimationError(0.0);
//        imgToModel->SetSubImageSize(25);
    linkProgressEventOf(imgToModel); //Keeps UI responsive
    imgToModel->Update();

    model.SetInput(imgToModel->GetOutput());
    loaded = true;

    triangulate();
    emit done(-1);

    refresh();
}

void milxQtModel::generateKMeansClustering(int numberOfClusters)
{
    if(loaded)
    {
        bool ok = false;

        if(numberOfClusters <= 0)
        {
            numberOfClusters = QInputDialog::getInt(this, tr("Please Provide the number of clusters to use"),
                                              tr("Clusters:"), 2, 0, 2147483647, 1, &ok);
        }
        else
            ok = true;

        if(!ok)
            return;

        emit working(-1);
        model.GenerateKMeansClustering(numberOfClusters);
        emit done(-1);

        refresh();
    }
}

void milxQtModel::generateQuantisedPoints(float quantiseFactor)
{
    if(loaded)
    {
        bool ok = false;

        if(quantiseFactor <= 0)
        {
            quantiseFactor = QInputDialog::getDouble(this, tr("Please Provide the Quantise Factor"),
                                             tr("Factor:"), 0.25, 0.0, 2147483647, 5, &ok);
        }
        else
            ok = true;

        if(!ok)
            return;

        emit working(-1);
        model.GenerateQuantisedPoints(quantiseFactor);
        emit done(-1);

        refresh();
    }
}

void milxQtModel::generateCappedBoundaries()
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateCappedBoundaries();
        emit done(-1);

        refresh();
    }
}

void milxQtModel::generateRegionLabels()
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateRegions();
        emit done(-1);

        refresh();
    }
}

void milxQtModel::generateElevation()
{
    if(loaded)
    {
        emit working(-1);
        model.GenerateElevation();
        emit done(-1);

        refresh();
    }
}

//~ void milxQtModel::generateReebGraphs()
//~ {
    //~ if(loaded)
    //~ {
        //~ emit working(-1);
        //~ model.GenerateReebGraph();
        //~ emit done(-1);

        //~ refresh();
    //~ }
//~ }

void milxQtModel::undoProcessing()
{
    if(!loaded || !model.PreviousResult())
        return;

    //Switch current and previous mesh pointers
    milx::Swap< vtkSmartPointer<vtkPolyData> >(model.Result(), model.PreviousResult());

    if(!undidProcess)
    {
        doAct->setText(QApplication::translate("Model", "&Redo Last Processing", 0, QApplication::UnicodeUTF8));
        undidProcess = true;
    }
    else
    {
        doAct->setText(QApplication::translate("Model", "&Undo Last Processing", 0, QApplication::UnicodeUTF8));
        undidProcess = false;
    }

    refresh();
    modelInfo();
}

void milxQtModel::clean()
{
    if(!loaded)
        return;

    emit working(-1);
    model.Clean();
    refresh();
    emit done(-1);
}

void milxQtModel::triangulate()
{
    if(!loaded)
        return;

    emit working(-1);
    model.Triangulate();
    refresh();
    emit done(-1);
}

void milxQtModel::decimate(double factor)
{
    if(!loaded)
        return;

    bool ok = false;

    if(factor == 0.0)
    {
        factor = QInputDialog::getDouble(this, tr("Please Provide the reduction factor"),
                                         tr("Factor:"), 0.5, 0.0, 1.0, 3, &ok);
    }
    else
        ok = true;

    if( factor > 0.0 && ok )
    {
        emit working(-1);
        printInfo("There were " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
        printInfo("There were " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");
        printInfo("--- Decimating with using DecimatePro");

        triangulate();
        model.Decimate(factor);

        printInfo("There are " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
        printInfo("There are " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");

        refresh();
        emit done(-1);
    }
}

void milxQtModel::quadricDecimate(double factor)
{
    if(!loaded)
        return;

    bool ok = false;

    if(factor == 0.0)
    {
        factor = QInputDialog::getDouble(this, tr("Please Provide the target reduction factor"),
                                         tr("Target Factor:"), 0.5, 0.0, 1.0, 3, &ok);
    }
    else
        ok = true;

    if( factor > 0.0 && ok )
    {
        emit working(-1);
        printInfo("There were " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
        printInfo("There were " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");
        printInfo("--- Decimating with using Quadric Decimation");

        model.Triangulate();
        model.QuadricDecimate(factor);

        printInfo("There are " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
        printInfo("There are " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");

        refresh();
        emit done(-1);
    }
}

void milxQtModel::clusterDecimate()
{
    if(!loaded)
        return;

    emit working(-1);
    printInfo("There were " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
    printInfo("There were " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");
    printInfo("--- Decimating with using Quadric Clustering Decimation");

    triangulate();
    model.ClusterDecimate();

    printInfo("There are " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
    printInfo("There are " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");

    refresh();
    emit done(-1);
}

void milxQtModel::renameScalars(QString newName)
{
    if(!loaded)
        return;

    if(newName.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        newName = QInputDialog::getText(this, tr("New Name for Scalars"),
                                          tr("New name:"), QLineEdit::Normal,
                                          "Scalars");
    }

    if(newName.isEmpty())
        return;

    model.Result()->GetPointData()->GetScalars()->SetName(newName.toStdString().c_str());
    ///The main reason you would rename scalars array is to label the scalar bar.
    if(scaleBefore)
        scale->SetTitle(newName.toStdString().c_str());

    refresh();
}

void milxQtModel::threshold(double lowValue, double upValue)
{
    if(!loaded)
        return;

    bool ok = false;

    if(upValue == 0.0 && lowValue == 0.0)
    {
        double range[2];
        model.Result()->GetScalarRange(range);

        lowValue = QInputDialog::getDouble(this, tr("Lower level of band to keep"),
                                         tr("Lower Value (inclusive):"), range[0], -2147483647, 2147483647, 12, &ok);
        upValue = QInputDialog::getDouble(this, tr("Upper level of band to keep (inclusive)"),
                                         tr("Upper Value (inclusive):"), range[1], -2147483647, 2147483647, 12, &ok);
    }
    else
        ok = true;

    if(ok)
    {
        emit working(-1);
        ///Threshold based on scalars
        printInfo("Thresholding Scalars and clipping mesh");
        model.Clip(lowValue, upValue);

        reset();
        emit done(-1);
    }
}

void milxQtModel::thresholdScalars(double lowValue, double upValue, double outsideVal)
{
    if(!loaded)
        return;

    bool ok = false;

    if(upValue == 0.0 && lowValue == 0.0 && outsideVal == 0.0)
    {
        double range[2];
        model.Result()->GetScalarRange(range);

        lowValue = QInputDialog::getDouble(this, tr("Lower level of band to keep"),
                                         tr("Lower Value (inclusive):"), range[0], -2147483647, 2147483647, 12, &ok);
        upValue = QInputDialog::getDouble(this, tr("Upper level of band to keep (inclusive)"),
                                         tr("Upper Value (inclusive):"), range[1], -2147483647, 2147483647, 12, &ok);
        outsideVal = QInputDialog::getDouble(this, tr("Value to give scalars outside the range"),
                                         tr("Outside Value:"), 0.0, -2147483647, 2147483647, 12, &ok);
    }
    else
        ok = true;

    if(ok)
    {
        emit working(-1);
        printInfo("Thresholding Scalars");
        model.Threshold(lowValue, upValue, outsideVal);
        reset();
        emit done(-1);
    }
}

void milxQtModel::thresholdScalarsBinary(double lowValue, double upValue, double insideVal, double outsideVal)
{
    if(!loaded)
        return;

    bool ok = false;

    if(upValue == 0.0 && lowValue == 0.0)
    {
        double range[2];
        model.Result()->GetScalarRange(range);

        lowValue = QInputDialog::getDouble(this, tr("Lower level of band to keep"),
                                         tr("Lower Value (inclusive):"), range[0], -2147483647, 2147483647, 12, &ok);
        upValue = QInputDialog::getDouble(this, tr("Upper level of band to keep"),
                                         tr("Upper Value (inclusive):"), range[1], -2147483647, 2147483647, 12, &ok);
        insideVal = QInputDialog::getDouble(this, tr("Value to give scalars inside the range"),
                                         tr("Inside Value:"), 1.0, -2147483647, 2147483647, 12, &ok);
        outsideVal = QInputDialog::getDouble(this, tr("Value to give scalars outside the range"),
                                         tr("Outside Value:"), 0.0, -2147483647, 2147483647, 12, &ok);
    }
    else
        ok = true;

    if(ok)
    {
        emit working(-1);
        printInfo("Binary Thresholding Scalars");
        model.Threshold(lowValue, upValue, insideVal, outsideVal);
        reset();
        emit done(-1);
    }
}

void milxQtModel::maskScalars(QString filename)
{
    if(!loaded)
        return;

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getOpenFileName(this,
                                               tr("Select File to Open"),
                                               path,
                                               tr(openModelExts.c_str()) );
    }

    if(filename.isEmpty())
        return;

    QPointer<milxQtModel> otherModel = new milxQtModel;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openModel(filename, otherModel);

    printInfo("Matching model information: ");
    otherModel->modelInfo();

    if(!success)
        return;

    emit working(-1);
    printInfo("Masking Scalars");
    model.Mask(otherModel->GetOutput());
    reset();
    emit done(-1);

    refresh();
}

void milxQtModel::gradient()
{
    if(!loaded || !GetScalars())
        return;

    emit working(-1);
    printInfo("Computing Gradient of Scalar Field");
    model.Gradient();
//    model.Result()->GetPointData()->SetActiveVectors("Gradient");
//    generateVectorField();
    refresh();
    emit done(-1);
}

void milxQtModel::smooth(int iterations)
{
    if(!loaded)
        return;

    bool ok = false;

    if(iterations == 0)
    {
        iterations = QInputDialog::getInt(this, tr("Please Provide the number of Iterations"),
                                          tr("Iterations:"), 50, 0, 10000, 1, &ok);
    }
    else
        ok = true;

    if(ok && iterations > 0)
    {
        emit working(-1);
        model.LaplacianSmoothing(iterations);

        printInfo("Smoothed for " + QString::number(iterations) + " iterations.");
        refresh();
        emit done(-1);
    }
}

void milxQtModel::smoothSinc(int iterations)
{
    if(!loaded)
        return;

    bool ok = false;

    if(iterations == 0)
    {
        iterations = QInputDialog::getInt(this, tr("Please Provide the number of Iterations"),
                                          tr("Iterations:"), 15, 0, 10000, 1, &ok);
    }
    else
        ok = true;

    if(ok && iterations > 0)
    {
        emit working(-1);
        model.WindowedSincSmoothing(iterations);
        printInfo("Smoothed for " + QString::number(iterations) + " iterations with a pass band of 0.1");
        refresh();
        emit done(-1);
    }
}

void milxQtModel::curvature()
{
    if(!loaded)
        return;

    emit working(-1);
    model.Curvature();
    scalarsSet = true;

    printInfo("Mean Curvature marked as scalars on mesh");

    refresh();
    emit done(-1);
}

void milxQtModel::transform(QString filename, bool inverse)
{
    if(!loaded)
        return;

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getOpenFileName(this,
                                               tr("Select File to Open"),
                                               path,
                                               tr("ITK Transform Files (*.txt *.tfm);;VTK Transform Files (*.trsf)") );

        QMessageBox msgBox;
        msgBox.setText("Invert Transform");
        msgBox.setInformativeText("Do you want to invert the transform?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        int ret = msgBox.exec();

        if(ret == QMessageBox::Yes)
            inverse = true;
    }

    if(filename.isEmpty())
        return;

    vtkSmartPointer<vtkTransform> transformer = vtkSmartPointer<vtkTransform>::New();
    vtkMatrix4x4* matrix;

    if(filename.contains(".trsf", Qt::CaseInsensitive))
        matrix = milx::File::OpenVTKTransform(filename.toStdString()); //!< Use VTK transform function
    else
        matrix = milx::File::OpenITKTransform(filename.toStdString()); //!< Use ITK to VTK transform function

    transformer->SetMatrix(matrix);

    if(inverse)
        transformer->Inverse();

    SetTransform(transformer);

    refresh();
}

void milxQtModel::removeScalars()
{
    if(!loaded)
        return;

    printInfo("Removing all arrays");
    model.RemoveScalars();

    refresh();
}

void milxQtModel::loadScalars(QString filename)
{
    if(!loaded)
        return;

    bool fileProvided = false;
    int ret = QMessageBox::No;
    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getOpenFileName(this,
                                               tr("Select File to Open"),
                                               path,
                                               tr(openModelExts.c_str()) );
        QMessageBox msgBox;
        msgBox.setText("Append Scalars into Model");
        msgBox.setInformativeText("Do you want to add scalars into model?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        ret = msgBox.exec();
    }
    else
        fileProvided = true;

    if(filename.isEmpty())
        return;

    QPointer<milxQtModel> otherModel = new milxQtModel;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openModel(filename, otherModel);

    if(!success)
    {
        printError("Failed reading model. Stopping.");
        return;
    }

    if(GetNumberOfPoints() != otherModel->GetNumberOfPoints())
    {
        printError("Models must have the same number of points! Stopping.");
        return;
    }

    if(!otherModel->GetScalars())
    {
        printError("Model " + filename + " does not have scalars");
        return;
    }

    //Check if multiple arrays present, if so ask if name wasnt provided initally
    int index = 0;
    if(otherModel->GetNumberOfArrays() > 1 && !fileProvided)
    {
        QDialog *arrayPossibilities = new QDialog(this);
        QComboBox *integerValues = new QComboBox(this);
        QPushButton *okButton = new QPushButton(this);
            okButton->setText("Ok");
        QLabel *lblMessage = new QLabel(this);
            lblMessage->setText("Choose Scalars from possibilities.");
        QFormLayout *formLayout = new QFormLayout(this);
            formLayout->addRow(lblMessage);
            formLayout->addRow(tr("&Select Scalars: "), integerValues);
            formLayout->addRow(okButton);
            arrayPossibilities->setLayout(formLayout);

        connect(okButton, SIGNAL(clicked(bool)), arrayPossibilities, SLOT(accept()));

        //Fill combox box
        for(int j = 0; j < otherModel->GetOutput()->GetPointData()->GetNumberOfArrays(); j ++)
            integerValues->addItem(otherModel->GetOutput()->GetPointData()->GetArray(j)->GetName());

        arrayPossibilities->exec();

        index = integerValues->currentIndex();
    }
    else if(otherModel->GetNumberOfArrays() > 1 && fileProvided) //copy all arrays
    {
        RemoveScalars();

        for(int j = 0; j < otherModel->GetOutput()->GetPointData()->GetNumberOfArrays(); j ++)
            AddArray(otherModel->GetOutput()->GetPointData()->GetArray(j));
    }

    if(ret == QMessageBox::Yes)
    {
        AddArray(otherModel->GetOutput()->GetPointData()->GetArray(index));
        SetActiveScalars(otherModel->GetOutput()->GetPointData()->GetArray(index)->GetName());
    }
    else
        SetScalars(otherModel->GetScalars());
    printInfo("Loaded Scalars from " + filename);

    if(modelled)
        updateLookupTable();
}

void milxQtModel::modelInfo()
{
    if(!loaded)
        return;

    cout << "Centroid: " << centroid() << " with Size: " << centroidSize() << endl;
    cout << "Covariance Matrix: " << endl << covarianceMatrix() << endl;
    printInfo("There are " + QString::number(model.Result()->GetNumberOfPoints()) + " points.");
    printInfo("There are " + QString::number(model.Result()->GetNumberOfPolys()) + " polygons.");
    printInfo("There are " + QString::number(model.Result()->GetNumberOfLines()) + " lines.");
    printInfo("There are " + QString::number(model.Result()->GetNumberOfStrips()) + " strips.");
    printInfo("There are " + QString::number(model.Result()->GetPointData()->GetNumberOfArrays()) + " arrays:");
    for(int j = 0; j < model.Result()->GetPointData()->GetNumberOfArrays(); j ++)
        printInfo(QString(model.Result()->GetPointData()->GetArray(j)->GetName()) + " ("
                  + QString::number(model.Result()->GetPointData()->GetArray(j)->GetNumberOfComponents())
                  + " component(s))");
}

void milxQtModel::matchInfo(QString filename, bool rescale, double factor)
{
    if(!loaded)
        return;

    printInfo("Initial model information: ");
    modelInfo();

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getOpenFileName(this,
                                               tr("Select File to Open"),
                                               path,
                                               tr(openModelExts.c_str()) );

        QMessageBox msgBox;
        msgBox.setText("Rescale");
        msgBox.setInformativeText("Do you want to match the scale also?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        int ret = msgBox.exec();

        if(ret == QMessageBox::No)
            rescale = false;
    }

    if(filename.isEmpty())
        return;

    QPointer<milxQtModel> otherModel = new milxQtModel;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openModel(filename, otherModel);

    printInfo("Matching model information: ");
    otherModel->modelInfo();

    if(!success)
        return;

    coordinate otherCentroid = otherModel->centroid();
    coordinate newCentroid = otherCentroid - centroid();
    double currentCentroidSize = centroidSize(true);
    double otherCentroidSize = otherModel->centroidSize(true);
    double scaling = otherCentroidSize/currentCentroidSize;
    if(rescale)
    {
        printInfo("Scale of Current model: " + QString::number(currentCentroidSize));
        printInfo("Scale of Other model: " + QString::number(otherCentroidSize));
        printInfo("Scale to be applied: " + QString::number(scaling));
        printInfo("Scale (with factor) to be applied: " + QString::number(factor*scaling));
    }

    vtkSmartPointer<vtkTransform> transformer = vtkSmartPointer<vtkTransform>::New();
    transformer->Translate(newCentroid[0], newCentroid[1], newCentroid[2]);
    if(rescale)
        transformer->Scale(factor*scaling, factor*scaling, factor*scaling);

    SetTransform(transformer);

    refresh();

    printInfo("Final model information: ");
    modelInfo();
}

void milxQtModel::registerICP(QString filename, bool similarity)
{
    if(!loaded)
        return;

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getOpenFileName(this,
                                               tr("Select File to Open"),
                                               path,
                                               tr(openModelExts.c_str()) );

        QMessageBox msgBox;
        msgBox.setText("Similarity Transform");
        msgBox.setInformativeText("Do you want to use similarity matching? Rigid Body will be used otherwise.");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        int ret = msgBox.exec();

        if(ret == QMessageBox::Yes)
            similarity = true;
    }

    if(filename.isEmpty())
        return;

    QPointer<milxQtModel> otherModel = new milxQtModel;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openModel(filename, otherModel);

    if(!success)
        return;

    emit working(-1);
    printInfo("MSE of points before ICP: " + QString::number( milx::Math<double>::MeanSquaredError(model.Result()->GetPoints(), otherModel->GetOutput()->GetPoints()) ));
    printInfo("Registering Meshes using ICP.");
    model.IterativeClosestPointsRegistration(otherModel->GetOutput(), similarity);
    refresh();
    emit done(-1);
    printInfo("Registration Done.");
    printInfo("MSE of points after ICP: " + QString::number( milx::Math<double>::MeanSquaredError(model.Result()->GetPoints(), otherModel->GetOutput()->GetPoints()) ));
}

void milxQtModel::registerLandmarks(QString filename, bool similarity)
{
    if(!loaded)
        return;

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getOpenFileName(this,
                                               tr("Select File to Open"),
                                               path,
                                               tr(openModelExts.c_str()) );

        QMessageBox msgBox;
        msgBox.setText("Similarity Transform");
        msgBox.setInformativeText("Do you want to use similarity matching? Rigid Body will be used otherwise.");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        int ret = msgBox.exec();

        if(ret == QMessageBox::Yes)
            similarity = true;
    }

    if(filename.isEmpty())
        return;

    QPointer<milxQtModel> otherModel = new milxQtModel;
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = reader->openModel(filename, otherModel);

    if(!success)
        return;

    if(otherModel->GetNumberOfPoints() != model.Result()->GetNumberOfPoints())
    {
        printError("Meshes must have same number of points to use this registration method.");
        return;
    }

    emit working(-1);
    printInfo("MSE of points before Landmark Transform: " + QString::number( milx::Math<double>::MeanSquaredError(model.Result()->GetPoints(), otherModel->GetOutput()->GetPoints()) ));
    printInfo("Registering Meshes using Landmark Transform.");
    model.LandmarkBasedRegistration(otherModel->GetOutput(), similarity);
    refresh();
    emit done(-1);
    printInfo("Registration Done.");
    printInfo("MSE of points after Landmark Transform: " + QString::number( milx::Math<double>::MeanSquaredError(model.Result()->GetPoints(), otherModel->GetOutput()->GetPoints()) ));
}

void milxQtModel::texture()
{
    if(!loaded)
        return;

    printInfo("Extracting Texture of Model as an Image.");
    if(modelActor->GetTexture())
    {
        vtkSmartPointer<vtkImageData> tex = modelActor->GetTexture()->GetInput();

        emit imageAvailable(tex, "Texture Map");
    }
    else
        printError("No Texture Found.");
}

void milxQtModel::voxelise()
{
    refresh();

    emit surfaceToImage(model.Result());
}

void milxQtModel::contour()
{
    if(!milxQtRenderWindow::contourWidget)
    {
        model.GenerateNormals(2); //for polygonal line placer (needs point norms) and interp (needs cell norms)

        milxQtRenderWindow::contourWidget = vtkSmartPointer<vtkContourWidget>::New();
        milxQtRenderWindow::contourWidget->SetInteractor(QVTKWidget::GetInteractor());

        printInfo("Contour Surface Mode enabled.\nLeft Click to place points, Right click to place end point.");
        printInfo("Delete key to delete point, Shift+Delete key to reset.");

        vtkSmartPointer<vtkOrientedGlyphContourRepresentation> rep = vtkOrientedGlyphContourRepresentation::SafeDownCast( milxQtRenderWindow::contourWidget->GetRepresentation() );
            rep->GetLinesProperty()->SetColor(1, 0.2, 0);
            rep->GetLinesProperty()->SetLineWidth(3.0);

        vtkSmartPointer<vtkPolygonalSurfacePointPlacer> pointPlacer = vtkSmartPointer<vtkPolygonalSurfacePointPlacer>::New();
            pointPlacer->AddProp(modelActor);
            pointPlacer->GetPolys()->AddItem(model.Result());
//        #if VTK_MAJOR_VERSION > 5
//            pointPlacer->SnapToClosestPointOn();
//        #endif
        rep->SetPointPlacer(pointPlacer);

        vtkSmartPointer<vtkPolygonalSurfaceContourLineInterpolator> interpolator = vtkSmartPointer<vtkPolygonalSurfaceContourLineInterpolator>::New();
            interpolator->GetPolys()->AddItem(model.Result());
        rep->SetLineInterpolator(interpolator);

        milxQtRenderWindow::contourWidget->SetRepresentation(rep);
    }

    if(milxQtRenderWindow::contourAct->isChecked())
    {
        milxQtRenderWindow::contourWidget->EnabledOn();
        printInfo("Enabled Contouring");
    }
    else
    {
        milxQtRenderWindow::contourWidget->EnabledOff();
        printInfo("Disabled Contouring");
    }

    refresh();
}

void milxQtModel::selectPointsInContour()
{
    if(!milxQtRenderWindow::contourAct->isChecked())
    {
        printInfo("Contouring not enabled. Stopping.");
        return;
    }

    vtkSmartPointer<vtkContourRepresentation> rep = milxQtRenderWindow::contourWidget->GetContourRepresentation();
    vtkSmartPointer<vtkPolyData> contourMesh = rep->GetContourRepresentationAsPolyData();

    model.DistanceMap(contourMesh->GetPoints());

    refresh();
}

void milxQtModel::weightGaussianFromContour(float stddev, float clampValue)
{
  if(!milxQtRenderWindow::contourAct->isChecked())
    {
      printInfo("Contouring not enabled. Stopping.");
      return;
    }

  //get stddev for Gaussian
  bool ok1 = false, ok2 = false;
  if(stddev < 0)
      stddev = QInputDialog::getDouble(this, tr("Please enter the standard deviation of the Gaussian"), tr("Standard Dev.:"), 3.0, 0.0, 1000.0, 3, &ok1);
  if(clampValue < 0)
    clampValue = QInputDialog::getDouble(this, tr("Please enter the clamp (lowest allowed) value"), tr("Clamp Value:"), 1e-6, 0, 1.0, 10, &ok2);

  if(!ok1 || !ok2)
      return;

  //get dmap first
  selectPointsInContour();

  double range[2];
  model.GetScalarRange(range);
  thresholdScalars(0.0, range[1], 0.0);

  model.GaussianWeightScalars(stddev, clampValue);

  updateLookupTable();

  refresh();
}

void milxQtModel::changeColour(float red, float green, float blue)
{
    if(red < 0 || green < 0 || blue < 0) //Already provided so change
    {
        QPointer<QColorDialog> newColours = new QColorDialog(this);

        QColor newColour = newColours->getColor(QColor(defaultChannelValue, defaultChannelValue, defaultChannelValue),
                                                this, "Choose New Colour", QColorDialog::ShowAlphaChannel);

        if(newColour.isValid())
        {
            printInfo("Colour selected was: " + QString::number(newColour.redF()) + ", " + QString::number(newColour.greenF()) + ", " + QString::number(newColour.blueF()));

            colourRed = newColour.redF();
            colourGreen = newColour.greenF();
            colourBlue = newColour.blueF();

            SetOpacity(newColour.alphaF());

            reset();
        }
    }
    else
    {
        colourRed = red;
        colourGreen = green;
        colourBlue = blue;

        reset();

        return;
    }
}

void milxQtModel::changeOpacity(float opacity)
{
    bool ok = false;

    if(opacity < 0)
        opacity = QInputDialog::getDouble(this, tr("Please enter the new opacity"),
                                          tr("Opacity:"), opacity, 0.0, 1.0, 3, &ok);

    if(ok)
    {
        SetOpacity(opacity);

        reset();
    }
}

void milxQtModel::toggleInterpolation()
{
    model.GenerateNormals();

    if(interpAct->isChecked())
        modelActor->GetProperty()->SetInterpolationToPhong();
    else
        modelActor->GetProperty()->SetInterpolationToFlat();

    QString shadingStr = modelActor->GetProperty()->GetInterpolationAsString();
    printInfo("Using " + shadingStr + " Shading");

    refresh();
}

void milxQtModel::normals(const bool turnOn)
{
    if(turnOn)
        normalsAct->setChecked(true);

    if(normalsAct->isChecked())
    {
        // Double point normals
        vtkSmartPointer<vtkFloatArray> normalData2 = vtkFloatArray::SafeDownCast(model.Result()->GetPointData()->GetNormals());

        if(!normalData2)
        {
            printInfo("No normals found. Generating point normals.");
            generateNormals();
            normalData2 = vtkFloatArray::SafeDownCast(model.Result()->GetPointData()->GetNormals());
        }

        printInfo("There are " + QString::number(normalData2->GetNumberOfComponents()) + " components in the current model (double Point Normals)");

        if(!normalsModel)
            normalsModel = new milxQtModel;
        normalsModel->SetPoints(model.Result()->GetPoints());
        normalsModel->SetVectors(normalData2);
        normalsModel->generateVectorField();

        milxQtRenderWindow::AddActor(normalsModel->GetActor());
        milxQtRenderWindow::generateRender();
    }
    else
    {
        if(normalsModel)
            milxQtRenderWindow::RemoveActor(normalsModel->GetActor());
    }
}

void milxQtModel::rotate(bool xAxis, bool yAxis, bool zAxis, float angle)
{
    if(!loaded)
        return;

    if(!xAxis && !yAxis && !zAxis)
    {
        QMessageBox msgBoxX, msgBoxY, msgBoxZ;
        msgBoxX.setText("The model will be rotated");
        msgBoxX.setInformativeText("Do you want to flip the x-axis?");
        msgBoxX.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxX.setDefaultButton(QMessageBox::No);
        msgBoxY.setText("The model will be rotated");
        msgBoxY.setInformativeText("Do you want to flip the y-axis?");
        msgBoxY.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxY.setDefaultButton(QMessageBox::Yes);
        msgBoxZ.setText("The model will be rotated");
        msgBoxZ.setInformativeText("Do you want to flip the z-axis?");
        msgBoxZ.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxZ.setDefaultButton(QMessageBox::No);
        int retX = msgBoxX.exec();
        int retY = msgBoxY.exec();
        int retZ = msgBoxZ.exec();

        bool ok1 = false;
        angle = QInputDialog::getDouble(this, tr("Enter angle to rotate by"),
                                         tr("Angle:"), 90, -2147483647, 2147483647, 5, &ok1);

        if(retX == QMessageBox::Yes)
            xAxis = true;
        if(retY == QMessageBox::Yes)
            yAxis = true;
        if(retZ == QMessageBox::Yes)
            zAxis = true;
        if(!ok1)
            return;
    }

    printInfo("Rotating Model along chosen axes from the model centroid.");
    emit working(-1);
    model.Rotate(xAxis, yAxis, zAxis, angle, centroid());
    refresh();
    emit done(-1);
}

void milxQtModel::flip(bool xAxis, bool yAxis, bool zAxis)
{
    if(!loaded)
        return;

    if(!xAxis && !yAxis && !zAxis)
    {
        QMessageBox msgBoxX, msgBoxY, msgBoxZ;
        msgBoxX.setText("The model will be flipped");
        msgBoxX.setInformativeText("Do you want to flip the x-axis?");
        msgBoxX.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxX.setDefaultButton(QMessageBox::No);
        msgBoxY.setText("The model will be flipped");
        msgBoxY.setInformativeText("Do you want to flip the y-axis?");
        msgBoxY.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxY.setDefaultButton(QMessageBox::Yes);
        msgBoxZ.setText("The model will be flipped");
        msgBoxZ.setInformativeText("Do you want to flip the z-axis?");
        msgBoxZ.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBoxZ.setDefaultButton(QMessageBox::No);
        int retX = msgBoxX.exec();
        int retY = msgBoxY.exec();
        int retZ = msgBoxZ.exec();

        if(retX == QMessageBox::Yes)
            xAxis = true;
        if(retY == QMessageBox::Yes)
            yAxis = true;
        if(retZ == QMessageBox::Yes)
            zAxis = true;
    }

    printInfo("Flipping/Reflecting Model along chosen axes.");
    emit working(-1);
    model.Flip(xAxis, yAxis, zAxis);
    refresh();
    emit done(-1);
}

void milxQtModel::background(bool white)
{
    if(milxQtRenderWindow::backgroundAct->isChecked() || white)
    {
        if(scaleBefore)
        {
            scale->GetLabelTextProperty()->SetColor(0, 0, 0);
            scale->GetTitleTextProperty()->SetColor(0, 0, 0);
        }
        if(outlineBefore)
            outlineActor->GetProperty()->SetColor(0, 0, 0);
        if(cubeAxesBefore)
            cubeAxesActor->GetProperty()->SetColor(0, 0, 0);
    }
    else
    {
        if(scaleBefore)
        {
            scale->GetLabelTextProperty()->SetColor(1.0, 1.0, 1.0);
            scale->GetTitleTextProperty()->SetColor(1.0, 1.0, 1.0);
        }
        if(outlineBefore)
            outlineActor->GetProperty()->SetColor(1, 1, 1);
        if(cubeAxesBefore)
            cubeAxesActor->GetProperty()->SetColor(1, 1, 1);
    }

    milxQtRenderWindow::background(white);
}

void milxQtModel::enableOutline(vtkDataObject *dataOutline)
{
    outlineMesh = vtkSmartPointer<vtkOutlineFilter>::New();
    if(dataOutline)
    #if VTK_MAJOR_VERSION <=5
        outlineMesh->SetInput(dataOutline);
    #else
        outlineMesh->SetInputData(dataOutline);
    #endif // VTK_MAJOR_VERSION
    else
    #if VTK_MAJOR_VERSION <=5
        outlineMesh->SetInput(model.Result());
    #else
        outlineMesh->SetInputData(model.Result());
    #endif
    outlineMesh->Update();

    if(!outlineModel)
        outlineModel = new milxQtModel;
    outlineModel->SetInput(outlineMesh->GetOutput());
    outlineModel->generateModel();

    outlineActor = outlineModel->GetActor();
    if(milxQtRenderWindow::backgroundAct->isChecked())
        outlineActor->GetProperty()->SetColor(0, 0, 0);
    else
        outlineActor->GetProperty()->SetColor(1, 1, 1);

    milxQtRenderWindow::AddActor(outlineActor);
    milxQtRenderWindow::generateRender();

    outlineBefore = true;
    outlineAct->setChecked(true);
}

void milxQtModel::disableOutline()
{
    if(!outlineActor)
        return; //nothing to do then

    RemoveActor(outlineActor);
    outlineBefore = false;
}

void milxQtModel::outlineDisplay()
{
    if(outlineAct->isChecked())
        enableOutline();
    else
        disableOutline();

    Render();
}

void milxQtModel::enableCubeAxes(double *range, double *bounds)
{
    cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
    if(bounds)
        cubeAxesActor->SetBounds(bounds);
    else
        cubeAxesActor->SetBounds(model.Result()->GetBounds());
    cubeAxesActor->SetCamera(milxQtRenderWindow::GetRenderer()->GetActiveCamera());
    if(range)
    {
        #if (VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 6))
            cubeAxesActor->SetXAxisRange(range[0], range[1]);
            cubeAxesActor->SetYAxisRange(range[2], range[3]);
            cubeAxesActor->SetZAxisRange(range[4], range[5]);
        #else
            printWarning("A scaling has been applied to the z-axis. Upgrade your VTK version to at least 5.8.x or above.");
        #endif
    }
    cubeAxesActor->SetFlyModeToOuterEdges();
    if(milxQtRenderWindow::backgroundAct->isChecked())
        cubeAxesActor->GetProperty()->SetColor(0, 0, 0);
    else
        cubeAxesActor->GetProperty()->SetColor(1, 1, 1);

    milxQtRenderWindow::AddActor(cubeAxesActor);
    milxQtRenderWindow::generateRender();

    cubeAxesBefore = true;
    cubeAxesAct->setChecked(true);
}

void milxQtModel::disableCubeAxes()
{
    if(!cubeAxesActor)
        return; //nothing to do then

    RemoveActor(cubeAxesActor);
    cubeAxesBefore = false;
}

void milxQtModel::cubeAxesDisplay(double *range)
{
    if(cubeAxesAct->isChecked())
        enableCubeAxes(range);
    else
        disableCubeAxes();

    Render();
}

void milxQtModel::removeOverlays()
{
    vtkSmartPointer<vtkActorCollection> collection = renderer->GetActors();

    const size_t n = collection->GetNumberOfItems();
    printDebug("Detected "+QString::number(n)+" overlays");

    collection->InitTraversal();
    for(size_t j = 0; j < n; j ++)
    {
        vtkSmartPointer<vtkActor> actor = collection->GetNextItem();
        if(modelActor != actor)
            removeModelActor(actor);
    }

    for(int j = 0; j < imageActors.size(); j ++)
        RemoveActor(imageActors[j].imageActor);
    imageActors.clear();

    Render();
}

void milxQtModel::enableScale(QString title, const bool quiet, double minRange, double maxRange, int noOfLabels)
{
    if(!milxQtRenderWindow::scale)
        milxQtRenderWindow::scale = vtkSmartPointer<vtkScalarBarActor>::New();
    if(!milxQtRenderWindow::scalarBar)
        milxQtRenderWindow::scalarBar = vtkSmartPointer<vtkScalarBarWidget>::New();

    ///Ask if use default scalar LUT
    QMessageBox msgBox;
    msgBox.setText("An auto adjusted bar is about to be created");
    msgBox.setInformativeText("Would you like to customise the bar?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);

    int ret = QMessageBox::No;
    if(!quiet)
      ret = msgBox.exec();

    const float barWidth = 0.1, barHeight = 0.7;

    vtkSmartPointer<vtkLogLookupTable> logLookupTable;
    if(milxQtRenderWindow::logScale)
    {
        printInfo("Detected log scale.");
        logLookupTable = vtkSmartPointer<vtkLogLookupTable>::New();
        logLookupTable->DeepCopy(vtkLookupTable::SafeDownCast(lookupTable));
    }

    if(ret == QMessageBox::Yes && !quiet)
    {
        bool ok1 = false, ok2 = false, ok3 = false, ok4 = false;

        minRange = QInputDialog::getDouble(this, tr("Enter Table Range of new Lookup Table"),
                                         tr("Minimum:"), 0, -2147483647, 2147483647, 5, &ok1);
        maxRange = QInputDialog::getDouble(this, tr("Enter Table Range of new Lookup Table"),
                                         tr("Maximum:"), 1, -2147483647, 2147483647, 5, &ok2);
        noOfLabels = QInputDialog::getInt(this, tr("How many labels to show"),
                                          tr("Labels:"), noOfLabels, 0, 99, 1, &ok3);
        title = QInputDialog::getText(this, tr("Title of Bar"),
                                          tr("Title:"), QLineEdit::Normal,
                                          title, &ok4);

        if(!ok1 || !ok2 || !ok3 || !ok4)
            return;

        if(milxQtRenderWindow::logScale)
            scale->SetLookupTable(logLookupTable);
        else
            scale->SetLookupTable(lookupTable);
        scale->SetNumberOfLabels(noOfLabels);
        modelMapper->SetLookupTable(lookupTable);
        modelMapper->SetScalarRange( minRange, maxRange );
        customScalarBar = true;
    }
    else if(quiet && minRange != maxRange)
    {
        printInfo("Using custom scalar range for scalars.");
        if(milxQtRenderWindow::logScale)
            scale->SetLookupTable(logLookupTable);
        else
            scale->SetLookupTable(lookupTable);
        scale->SetNumberOfLabels(noOfLabels);
        modelMapper->SetLookupTable(lookupTable);
        modelMapper->SetScalarRange( minRange, maxRange );
        customScalarBar = true;
    }
    else
    {
        modelMapper->SetScalarRange( model.Result()->GetScalarRange() );
        if(milxQtRenderWindow::logScale)
            scale->SetLookupTable(logLookupTable);
        else
            scale->SetLookupTable(lookupTable);
        scale->SetNumberOfLabels(3);
        customScalarBar = false;
    }

    scale->SetTitle(title.toStdString().c_str());
    scale->GetLabelTextProperty()->SetFontFamilyToArial();
    scale->GetLabelTextProperty()->SetFontSize(8);
    scale->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scale->GetPositionCoordinate()->SetValue(.2,.05);
    scale->SetWidth( barWidth );
    scale->SetHeight( barHeight );
    scale->SetPosition( 0.99 - barWidth, 0.1 );
    scale->SetLabelFormat("%-#6.3f");
    scale->GetTitleTextProperty()->SetFontFamilyToArial();
    scale->GetTitleTextProperty()->SetFontSize(8);
    scale->GetLabelTextProperty()->SetJustificationToCentered();

    if(milxQtRenderWindow::backgroundAct->isChecked())
    {
        scale->GetLabelTextProperty()->SetColor(0, 0, 0);
        scale->GetTitleTextProperty()->SetColor(0, 0, 0);
    }

    //Add scale to scale widget
    scalarBar->SetInteractor(QVTKWidget::GetInteractor());
    scalarBar->SetScalarBarActor(scale);
    scalarBar->EnabledOn();

    scaleBefore = true;
    scaleAct->setChecked(true);
}

void milxQtModel::scaleDisplay(const bool forceDisplay)
{
    if(model.Result()->GetPointData()->GetScalars() == NULL)
    {
        printError("No scalars present.");
        scaleAct->setChecked(false); //reset
        return;
    }

    if(scaleAct->isChecked() || forceDisplay)
        enableScale(model.Result()->GetPointData()->GetScalars()->GetName(), forceDisplay);
    else
        disableScale();

    Render();
}

void milxQtModel::outputScalars(QString filename)
{
    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileOpener = new QFileDialog(this);
        filename = fileOpener->getSaveFileName(this,
                                               tr("Select File to Save"),
                                               path,
                                               tr("CSV Files (*.txt *.csv)") );
    }

    if(filename.isEmpty())
        return;

    QPointer<milxQtFile> scalarsWriter = new milxQtFile;

    scalarsWriter->saveScalarsOfModel(filename, this);
}

void milxQtModel::viewToXYPlane()
{
    if(modelled)
    {
        if(!computedCentroid)
            centroid();

        //viewToAxial
        vtkCamera *cam = milxQtRenderWindow::GetRenderer()->GetActiveCamera();
        if(!milxQtRenderWindow::GetRenderer()->IsActiveCameraCreated())
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            cam->SetFocalPoint(modelCentroid.data_block());
            if(orientationView == RADIOLOGICAL)
            {
                cam->SetPosition(0,0,centroidSize());
                cam->SetViewUp(0,-centroidSize(),0);
            }
            else //Neurological
            {
                cam->SetPosition(0,0,-centroidSize());
                cam->SetViewUp(0,centroidSize(),0);
            }
            milxQtRenderWindow::GetRenderer()->ResetCamera();
            milxQtRenderWindow::Render();

            currentView = AXIAL;
    }
}

void milxQtModel::viewToZXPlane()
{
    if(modelled)
    {
        if(!computedCentroid)
            centroid();

        //viewToCoronal
        vtkCamera *cam = milxQtRenderWindow::GetRenderer()->GetActiveCamera();
        if(!milxQtRenderWindow::GetRenderer()->IsActiveCameraCreated())
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            cam->SetFocalPoint(modelCentroid.data_block());
            if(orientationView == RADIOLOGICAL)
            {
                cam->SetPosition(0,-centroidSize(),0);
                cam->SetViewUp(0,0,centroidSize());
            }
            else //Neurological
            {
                cam->SetPosition(0,centroidSize(),0);
                cam->SetViewUp(0,0,centroidSize());
            }
            milxQtRenderWindow::GetRenderer()->ResetCamera();
            milxQtRenderWindow::Render();

            currentView = CORONAL;
    }
}

void milxQtModel::viewToZYPlane()
{
    if(modelled)
    {
        if(!computedCentroid)
            centroid();

        //viewToSagittal
        vtkCamera *cam = milxQtRenderWindow::GetRenderer()->GetActiveCamera();
        if(!milxQtRenderWindow::GetRenderer()->IsActiveCameraCreated())
            milxQtRenderWindow::GetRenderer()->ResetCamera(); //just in case if new camera
            cam->SetFocalPoint(modelCentroid.data_block());
            cam->SetPosition(-centroidSize(),0,0);
            cam->SetViewUp(0,0,centroidSize());
            milxQtRenderWindow::GetRenderer()->ResetCamera();
            milxQtRenderWindow::Render();

            currentView = SAGITTAL;
    }
}

void milxQtModel::showArray(QAction * action)
{
    showArray(action->text());
}

void milxQtModel::showArray(const QString arrayName)
{
    ///Check if array has scalars or vectors
    if(model.Result()->GetPointData()->GetArray(arrayName.toStdString().c_str())->GetNumberOfComponents() == 1)
    {
        printInfo("Loaded Array " + arrayName + " as scalars");
        model.Result()->GetPointData()->SetActiveAttribute( arrayName.toStdString().c_str(), vtkDataSetAttributes::SCALARS );
    }
    else if(model.Result()->GetPointData()->GetArray(arrayName.toStdString().c_str())->GetNumberOfComponents() == 3)
    {
        printInfo("Loaded Array " + arrayName + " as vectors");
        model.Result()->GetPointData()->SetActiveAttribute( arrayName.toStdString().c_str(), vtkDataSetAttributes::VECTORS );
        printInfo("Use 'Generate Vector Field' to see the vectors.");
    }
    else
        printInfo("Cannot Load Array " + arrayName + ". Too many components.");

    if(modelled)
        updateLookupTable();
}

void milxQtModel::updateLookupTable()
{
    double range[2] = {0, 1.0};

    if(GetScalars())
    {
        model.Result()->GetScalarRange(range);
        if(GetScalars()->GetNumberOfComponents() > 1)
            printInfo("Found multi-component scalars. Not using colour map and using colours directly.");
    }

    ///\todo Force Opacity to 1.0 here, better more clearer way?
    ///Alpha channel for colormap not allow for models
    double zeroRGBA[4];
    lookupTable->GetTableValue(0, zeroRGBA);
    if(zeroRGBA[3] == 0.0)
        lookupTable->SetTableValue(0, zeroRGBA[0], zeroRGBA[1], zeroRGBA[2], 1.0); //opaque

    modelMapper->SetLookupTable(lookupTable);
    modelMapper->SetScalarRange( range[0], range[1] );

    scaleDisplay();
}

void milxQtModel::createMenu(QMenu *menu)
{
    if(!menu)
        return;

    menu->clear();
    menu->addMenu(basicContextMenu());

    if(!extActionsToAdd.empty())
        menu->addSeparator()->setText(tr("Extensions"));
    foreach(QAction *currAct, extActionsToAdd)
    {
        menu->addAction(currAct);
    }
    menu->addSeparator();
    menu->addAction(scaleAct);
    menu->addMenu(milxQtRenderWindow::contourMenu);
    menu->addMenu(milxQtRenderWindow::windowPropertiesMenu);

    menu->addSeparator();
    foreach(QAction *currAct, milxQtWindow::actionsToAppend)
    {
      menu->addAction(currAct);
    }
    foreach(QMenu *currMenu, milxQtWindow::menusToAppend)
    {
      menu->addMenu(currMenu);
    }
    menu->addAction(milxQtRenderWindow::refreshAct);
    menu->addAction(milxQtRenderWindow::resetAct);

    milxQtRenderWindow::contourMenu->addAction(selectPointsAct);
    milxQtRenderWindow::contourMenu->addAction(weightAct);

    //disabling
    milxQtRenderWindow::contourPolyDataAct->setDisabled(!contourAct->isChecked());
    milxQtRenderWindow::contourNodePolyDataAct->setDisabled(!contourAct->isChecked());
    milxQtRenderWindow::contourInitAct->setDisabled(!contourAct->isChecked());
    selectPointsAct->setDisabled(!contourAct->isChecked());
    weightAct->setDisabled(!contourAct->isChecked());
    milxQtRenderWindow::contourInitAct->setDisabled(true); ///\todo contour init doesn't work. VTK bug?
}

void milxQtModel::updateCoords(vtkObject *obj)
{
    ///Get interactor
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    QString message = "";

    ///Get event position
    ///Code initial by Mark Wyszomierski 2003-2007 @ devsample
    ///Modified by Shekhar Chandra
    ///Do the pick. It will return a non-zero value if we intersected the bounding box.
    if (dataPicker->Pick(iren->GetEventPosition()[0],
                         iren->GetEventPosition()[1],
                         0,  // always zero.
                         renderer))
    {
        // Get the volume index within the entire volume now.
        vtkIdType nVolIdx = dataPicker->GetPointId();

        if(nVolIdx >= 0) //-1 means no point picked
        {
            double *coords = model.Result()->GetPoint(nVolIdx);

            if(model.Result()->GetPointData()->GetScalars())
            {
                double *scalar = model.Result()->GetPointData()->GetScalars()->GetTuple(nVolIdx);

                if(model.Result()->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
                    message = "Point " + QString::number(nVolIdx) + ": (" + QString::number(coords[0]) + ", " + QString::number(coords[1]) + ", " + QString::number(coords[2]) + ") = "
                              + "[" + QString::number(scalar[0]) + ", " + QString::number(scalar[1]) + ", " + QString::number(scalar[2]) + "]";
                else if(model.Result()->GetPointData()->GetScalars()->GetNumberOfComponents() == 2)
                    message = "Point " + QString::number(nVolIdx) + ": (" + QString::number(coords[0]) + ", " + QString::number(coords[1]) + ", " + QString::number(coords[2]) + ") = "
                              + "[" + QString::number(scalar[0]) + ", " + QString::number(scalar[1]) + "]";
                else
                    message = "Point " + QString::number(nVolIdx) + ": (" + QString::number(coords[0]) + ", " + QString::number(coords[1]) + ", " + QString::number(coords[2]) + ") = "
                              + QString::number(scalar[0]);
            }
            else
            {
                message = "Point " + QString::number(nVolIdx) + ": (" + QString::number(coords[0]) + ", " + QString::number(coords[1]) + ", " + QString::number(coords[2]) + ")";
            }
        }
    }

    ///Write message to status bar
    updateBar->showMessage(message);
}

void milxQtModel::copyToContextMenu(QMenu *copyMenu)
{
    milxQtWindow::copyToContextMenu(copyMenu);
    createCustomConnections(copyMenu->actions());
    printDebug("Model Copy to Context");
}

//Helper Members
void milxQtModel::resetFlags()
{
    computedCentroid = false;
    computedCovariance = false;
    if(!modelled)
    {
        appended = false;
        labelled = false;
        scalarsSet = false;
        customScalarBar = false;
        transformed = false;
        largeMode = false;
        scaleBefore = false;
        outlineBefore = false;
        cubeAxesBefore = false;
        undidProcess = false;
    }
}

void milxQtModel::createActions()
{
    //Generate
    generateMenu = new QMenu(this);
    generateMenu->setTitle(QApplication::translate("Model", "Generate", 0, QApplication::UnicodeUTF8));
    genModelAct = new QAction(this);
    genModelAct->setText(QApplication::translate("Model", "Model", 0, QApplication::UnicodeUTF8));
    genModelAct->setShortcut(tr("Shift+Alt+m"));
    genVerticesAct = new QAction(this);
    genVerticesAct->setText(QApplication::translate("Model", "Vertices (Vertex Glyphs)", 0, QApplication::UnicodeUTF8));
    genVerticesAct->setShortcut(tr("Shift+Alt+v"));
    genNormalsAct = new QAction(this);
    genNormalsAct->setText(QApplication::translate("Model", "Normals", 0, QApplication::UnicodeUTF8));
    genNormalsAct->setShortcut(tr("Shift+Alt+n"));
    genVectorsAct = new QAction(this);
    genVectorsAct->setText(QApplication::translate("Model", "Vector Field", 0, QApplication::UnicodeUTF8));
    genVectorsAct->setShortcut(tr("Shift+Alt+v"));
    genTensorsAct = new QAction(this);
    genTensorsAct->setText(QApplication::translate("Model", "Tensor Field", 0, QApplication::UnicodeUTF8));
    genTensorsAct->setShortcut(tr("Shift+Alt+t"));
    genHedgehogAct = new QAction(this);
    genHedgehogAct->setText(QApplication::translate("Model", "Hedgehog Field", 0, QApplication::UnicodeUTF8));
    genHedgehogAct->setShortcut(tr("Shift+Alt+h"));
    genDelaunayAct = new QAction(this);
    genDelaunayAct->setText(QApplication::translate("Model", "Delaunay Graph", 0, QApplication::UnicodeUTF8));
    genDelaunayAct->setShortcut(tr("Shift+Alt+g"));
    genDelaunayTri2DAct = new QAction(this);
    genDelaunayTri2DAct->setText(QApplication::translate("Model", "Delaunay 2D Triangulation", 0, QApplication::UnicodeUTF8));
    genDelaunayTri2DAct->setShortcut(tr("Shift+Ctrl+t"));
    genDelaunayTriAct = new QAction(this);
    genDelaunayTriAct->setText(QApplication::translate("Model", "Delaunay Triangulation", 0, QApplication::UnicodeUTF8));
    genDelaunayTriAct->setShortcut(tr("Shift+Alt+t"));
    genDelaunayTriAct->setCheckable(true);
    genDelaunayTriAct->setChecked(false);
    genPointModelAct = new QAction(this);
    genPointModelAct->setText(QApplication::translate("Model", "Point Model", 0, QApplication::UnicodeUTF8));
    genPointModelAct->setShortcut(tr("Shift+Alt+p"));
    genLabelsAct = new QAction(this);
    genLabelsAct->setText(QApplication::translate("Model", "Vertex Labels", 0, QApplication::UnicodeUTF8));
    genLabelsAct->setShortcut(tr("Shift+Alt+l"));
    genPointIDsAct = new QAction(this);
    genPointIDsAct->setText(QApplication::translate("Model", "Point ID Scalars", 0, QApplication::UnicodeUTF8));
    genPointIDsAct->setShortcut(tr("Shift+Alt+i"));
    genSamplesAct = new QAction(this);
    genSamplesAct->setText(QApplication::translate("Model", "Sampled Points", 0, QApplication::UnicodeUTF8));
    genSamplesAct->setShortcut(tr("Shift+Alt+s"));
    genKmeansAct = new QAction(this);
    genKmeansAct->setText(QApplication::translate("Model", "K-Means Point Clustering", 0, QApplication::UnicodeUTF8));
    genKmeansAct->setShortcut(tr("Shift+Alt+k"));
    genQuantiseAct = new QAction(this);
    genQuantiseAct->setText(QApplication::translate("Model", "Snap/Quantise Points to a Grid", 0, QApplication::UnicodeUTF8));
    genQuantiseAct->setShortcut(tr("Shift+Alt+q"));
    genCapBoundaries = new QAction(this);
    genCapBoundaries->setText(QApplication::translate("Model", "Cap the Boundaries of Model", 0, QApplication::UnicodeUTF8));
    genCapBoundaries->setShortcut(tr("Shift+Alt+c"));
    genRegionLabels = new QAction(this);
    genRegionLabels->setText(QApplication::translate("Model", "Label Unconnected Regions", 0, QApplication::UnicodeUTF8));
    genRegionLabels->setShortcut(tr("Shift+Alt+u"));
    genElevation = new QAction(this);
    genElevation->setText(QApplication::translate("Model", "Elevation", 0, QApplication::UnicodeUTF8));
    genElevation->setShortcut(tr("Shift+Alt+e"));
    genReebGraph = new QAction(this);
    genReebGraph->setText(QApplication::translate("Model", "Reeb Graph", 0, QApplication::UnicodeUTF8));
    genReebGraph->setShortcut(tr("Shift+Alt+r"));
    genReebGraph->setDisabled(true); //Not working so disable tmp

    //Operations
    operateMenu = new QMenu(this);
    operateMenu->setTitle(QApplication::translate("Model", "Operations", 0, QApplication::UnicodeUTF8));
    cleanAct = new QAction(this);
    cleanAct->setText(QApplication::translate("Model", "&Clean Mesh", 0, QApplication::UnicodeUTF8));
    cleanAct->setShortcut(tr("Alt+c"));
    triAct = new QAction(this);
    triAct->setText(QApplication::translate("Model", "&Triangulate", 0, QApplication::UnicodeUTF8));
    triAct->setShortcut(tr("Alt+t"));
    decimateAct = new QAction(this);
    decimateAct->setText(QApplication::translate("Model", "&Decimate Mesh", 0, QApplication::UnicodeUTF8));
    decimateAct->setShortcut(tr("Alt+d"));
    quadricDecimateAct = new QAction(this);
    quadricDecimateAct->setText(QApplication::translate("Model", "&Quadric Decimate Mesh", 0, QApplication::UnicodeUTF8));
    quadricDecimateAct->setShortcut(tr("Alt+q"));
    clusterDecimateAct = new QAction(this);
    clusterDecimateAct->setText(QApplication::translate("Model", "&Quadric Cluster Decimate Mesh", 0, QApplication::UnicodeUTF8));
    clusterDecimateAct->setShortcut(tr("Shift+Alt+q"));
    smoothAct = new QAction(this);
    smoothAct->setText(QApplication::translate("Model", "Smooth/&Fair (Laplacian) Mesh", 0, QApplication::UnicodeUTF8));
    smoothAct->setShortcut(tr("Alt+f"));
    smoothSincAct = new QAction(this);
    smoothSincAct->setText(QApplication::translate("Model", "Smooth/&Fair (Windowed Sinc) Mesh", 0, QApplication::UnicodeUTF8));
    smoothSincAct->setShortcut(tr("Shift+Alt+f"));
    curvatureAct = new QAction(this);
    curvatureAct->setText(QApplication::translate("Model", "Mean Curvature", 0, QApplication::UnicodeUTF8));
    curvatureAct->setShortcut(tr("Shift+Alt+c"));

    //Transform
    transformMenu = new QMenu(this);
    transformMenu->setTitle(QApplication::translate("Model", "Transforms", 0, QApplication::UnicodeUTF8));
    transformAct = new QAction(this);
    transformAct->setText(QApplication::translate("Model", "Transform (via File) from ...", 0, QApplication::UnicodeUTF8));
    transformAct->setShortcut(tr("Alt+i"));
    matchAct = new QAction(this);
    matchAct->setText(QApplication::translate("Model", "&Match Info to ...", 0, QApplication::UnicodeUTF8));
    matchAct->setShortcut(tr("Alt+a"));
    registerAct = new QAction(this);
    registerAct->setText(QApplication::translate("Model", "&Register (via ICP) to ...", 0, QApplication::UnicodeUTF8));
    registerAct->setShortcut(tr("Alt+r"));
    registerLandmarkAct = new QAction(this);
    registerLandmarkAct->setText(QApplication::translate("Model", "&Register (via Landmarks) to ...", 0, QApplication::UnicodeUTF8));
    registerLandmarkAct->setShortcut(tr("Shift+Alt+r"));
    voxeliseAct = new QAction(this);
    voxeliseAct->setText(QApplication::translate("Model", "&Voxelise Model", 0, QApplication::UnicodeUTF8));
    voxeliseAct->setShortcut(tr("Alt+v"));

    //Scalars
    scalarMenu = new QMenu(this);
    scalarMenu->setTitle(QApplication::translate("Model", "Scalar Operations", 0, QApplication::UnicodeUTF8));
    renameScalarsAct = new QAction(this);
    renameScalarsAct->setText(QApplication::translate("Model", "Rename Scalars", 0, QApplication::UnicodeUTF8));
    renameScalarsAct->setShortcut(tr("Shift+Alt+r"));
    thresholdAct = new QAction(this);
    thresholdAct->setText(QApplication::translate("Model", "Clip Mesh based on Scalar Threshold", 0, QApplication::UnicodeUTF8));
    thresholdAct->setShortcut(tr("Shift+Alt+t"));
    thresholdScalarsAct = new QAction(this);
    thresholdScalarsAct->setText(QApplication::translate("Model", "Threshold Scalar Values", 0, QApplication::UnicodeUTF8));
    thresholdScalarsAct->setShortcut(tr("Shift+Ctrl+t"));
    thresholdScalarsBinaryAct = new QAction(this);
    thresholdScalarsBinaryAct->setText(QApplication::translate("Model", "Binary Threshold Scalar Values", 0, QApplication::UnicodeUTF8));
    thresholdScalarsBinaryAct->setShortcut(tr("Shift+Ctrl+b"));
    maskScalarsAct = new QAction(this);
    maskScalarsAct->setText(QApplication::translate("Model", "Mask Scalar Values", 0, QApplication::UnicodeUTF8));
    maskScalarsAct->setShortcut(tr("Shift+Ctrl+m"));
    gradientAct = new QAction(this);
    gradientAct->setText(QApplication::translate("Model", "Scalar Gradient Field", 0, QApplication::UnicodeUTF8));
    gradientAct->setShortcut(tr("Shift+Alt+c"));
    removeScalarsAct = new QAction(this);
    removeScalarsAct->setText(QApplication::translate("Model", "&Remove All Scalars", 0, QApplication::UnicodeUTF8));
    removeScalarsAct->setShortcut(tr("Shift+Alt+x"));
    scalarsAct = new QAction(this);
    scalarsAct->setText(QApplication::translate("Model", "&Load Scalars from ...", 0, QApplication::UnicodeUTF8));
    scalarsAct->setShortcut(tr("Shift+Alt+s"));
    outScalarsAct = new QAction(this);
    outScalarsAct->setText(QApplication::translate("Model", "&Output Scalars as CSV File", 0, QApplication::UnicodeUTF8));
    outScalarsAct->setShortcut(tr("Shift+Alt+o"));
    scalarsGroup = new QActionGroup(this);
    scalarsGroup->addAction(renameScalarsAct);
    scalarsGroup->addAction(thresholdAct);
    scalarsGroup->addAction(thresholdScalarsAct);
    scalarsGroup->addAction(thresholdScalarsBinaryAct);
    scalarsGroup->addAction(gradientAct);
    scalarsGroup->addAction(removeScalarsAct);
    scalarsGroup->addAction(outScalarsAct);

    //Misc
    doAct = new QAction(this);
    doAct->setText(QApplication::translate("Model", "&Undo Last Processing", 0, QApplication::UnicodeUTF8));
    doAct->setShortcut(tr("Alt+u"));
    infoAct = new QAction(this);
    infoAct->setText(QApplication::translate("Model", "&Model Information", 0, QApplication::UnicodeUTF8));
    infoAct->setShortcut(tr("Alt+m"));
    textureAct = new QAction(this);
    textureAct->setText(QApplication::translate("Model", "&Extract Texture", 0, QApplication::UnicodeUTF8));
    textureAct->setShortcut(tr("Alt+t"));
    flipAct = new QAction(this);
    flipAct->setText(QApplication::translate("Model", "Flip Model", 0, QApplication::UnicodeUTF8));
    flipAct->setShortcut(tr("Ctrl+Alt+f"));
    rotateAct = new QAction(this);
    rotateAct->setText(QApplication::translate("Model", "Rotate Model", 0, QApplication::UnicodeUTF8));
    rotateAct->setShortcut(tr("Ctrl+Alt+r"));

    //View as
    pointsAct = new QAction(this);
    pointsAct->setText(QApplication::translate("Model", "&Points", 0, QApplication::UnicodeUTF8));
    pointsAct->setShortcut(tr("Alt+p"));
    pointsAct->setCheckable(true);
    wireframeAct = new QAction(this);
    wireframeAct->setText(QApplication::translate("Model", "&Wireframe", 0, QApplication::UnicodeUTF8));
    wireframeAct->setShortcut(tr("Alt+w"));
    wireframeAct->setCheckable(true);
    surfaceAct = new QAction(this);
    surfaceAct->setText(QApplication::translate("Model", "&Surface", 0, QApplication::UnicodeUTF8));
    surfaceAct->setShortcut(tr("Alt+s"));
    surfaceAct->setCheckable(true);
    displayGroup = new QActionGroup(this);
    displayGroup->addAction(pointsAct);
    displayGroup->addAction(wireframeAct);
    displayGroup->addAction(surfaceAct);

    colourAct = new QAction(this);
    colourAct->setText(QApplication::translate("Model", "Change Mesh &Colour/Transparency", 0, QApplication::UnicodeUTF8));
    colourAct->setShortcut(tr("Alt+c"));
    interpAct = new QAction(this);
    interpAct->setText(QApplication::translate("Model", "Best (Phong) Shading", 0, QApplication::UnicodeUTF8));
    interpAct->setShortcut(tr("Alt+g"));
    interpAct->setCheckable(true);
    interpAct->setChecked(false);

    normalsAct = new QAction(this);
    normalsAct->setText(QApplication::translate("Model", "Show Normals", 0, QApplication::UnicodeUTF8));
    normalsAct->setShortcut(tr("Alt+n"));
    normalsAct->setCheckable(true);
    normalsAct->setChecked(false);
    centroidAct = new QAction(this);
    centroidAct->setText(QApplication::translate("Model", "Show Centroid", 0, QApplication::UnicodeUTF8));
    centroidAct->setShortcut(tr("Alt+m"));
    centroidAct->setCheckable(true);
    centroidAct->setChecked(false);
    outlineAct = new QAction(this);
    outlineAct->setText(QApplication::translate("Model", "Show Outline Box", 0, QApplication::UnicodeUTF8));
    outlineAct->setShortcut(tr("Alt+b"));
    outlineAct->setCheckable(true);
    outlineAct->setChecked(false);
    cubeAxesAct = new QAction(this);
    cubeAxesAct->setText(QApplication::translate("Model", "Show Cube (Plot) Axes", 0, QApplication::UnicodeUTF8));
    cubeAxesAct->setShortcut(tr("Alt+a"));
    cubeAxesAct->setCheckable(true);
    cubeAxesAct->setChecked(false);
    overlaysAct = new QAction(this);
    overlaysAct->setText(QApplication::translate("Model", "Remove Overlays/Actors", 0, QApplication::UnicodeUTF8));
    overlaysAct->setShortcut(tr("Alt+o"));

    //Contour
    selectPointsAct = new QAction(this);
    selectPointsAct->setText(QApplication::translate("Model", "Select Points within Contour", 0, QApplication::UnicodeUTF8));
    weightAct = new QAction(this);
    weightAct->setText(QApplication::translate("Model", "Gaussian Weights from Contour", 0, QApplication::UnicodeUTF8));
}

void milxQtModel::createConnections()
{
    //Generation
    connect(genModelAct, SIGNAL(triggered()), this, SLOT(generateModel()));
    connect(genVerticesAct, SIGNAL(triggered()), this, SLOT(generateVertices()));
    connect(genNormalsAct, SIGNAL(triggered()), this, SLOT(generateNormals()));
    connect(genVectorsAct, SIGNAL(triggered()), this, SLOT(generateVectorField()));
    connect(genTensorsAct, SIGNAL(triggered()), this, SLOT(generateTensorField()));
    connect(genHedgehogAct, SIGNAL(triggered()), this, SLOT(generateHedgehog()));
    connect(genDelaunayAct, SIGNAL(triggered()), this, SLOT(generateDelaunayGraph()));
    connect(genDelaunayTri2DAct, SIGNAL(triggered()), this, SLOT(generateDelaunay2DTriangulation()));
    connect(genDelaunayTriAct, SIGNAL(triggered()), this, SLOT(generateDelaunayTriangulation()));
    connect(genPointModelAct, SIGNAL(triggered()), this, SLOT(generatePointModel()));
    connect(genLabelsAct, SIGNAL(triggered()), this, SLOT(generateLabels()));
    connect(genPointIDsAct, SIGNAL(triggered()), this, SLOT(generatePointIDsScalars()));
    connect(genSamplesAct, SIGNAL(triggered()), this, SLOT(generateSampledPoints()));
    connect(genKmeansAct, SIGNAL(triggered()), this, SLOT(generateKMeansClustering()));
    connect(genQuantiseAct, SIGNAL(triggered()), this, SLOT(generateQuantisedPoints()));
    connect(genCapBoundaries, SIGNAL(triggered()), this, SLOT(generateCappedBoundaries()));
    connect(genRegionLabels, SIGNAL(triggered()), this, SLOT(generateRegionLabels()));
    connect(genElevation, SIGNAL(triggered()), this, SLOT(generateElevation()));
    //~ connect(genReebGraph, SIGNAL(triggered()), this, SLOT(generateReebGraphs()));
    //Operations
    connect(doAct, SIGNAL(triggered()), this, SLOT(undoProcessing()));
    connect(cleanAct, SIGNAL(triggered()), this, SLOT(clean()));
    connect(triAct, SIGNAL(triggered()), this, SLOT(triangulate()));
    connect(decimateAct, SIGNAL(triggered()), this, SLOT(decimate()));
    connect(quadricDecimateAct, SIGNAL(triggered()), this, SLOT(quadricDecimate()));
    connect(clusterDecimateAct, SIGNAL(triggered()), this, SLOT(clusterDecimate()));
    connect(smoothAct, SIGNAL(triggered()), this, SLOT(smooth()));
    connect(smoothSincAct, SIGNAL(triggered()), this, SLOT(smoothSinc()));
    connect(curvatureAct, SIGNAL(triggered()), this, SLOT(curvature()));
    //Transformations
    connect(transformAct, SIGNAL(triggered()), this, SLOT(transform()));
    connect(matchAct, SIGNAL(triggered()), this, SLOT(matchInfo()));
    connect(registerAct, SIGNAL(triggered()), this, SLOT(registerICP()));
    connect(registerLandmarkAct, SIGNAL(triggered()), this, SLOT(registerLandmarks()));
    connect(voxeliseAct, SIGNAL(triggered()), this, SLOT(voxelise()));
    //Scalars
    connect(renameScalarsAct, SIGNAL(triggered()), this, SLOT(renameScalars()));
    connect(thresholdAct, SIGNAL(triggered()), this, SLOT(threshold()));
    connect(thresholdScalarsAct, SIGNAL(triggered()), this, SLOT(thresholdScalars()));
    connect(thresholdScalarsBinaryAct, SIGNAL(triggered()), this, SLOT(thresholdScalarsBinary()));
    connect(maskScalarsAct, SIGNAL(triggered()), this, SLOT(maskScalars()));
    connect(gradientAct, SIGNAL(triggered()), this, SLOT(gradient()));
    connect(removeScalarsAct, SIGNAL(triggered()), this, SLOT(removeScalars()));
    connect(scalarsAct, SIGNAL(triggered()), this, SLOT(loadScalars()));
    connect(outScalarsAct, SIGNAL(triggered()), this, SLOT(outputScalars()));
    //Misc
    connect(infoAct, SIGNAL(triggered()), this, SLOT(modelInfo()));
    connect(textureAct, SIGNAL(triggered()), this, SLOT(texture()));
    connect(flipAct, SIGNAL(triggered()), this, SLOT(flip()));
    connect(rotateAct, SIGNAL(triggered()), this, SLOT(rotate()));
    //View as
    connect(pointsAct, SIGNAL(triggered()), this, SLOT(generatePoints()));
    connect(wireframeAct, SIGNAL(triggered()), this, SLOT(generateWireframe()));
    connect(surfaceAct, SIGNAL(triggered()), this, SLOT(generateSurface()));
    //Properties
    connect(colourAct, SIGNAL(triggered()), this, SLOT(changeColour()));
    connect(interpAct, SIGNAL(triggered()), this, SLOT(toggleInterpolation()));
    //Show
    connect(normalsAct, SIGNAL(triggered()), this, SLOT(normals()));
    connect(centroidAct, SIGNAL(triggered()), this, SLOT(centroid()));
    connect(outlineAct, SIGNAL(triggered()), this, SLOT(outlineDisplay()));
    connect(cubeAxesAct, SIGNAL(triggered()), this, SLOT(cubeAxesDisplay()));
    connect(overlaysAct, SIGNAL(triggered()), this, SLOT(removeOverlays()));
    //Contouring
    connect(selectPointsAct, SIGNAL(triggered()), this, SLOT(selectPointsInContour()));
    connect(weightAct, SIGNAL(triggered()), this, SLOT(weightGaussianFromContour()));

//    connect(milxQtRenderWindow::backgroundAct, SIGNAL(triggered()), this, SLOT(background())); //Also connected to milxQtRenderWindow background member
    connect(milxQtRenderWindow::refreshAct, SIGNAL(triggered()), this, SLOT(refresh()));
    connect(milxQtRenderWindow::resetAct, SIGNAL(triggered()), this, SLOT(reset()));

    milxQtWindow::createConnections(); //consume right click events etc.
}

void milxQtModel::setupTooltips()
{
    doAct->setToolTip("Undo/Redo last processing computed.");
    doAct->setStatusTip("Undo/Redo last processing computed.");
    cleanAct->setToolTip("Clean the mesh removing and merging duplicate points");
    cleanAct->setStatusTip("Clean the mesh removing and merging duplicate points");
    triAct->setToolTip("Triangulate the mesh. Often a prerequisite to many algorithms");
    triAct->setStatusTip("Triangulate the mesh. Often a prerequisite to many algorithms");
    decimateAct->setToolTip("Reduce the number of polygon and vertices of the mesh by factor");
    decimateAct->setStatusTip("Reduce the number of polygon and vertices of the mesh by factor");
    smoothAct->setToolTip("Denoise the mesh using the standard Laplacian smoothing algorithm");
    smoothAct->setStatusTip("Denoise the mesh using the standard Laplacian smoothing algorithm");
    smoothSincAct->setToolTip("Denoise the mesh using the windowed sinc smoothing algorithm");
    smoothSincAct->setStatusTip("Denoise the mesh using the windowed sinc smoothing algorithm");
    transformAct->setToolTip("Affine transform the model using a ITK transform file");
    transformAct->setStatusTip("Affine transform the model using a ITK transform file");
    infoAct->setToolTip("Show the centroid, centroid size and covariance matrix of the model");
    infoAct->setStatusTip("Show the centroid, centroid size and covariance matrix of the model");
    matchAct->setToolTip("Match the information (centroid and scale) of another model to this model");
    matchAct->setStatusTip("Match the information (centroid and scale) of another model to this model");
    registerAct->setToolTip("Register or align this model to another model using Iterative Closest Point algorithm");
    registerAct->setStatusTip("Register or align this model to another model using Iterative Closest Point algorithm");
    voxeliseAct->setToolTip("Convert the surface to an image volume by voxelising it.");
    voxeliseAct->setStatusTip("Convert the surface to an image volume by voxelising it.");
    colourAct->setToolTip("Change the mesh colour and transparency (colour is only applicable for meshes with no scalars)");
    colourAct->setStatusTip("Change the mesh colour and transparency (colour is only applicable for meshes with no scalars)");
    interpAct->setToolTip("Change the interpolation of the mesh");
    interpAct->setStatusTip("Change the interpolation of the mesh");
    scaleAct->setToolTip("Show a colour bar of the mesh scalars");
    scaleAct->setStatusTip("Show a colour bar of the mesh scalars");
}

void milxQtModel::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    createMenu(contextMenu);

    contextMenu->exec(currentEvent->globalPos());
}

QMenu* milxQtModel::basicContextMenu()
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection
    contextMenu->setTitle(QApplication::translate("MainWindow", "Modelling", 0, QApplication::UnicodeUTF8));

    foreach(QAction *currAct, milxQtWindow::actionsToAdd)
    {
        contextMenu->addAction(currAct);
    }
    foreach(QMenu *currMenu, milxQtWindow::menusToAdd)
    {
        contextMenu->addMenu(currMenu);
    }
    contextMenu->addMenu(generationMenu());
    contextMenu->addMenu(operationsMenu());
    contextMenu->addMenu(transformsMenu());
    contextMenu->addMenu(scalarsMenu());
    contextMenu->addSeparator()->setText(tr("Options"));
    contextMenu->addAction(doAct);
    contextMenu->addAction(infoAct);
    contextMenu->addAction(textureAct);
    contextMenu->addAction(colourAct);
    contextMenu->addAction(interpAct);
    contextMenu->addSeparator()->setText(tr("Display"));
    contextMenu->addAction(pointsAct);
    contextMenu->addAction(wireframeAct);
    contextMenu->addAction(surfaceAct);
    contextMenu->addMenu(milxQtRenderWindow::viewMenu);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewXY);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZX);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::viewZY);
    milxQtRenderWindow::enableActionBasedOnView();
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::saveViewFileAct);
    milxQtRenderWindow::viewMenu->addAction(milxQtRenderWindow::loadViewFileAct);
    contextMenu->addMenu(milxQtRenderWindow::colourMapMenu);
    showMenu = contextMenu->addMenu("Show"); //!< Only exists for the duration of the context selection
    showMenu->addAction(normalsAct);
    showMenu->addAction(centroidAct);
    showMenu->addAction(outlineAct);
    showMenu->addAction(cubeAxesAct);
    showMenu->addAction(overlaysAct);
    showMenu->addAction(milxQtRenderWindow::humanAct);
    showMenu->addMenu(arraysMenu());

    if(GetScalars() && GetScalars()->GetNumberOfComponents() > 1 && vtkUnsignedCharArray::SafeDownCast(GetScalars()))
        milxQtRenderWindow::colourMapMenu->setDisabled(true);
    else if(!GetScalars() && !GetVectors() && !GetTensors())
        milxQtRenderWindow::colourMapMenu->setDisabled(true);
    else
        milxQtRenderWindow::colourMapMenu->setDisabled(false);

    return contextMenu;
}

QMenu* milxQtModel::generationMenu()
{
    generateMenu->addAction(genModelAct);
    generateMenu->addAction(genVerticesAct);
    generateMenu->addAction(genNormalsAct);
    generateMenu->addAction(genVectorsAct);
    generateMenu->addAction(genTensorsAct);
    generateMenu->addAction(genHedgehogAct);
    generateMenu->addAction(genDelaunayAct);
    generateMenu->addAction(genDelaunayTri2DAct);
    generateMenu->addAction(genDelaunayTriAct);
    generateMenu->addAction(genPointModelAct);
    generateMenu->addAction(genLabelsAct);
    generateMenu->addAction(genPointIDsAct);
    generateMenu->addAction(genSamplesAct);
    generateMenu->addAction(genKmeansAct);
    generateMenu->addAction(genQuantiseAct);
    generateMenu->addAction(genCapBoundaries);
    generateMenu->addAction(genRegionLabels);
    generateMenu->addAction(genElevation);
    generateMenu->addAction(genReebGraph);

    return generateMenu;
}

QMenu* milxQtModel::operationsMenu()
{
    operateMenu->addAction(cleanAct);
    operateMenu->addAction(triAct);
    operateMenu->addAction(decimateAct);
    operateMenu->addAction(quadricDecimateAct);
    operateMenu->addAction(clusterDecimateAct);
    operateMenu->addAction(smoothAct);
    operateMenu->addAction(smoothSincAct);
    operateMenu->addAction(curvatureAct);

    return operateMenu;
}

QMenu* milxQtModel::transformsMenu()
{
    transformMenu->addAction(flipAct);
    transformMenu->addAction(rotateAct);
    transformMenu->addAction(transformAct);
    transformMenu->addAction(matchAct);
    transformMenu->addAction(registerAct);
    transformMenu->addAction(registerLandmarkAct);
    transformMenu->addAction(voxeliseAct);

    return transformMenu;
}

QMenu* milxQtModel::scalarsMenu()
{
    scalarMenu->addAction(renameScalarsAct);
    scalarMenu->addAction(thresholdAct);
    scalarMenu->addAction(thresholdScalarsAct);
    scalarMenu->addAction(thresholdScalarsBinaryAct);
    scalarMenu->addAction(maskScalarsAct);
    scalarMenu->addAction(gradientAct);
    scalarMenu->addAction(removeScalarsAct);
    scalarMenu->addAction(scalarsAct);
    scalarMenu->addAction(outScalarsAct);

    if(!model.Result()->GetPointData()->GetScalars())
        scalarsGroup->setDisabled(true);
    else
        scalarsGroup->setDisabled(false);

    return scalarMenu;
}

QMenu* milxQtModel::arraysMenu()
{
    arrayGroup = new QActionGroup(this);
    arrayMenu = new QMenu(this);
    arrayMenu->setTitle(QApplication::translate("Model", "&Load Array as Scalars/Vectors", 0, QApplication::UnicodeUTF8));

    for(int j = 0; j < model.Result()->GetPointData()->GetNumberOfArrays(); j ++)
    {
        QAction *act = new QAction(this);
        act->setText(QApplication::translate("Model", model.Result()->GetPointData()->GetArrayName(j), 0, QApplication::UnicodeUTF8));
        arrayMenu->addAction(act);
        arrayGroup->addAction(act);
    }
    connect(arrayGroup, SIGNAL(triggered(QAction *)), this, SLOT(showArray(QAction *)));

    if(model.Result()->GetPointData()->GetNumberOfArrays() == 0)
        arrayMenu->setDisabled(true);
    else
        arrayMenu->setDisabled(false);

    return arrayMenu;
}
