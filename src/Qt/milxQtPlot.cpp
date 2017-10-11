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
#include "milxQtPlot.h"

#include <QMenu>
#include <QDialog>
#include <QInputDialog>
#include <QMessageBox>

#include <vtkContextView.h>
#include <vtkPlotPoints.h>
#include <vtkContextScene.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkWarpScalar.h>
#include <vtkPolyDataNormals.h>
#include <vtkAxis.h>
#include <vtkTextProperty.h>
//Volume Rendering
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>

milxQtPlot::milxQtPlot(QWidget *theParent, bool contextSystem) : milxQtModel(theParent, contextSystem)
{
    sourceLoaded = false;
    sourceEightbit = false;
    plotType2D = false;
    plotType3D = false;
    plotTypeSurface = false;
    plotTypeVolume = false;
    plotTypeGPU = false;
    plotted = false;
    xIndex = 0;
    yIndex = 1;
    zIndex = 2;
    tickNumber = 100;
    slideInterval = 1;

    table = vtkSmartPointer<vtkTable>::New();

    ///Set strings
    milxQtWindow::prefix = "Plot: ";

    ///Allocate critical aspects
    milxQtWindow::setDeletableOnClose(true);

    milxQtRenderWindow::humanAct->setChecked(false);
    milxQtRenderWindow::humanAct->setEnabled(false);
    humanGlyph->Off();

    createActions();
    createConnections();
}

milxQtPlot::~milxQtPlot()
{
    //dtor
}

void milxQtPlot::createMenu(QMenu *menu)
{
    cout << "Creating Plot Menu 1" << std::endl;
    if(!menu)
        return;

    menu->clear();
    if(plotTypeSurface || plotType3D)
        menu->addMenu(milxQtModel::basicContextMenu()); ///Have all the basic model options

    cout << "Creating Plot Menu 2" << std::endl;
    if(plotTypeVolume)
    {
      foreach(QAction *currAct, actionsToAdd)
      {
          menu->addAction(currAct);
      }
      foreach(QMenu *currMenu, menusToAdd)
      {
          menu->addMenu(currMenu);
      }
    }
    menu->addSeparator()->setText(tr("Plotting"));
    menu->addAction(xAxisName);
    if(plotType2D)
    {
        menu->addAction(titleName);
        menu->addAction(pointsAct);
        menu->addAction(logScaleAct);
        menu->addAction(legendAct);
    }
    if(!plotType2D)
    {
        menu->addAction(milxQtModel::outlineAct);
        menu->addAction(milxQtModel::cubeAxesAct);
        menu->addAction(milxQtModel::scaleAct);
        menu->addSeparator();
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
        menu->addMenu(milxQtRenderWindow::contourMenu);
        menu->addMenu(milxQtRenderWindow::windowPropertiesMenu);
    }
    menu->addAction(milxQtRenderWindow::refreshAct);
    if(!plotType2D)
        menu->addAction(milxQtRenderWindow::resetAct);
}

void milxQtPlot::generatePlot()
{
    if(!sourceLoaded)
        return;

    if(plotTypeVolume)
    {
        printInfo("Volume Plot mode set. Trying Volume plot.");
        volumePlot(imageData, sourceEightbit, plotTypeGPU);
        return;
    }
    if( plotTypeSurface && (table->GetNumberOfColumns() > 0 && !imageData) )
    {
        printInfo("Surface Plot mode set. Trying Table Surface plot.");
        surfacePlot(table);
        return;
    }
    if(plotType2D && table->GetNumberOfColumns() > 0)
    {
        printInfo("Scatter 2D Plot mode set. Trying Scatter xy plot.");
        scatterPlot(table, xIndex, yIndex);
        return;
    }
    if(plotType3D && table->GetNumberOfColumns() > 2)
    {
        printInfo("Scatter 3D Plot mode set. Trying Scatter xyz plot.");
        scatterPlot(table, xIndex, yIndex, zIndex);
        return;
    }
    else if(plotType3D && table->GetNumberOfColumns() < 3)
    {
        printError("Insufficient number of columns for 3D scatter plot. Try adding more columns.");
        return;
    }

    ///Not told anything, detect
    if(imageData)
    {
        printInfo("Image data found. Trying Image Surface plot.");
        surfacePlot(imageData, displaceAxis);
        return;
    }

    if(table->GetNumberOfColumns() == 1) //bar chart
    {
        cerr << "One column is not supported yet." << std::endl;
    }
    else if(table->GetNumberOfColumns() == 2) //2D scatter
    {
        printInfo("Two columns found. Trying 2D scatter plot.");
        plotType2D = true;
        scatterPlot(table, xIndex, yIndex);
    }
    else if(table->GetNumberOfColumns() > 3 && table->GetNumberOfColumns() < 10) //3D scatter
    {
        printInfo("Less than 10 columns found. Trying 3D scatter plot.");
        plotType3D = true;
        scatterPlot(table, xIndex, yIndex, zIndex);
    }
    else //surface plot
    {
        printInfo("More than 10 columns found. Trying Surface plot.");
        plotTypeSurface = true;
        surfacePlot(table);
    }
}

void milxQtPlot::scatterPlot(vtkSmartPointer<vtkTable> table, const int xColumn, const int yColumn)
{
    printDebug("Plotting... ");
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();

    templateChart = vtkSmartPointer<vtkChartXY>::New();
//        templateChart->SetRenderEmpty(true);
        templateChart->SetShowLegend(true);

    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
    view->GetRenderWindow()->SetSize(400, 300);
    view->GetScene()->AddItem(templateChart);
    //view->GetRenderWindow()->SetMultiSamples(0);
#if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION < 8)
    view->SetLabelRenderModeToQt();
#endif

    printInfo("Creating Scatter Plot ");
    vtkPlot *templatePoints = templateChart->AddPlot(vtkChart::LINE);
        printDebug("Source Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));
    #if VTK_MAJOR_VERSION <=5
        templatePoints->SetInput(table, xColumn, yColumn);
    #else
        templatePoints->SetInputData(table, xColumn, yColumn);
    #endif // VTK_MAJOR_VERSION
        templatePoints->SetColor(0, 0, 0, 255);
        templatePoints->SetWidth(1.0);
        vtkPlotPoints::SafeDownCast(templatePoints)->SetMarkerStyle(vtkPlotPoints::DIAMOND);
        //vtkPlotPoints::SafeDownCast(templatePoints)->SetMarkerStyle(vtkPlotPoints::CROSS);
        //vtkPlotPoints::SafeDownCast(templatePoints)->SetMarkerStyle(vtkPlotPoints::PLUS);
        templatePoints->Update();

    view->SetInteractor(QVTKWidget::GetInteractor()); //order important, also do not move up
    QVTKWidget::SetRenderWindow(view->GetRenderWindow());
    milxQtRenderWindow::SetRenderer(view->GetRenderer());
    milxQtRenderWindow::backgroundAct->setChecked(true);
    milxQtRenderWindow::generateRender();
    view->ResetCamera();
    view->Update();
    refresh();
    plotted = true;

    QString xAxisTitle(table->GetColumnName(xColumn));
    QString yAxisTitle(table->GetColumnName(yColumn));
    if(!xAxisTitle.isEmpty() && !yAxisTitle.isEmpty())
        renameAxes(xAxisTitle, yAxisTitle);
    renameTitle(this->name);
    legend(legendAct->isChecked());

    emit resultAvailable( qobject_cast<milxQtRenderWindow *>(this) ); //only using render class
}

void milxQtPlot::scatterPlot(vtkSmartPointer<vtkTable> table, const int xColumn, const int yColumn, const int zColumn)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(int j = 0; j < table->GetNumberOfRows(); j ++)
    {
        vtkVariant x = table->GetValue(j, xColumn);
        vtkVariant y = table->GetValue(j, yColumn);
        vtkVariant z = table->GetValue(j, zColumn);
        points->InsertNextPoint(x.ToDouble(), y.ToDouble(), z.ToDouble());
    }

    setName("Plot");
    SetPoints(points);
//    double bounds[6];
//    GetOutput()->GetBounds(bounds);
//    double scaling = ( (bounds[1]-bounds[0]) + (bounds[3]-bounds[2]) + (bounds[5]-bounds[4]) )/150;
//    if(points->GetNumberOfPoints() > 20000)
        generateVertices(0,0,0);
//        generateLabels();
//        generateVerticesAs(milxQtModel::Cross);
//    else
//        generatePointModel(scaling);
//    milxQtRenderWindow::backgroundAct->setChecked(true);
//    background(true);
    enableCubeAxes();
    printInfo("To generate larger points for scatter plot. Use Generate->Point Model in context menu.");
    plotted = true;

    QString xAxisTitle(table->GetColumnName(xColumn));
    QString yAxisTitle(table->GetColumnName(yColumn));
    QString zAxisTitle(table->GetColumnName(zColumn));
    printInfo("Column Names: " + xAxisTitle + ", " + yAxisTitle + ", " + zAxisTitle);
    if(!xAxisTitle.isEmpty() && !yAxisTitle.isEmpty() && !zAxisTitle.isEmpty())
        renameAxes(xAxisTitle, yAxisTitle, zAxisTitle);

    emit resultAvailable(this);
}

void milxQtPlot::surfacePlot(vtkSmartPointer<vtkImageData> img, const int zDirection)
{
    if(!imageData)
      imageData->DeepCopy(img);

    printDebug("Creating surface geometry");
    vtkSmartPointer<vtkImageDataGeometryFilter> geometry = vtkSmartPointer<vtkImageDataGeometryFilter>::New();
    #if VTK_MAJOR_VERSION <=5
        geometry->SetInput(img);
    #else
        geometry->SetInputData(img);
    #endif // VTK_MAJOR_VERSION
        geometry->Update();

    double range[2], bounds[6];
    img->GetScalarRange(range);
    img->GetBounds(bounds);
    const double scaling = 0.25*( 1.0 / ( 2.0*(range[1]-range[0]) / (bounds[1]+bounds[3]) ) );
    printInfo("Min/Max for z-axis: " + QString::number(range[0]) + "/" + QString::number(range[1]));
    printInfo("Scaling for z-axis: " + QString::number(scaling));

    printDebug("Warping based on scalars");
    vtkSmartPointer<vtkWarpScalar> warpScalar = vtkSmartPointer<vtkWarpScalar>::New();
        warpScalar->SetInputConnection(geometry->GetOutputPort());
        if(zDirection == 0) //x
            warpScalar->SetNormal(1, 0, 0);
        else if(zDirection == 1) //y
            warpScalar->SetNormal(0, 1, 0);
        else //z default
            warpScalar->SetNormal(0, 0, 1);
        warpScalar->UseNormalOn();
        warpScalar->SetScaleFactor(scaling);
        warpScalar->Update();

    ///Rewrite bounds since z axis is scaled for display
    printDebug("Rescaling z-axis for Display");
    double axesRange[6];
    axesRange[0] = bounds[0];
    axesRange[1] = bounds[1];
    axesRange[2] = bounds[2];
    axesRange[3] = bounds[3];
    axesRange[4] = range[0];
    axesRange[5] = range[1];

    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        normals->SetInputConnection(warpScalar->GetOutputPort());
        normals->Update();

    printDebug("Creating Surface Model for Display");
    setName("Surface Plot");
//        SetInputPointSet(warpScalar->GetOutput()); //doesnt require generate model call
    SetInput(normals->GetOutput()); //doesnt require generate model call
    generateModel();
    enableCubeAxes(axesRange);
//    enableCubeAxes();

    plotted = true;
    legendAct->setDisabled(true);
    emit resultAvailable(this);
}

void milxQtPlot::surfacePlot(vtkSmartPointer<vtkTable> table)
{
    int dim[3];
    dim[0] = table->GetNumberOfRows();
    dim[1] = table->GetNumberOfColumns();
    dim[2] = 0;

    imageData = vtkSmartPointer<vtkImageData>::New();
        imageData->SetDimensions(dim);
        imageData->SetExtent(0, table->GetNumberOfRows()-1, 0, table->GetNumberOfColumns()-1, 0, 0);
    #if VTK_MAJOR_VERSION <=5
        imageData->SetNumberOfScalarComponents(1);
        imageData->SetScalarTypeToDouble();
        imageData->AllocateScalars();
    #else
        imageData->AllocateScalars(VTK_DOUBLE,1);
    #endif

    ///Copy data into imagedata
    for(int j = 0; j < table->GetNumberOfRows(); j ++)
        for(int k = 0; k < table->GetNumberOfColumns(); k ++)
            imageData->SetScalarComponentFromDouble(j, k, 0, 0, table->GetValue(j, k).ToDouble());

    surfacePlot(imageData);
}

void milxQtPlot::volumePlot(vtkSmartPointer<vtkImageData> img, const bool eightbit, const bool quiet)
{
    if(!sourceLoaded)
        return;

    double range[2], bounds[6];
    img->GetScalarRange(range);
    img->GetBounds(bounds);

    //Create slicers and their dialogs
    QDialog *slidersDlg = new QDialog(this);
    QLabel *lowValueLbl = new QLabel(this);
    lowSldr = new QSlider(this);
    QLabel *highValueLbl = new QLabel(this);
    highSldr = new QSlider(this);
    QVBoxLayout *sliderLayout = new QVBoxLayout(this);
    slideInterval = (range[1]-range[0])/tickNumber;

    lowValueLbl->setText("Lower Value");
    lowSldr->setMinimum(0);
    lowSldr->setMaximum(tickNumber);
    lowSldr->setValue(0);
    highValueLbl->setText("Upper Value");
    highSldr->setMinimum(0);
    highSldr->setMaximum(tickNumber);
    highSldr->setValue(tickNumber);

    //Connect sliders
    connect(lowSldr, SIGNAL(valueChanged(int)), this, SLOT(updateVolumePlot(int)));
    connect(highSldr, SIGNAL(valueChanged(int)), this, SLOT(updateVolumePlot(int)));

    sliderLayout->addWidget(lowValueLbl);
    sliderLayout->addWidget(lowSldr);
    sliderLayout->addWidget(highValueLbl);
    sliderLayout->addWidget(highSldr);
    slidersDlg->setLayout(sliderLayout);
    if(!quiet)
        slidersDlg->show();

    //Construct rendering
    emit working(-1);
    generateRender();
    milxQtRenderWindow::generateRender(); // make sure we have an OpenGL context.

    printDebug("Colour Function Setup");
    vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
    if(opacityTransferValues.empty())
    {
        compositeOpacity->AddPoint(range[0], 0.0);
        compositeOpacity->AddPoint(range[1], 1.0);
    }
    else
    {
        for(size_t j = 0; j < opacityTransferValues.size(); j ++)
            compositeOpacity->AddPoint(colourTransferValues[j], opacityTransferValues[j]);
    }

    vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
    if(colourTransferValues.empty())
    {
        /*color->AddRGBPoint(range[0], 0.0,0.0,0.0);
        color->AddRGBPoint(range[1]/3, 1.0,0.0,0.0);
        color->AddRGBPoint(2*range[1]/3, 0.0,1.0,0.0);
        color->AddRGBPoint(range[1], 0.0,0.0,1.0);*/
        double incr = (range[1]-range[0])/100;
        for(double j = range[0]; j <= range[1]; j += incr)
        {
            ///Get Colour
            double colour[3];
            lookupTable->GetColor(j, colour); //!< Pull colour for data
            color->AddRGBPoint(j, colour[0], colour[1], colour[2]);
        }
    }
    else
    {
        /*vtkSmartPointer<vtkLookupTable> transferLookupTable = vtkSmartPointer<vtkLookupTable>::New();
//            transferLookupTable->SetTableRange(0.0, colourTransferValues.size());
            transferLookupTable->SetTableRange(range[0], range[1]);
            transferLookupTable->Build();*/
//        float minVolumeValue = numeric_limits<float>::max();
//        float maxVolumeValue = numeric_limits<float>::min();

        for(size_t j = 0; j < colourTransferValues.size(); j ++)
        {
            ///Get Colour
            double colour[3];
//            transferLookupTable->GetColor(j, colour); //!< Pull colour for data
//            transferLookupTable->GetColor(colourTransferValues[j], colour); //!< Pull colour for data
            lookupTable->GetColor(colourTransferValues[j], colour); //!< Pull colour for data

            color->AddRGBPoint(colourTransferValues[j], colour[0], colour[1], colour[2]);
        }
    }

    printDebug("Volume Property Setup");
    vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        volumeProperty->SetScalarOpacity(compositeOpacity);
        volumeProperty->SetColor(color);
//        volumeProperty->SetScalarOpacityUnitDistance(0.001);
        volumeProperty->ShadeOn();
        volumeProperty->SetAmbient(0.6);
        volumeProperty->SetDiffuse(0.4);
        volumeProperty->SetSpecular(0.2);
    if(eightbit)
        volumeProperty->SetInterpolationTypeToNearest();
    else
        volumeProperty->SetInterpolationTypeToLinear();

    printDebug("Volume Setup");
    volume = vtkSmartPointer<vtkVolume>::New();
        volume->SetProperty(volumeProperty);

    printDebug("Smart Volume Rendering Setup");
    vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
        volumeMapper->SetBlendModeToComposite(); // composite first
    #if VTK_MAJOR_VERSION <= 5
        volumeMapper->SetInputConnection(img->GetProducerPort());
    #else
        volumeMapper->SetInputData(img);
    #endif
    if(eightbit)
        volumeMapper->SetInterpolationModeToNearestNeighbor();
    else
        volumeMapper->SetInterpolationModeToLinear();
        volumeMapper->SetRequestedRenderModeToDefault(); //choose best automatically
//        volumeMapper->SetRequestedRenderModeToRayCastAndTexture();
//        volumeMapper->SetRequestedRenderModeToRayCast();
        linkProgressEventOf(volumeMapper);

    volume->SetMapper(volumeMapper);
    linkProgressEventOf(volume);
    volume->Update();

    AddVolume(volume);
    milxQtRenderWindow::renderer->ResetCamera();
//    if(colourTransferValues.empty())
        colourMapToJet();
    refresh();
    emit done(-1);

    enableCubeAxes(bounds, bounds);
    plotted = true;
    legendAct->setDisabled(true);
    emit resultAvailable( qobject_cast<milxQtRenderWindow *>(this) );
}

void milxQtPlot::renameAxes(QString xLabel, QString yLabel, QString zLabel)
{
    if(!plotted)
        return;

    if(xLabel.isEmpty() || yLabel.isEmpty())
    {
        bool ok1, ok2;
        xLabel = QInputDialog::getText(this, tr("Please Provide the name of the x-axis"),
                                              tr("x-Axis Name:"), QLineEdit::Normal, "x", &ok1);
        yLabel = QInputDialog::getText(this, tr("Please Provide the name of the y-axis"),
                                              tr("y-Axis Name:"), QLineEdit::Normal, "y", &ok2);

        if(!ok1 || !ok2)
            return;
    }

    if(zLabel.isEmpty() && (plotType3D || plotTypeSurface))
    {
        bool ok3;
        zLabel = QInputDialog::getText(this, tr("Please Provide the name of the z-axis"),
                                          tr("z-Axis Name:"), QLineEdit::Normal, "z", &ok3);

        if(!ok3)
            return;
    }

    if(plotType2D)
    {
        templateChart->GetPlot(0)->GetXAxis()->SetTitle(xLabel.toStdString().c_str());
        templateChart->GetPlot(0)->GetYAxis()->SetTitle(yLabel.toStdString().c_str());
    }
    else if(plotType3D || plotTypeSurface)
    {
        cubeAxesActor->SetXTitle(xLabel.toStdString().c_str());
        cubeAxesActor->SetYTitle(yLabel.toStdString().c_str());
        cubeAxesActor->SetZTitle(zLabel.toStdString().c_str());
    }
}

void milxQtPlot::renameTitle(QString newTitle)
{
    if(!plotted)
        return;

    if(newTitle.isEmpty())
    {
      bool ok1;
      newTitle = QInputDialog::getText(this, tr("Please Provide the name of the title"),
                                            tr("Title Name:"), QLineEdit::Normal, "", &ok1);

      if(!ok1)
          return;
    }

    if(plotType2D)
        templateChart->SetTitle(newTitle.toStdString().c_str());
//    else
//    {
//        cubeAxesActor->SetXTitle(xLabel.toStdString().c_str());
//        cubeAxesActor->SetYTitle(yLabel.toStdString().c_str());
//        cubeAxesActor->SetZTitle(zLabel.toStdString().c_str());
//    }
}

void milxQtPlot::legend(const bool showIt)
{
    if(!showIt)
        legendAct->setChecked(false);
    else
        legendAct->setDisabled(false);

    if(!plotted)
        return;

    if(plotType2D)
    {
        if(legendAct->isChecked())
            templateChart->SetShowLegend(true);
        else
            templateChart->SetShowLegend(false);
    }
}

void milxQtPlot::points(const bool showIt)
{
    if(!showIt)
        pointsAct->setChecked(false);
    else
        pointsAct->setDisabled(false);

    if(!plotted)
        return;

    if(plotType2D)
    {
        vtkPlot *plot = templateChart->GetPlot(0);
        if(pointsAct->isChecked())
            vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::DIAMOND);
        else
            vtkPlotPoints::SafeDownCast(plot)->SetMarkerStyle(vtkPlotPoints::NONE);
    }
}

void milxQtPlot::logScale(const bool showIt)
{
    if(!showIt)
        logScaleAct->setChecked(false);
    else
        logScaleAct->setDisabled(false);

    if(!plotted)
        return;

    if(plotType2D)
    {
        double range[2];
        vtkPlot *plot = templateChart->GetPlot(0);
        vtkDataArray *data = vtkDataArray::SafeDownCast(plot->GetInput()->GetColumn(1));
          data->GetRange(range);

        vtkAxis *axis = plot->GetYAxis();
        if(logScaleAct->isChecked())
        {
            axis->SetLogScale(true);
            #if ( (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION >= 8) || (VTK_MAJOR_VERSION > 5) )
              double scaleMin = log(range[0]);
              double scaleMax = log(range[1]);
              axis->SetRange(scaleMin, scaleMax);
            #endif
        }
        else
        {
            axis->SetLogScale(false);
            #if ( (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION >= 8) || (VTK_MAJOR_VERSION > 5) )
              double scaleMin = range[0];
              double scaleMax = range[1];
              axis->SetRange(scaleMin, scaleMax);
            #endif
        }
        //axis->RecalculateTickSpacing();
        //axis->AutoScale();
        axis->Update();
        templateChart->Update();
    }
}

void milxQtPlot::enableScale(QString title)
{
    if(!scale)
        scale = vtkSmartPointer<vtkScalarBarActor>::New();
    if(!scalarBar)
        scalarBar = vtkSmartPointer<vtkScalarBarWidget>::New();

    ///Ask if use default scalar LUT
    QMessageBox msgBox;
    msgBox.setText("An auto adjusted bar is about to be created");
    msgBox.setInformativeText("Would you like to customise the bar?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();

    const float barWidth = 0.1, barHeight = 0.7;

    if(ret == QMessageBox::Yes)
    {
        bool ok1, ok2, ok3, ok4;

        double minRange = QInputDialog::getDouble(this, tr("Enter Table Range of new Lookup Table"),
                                         tr("Minimum:"), 0, -2147483647, 2147483647, 5, &ok1);
        double maxRange = QInputDialog::getDouble(this, tr("Enter Table Range of new Lookup Table"),
                                         tr("Maximum:"), 1.0, -2147483647, 2147483647, 5, &ok2);
        int noOfLabels = QInputDialog::getInt(this, tr("How many labels to show"),
                                          tr("Labels:"), 3, 0, 99, 1, &ok3);
        title = QInputDialog::getText(this, tr("Title of Bar"),
                                          tr("Title:"), QLineEdit::Normal,
                                          title, &ok4);

        if(!ok1 || !ok2 || !ok3 || !ok4)
            return;

        scale->SetLookupTable(lookupTable);
        scale->SetNumberOfLabels(noOfLabels);
        if(plotTypeSurface)
        {
            modelMapper->SetLookupTable(lookupTable);
            modelMapper->SetScalarRange( minRange, maxRange );
        }
        else if(plotTypeVolume)
        {
            slideInterval = (maxRange-minRange)/tickNumber;
            vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = volume->GetProperty()->GetScalarOpacity();
                compositeOpacity->RemoveAllPoints();
                compositeOpacity->AddPoint(minRange, 0.0);
                compositeOpacity->AddPoint(maxRange, 0.2);

//            vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
//                color->DeepCopy(lookupTable.GetPointer());

//            volume->GetProperty()->SetColor(color);
            volume->Update();
        }
        customScalarBar = true;
    }
    else
    {
        if(plotTypeSurface)
            modelMapper->SetScalarRange( model.Result()->GetScalarRange() );
        else if(plotTypeVolume)
        {
            double range[2];
            imageData->GetScalarRange(range);

            slideInterval = (range[1]-range[0])/tickNumber;
            vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = volume->GetProperty()->GetScalarOpacity();
                compositeOpacity->RemoveAllPoints();
                compositeOpacity->AddPoint(range[0], 0.0);
                compositeOpacity->AddPoint(range[1], 0.2);

//            vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
//                color->DeepCopy(lookupTable.GetPointer());

//            volume->GetProperty()->SetColor(color);
            volume->Update();
        }
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

void milxQtPlot::scaleDisplay()
{
    if(plotTypeSurface)
        milxQtModel::scaleDisplay();
    else if(plotTypeVolume)
    {
      if(scaleAct->isChecked())
          enableScale("Voxels");
      else
          disableScale();

      Render();
    }
}

void milxQtPlot::updateLookupTable()
{
    double range[2];
    GetDataSet()->GetScalarRange(range); //This will propagate upwards to get the range for images or meshes

    if(plotTypeSurface)
        milxQtModel::updateLookupTable();
    else if(plotTypeVolume)
    {
        vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
        if(colourTransferValues.empty())
        {
            double incr = (range[1]-range[0])/100;
            for(double j = range[0]; j <= range[1]; j += incr)
            {
                ///Get Colour
                double colour[3];
                lookupTable->GetColor(j, colour); //!< Pull colour for data
                color->AddRGBPoint(j, colour[0], colour[1], colour[2]);
            }
        }
        else
        {
            for(size_t j = 0; j < colourTransferValues.size(); j ++)
            {
                ///Get Colour
                double colour[3];
                lookupTable->GetColor(colourTransferValues[j], colour); //!< Pull colour for data
                color->AddRGBPoint(colourTransferValues[j], colour[0], colour[1], colour[2]);
            }
        }

        volume->GetProperty()->SetColor(color);
        volume->Update();

        if(scaleAct->isChecked())
            enableScale("Voxels");
        else
            disableScale();

        Render();
    }
}

void milxQtPlot::updateCoords(vtkObject *obj)
{
    if(plotTypeSurface || plotType3D)
    {
        milxQtModel::updateCoords(obj);
        return;
    }

    ///Get interactor
//    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    QString message = "Coordinates Unsupported at the moment.";
    /*
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
            double *coords = currentMesh->GetPoint(nVolIdx);

            if(currentMesh->GetPointData()->GetScalars())
            {
                double *scalar = currentMesh->GetPointData()->GetScalars()->GetTuple(nVolIdx);

                if(currentMesh->GetPointData()->GetScalars()->GetNumberOfComponents() == 3)
                    message = "Point " + QString::number(nVolIdx) + ": (" + QString::number(coords[0]) + ", " + QString::number(coords[1]) + ", " + QString::number(coords[2]) + ") = "
                              + "[" + QString::number(scalar[0]) + ", " + QString::number(scalar[1]) + ", " + QString::number(scalar[2]) + "]";
                else if(currentMesh->GetPointData()->GetScalars()->GetNumberOfComponents() == 2)
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
    }*/

    ///Write message to status bar
    updateBar->showMessage(message);
}

void milxQtPlot::updateVolumePlot(int value)
{
    const float value1 = lowSldr->minimum() + lowSldr->value()*slideInterval;
    const float value2 = highSldr->minimum() + highSldr->value()*slideInterval;
    printDebug("Slider Min Value: " + QString::number(value1));
    printDebug("Slider Max Value: " + QString::number(value2));

    vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = volume->GetProperty()->GetScalarOpacity();
        compositeOpacity->RemoveAllPoints();
        compositeOpacity->AddPoint(value1, 0.0);
        compositeOpacity->AddPoint(value2, 0.2);

    volume->Update();
    Render();
}

void milxQtPlot::createActions()
{
    //axes
    xAxisName = new QAction(this);
    xAxisName->setText(tr("Rename &Axes", 0));
    xAxisName->setShortcut(tr("Alt+x"));
    titleName = new QAction(this);
    titleName->setText(tr("Rename &Title", 0));
    titleName->setShortcut(tr("Alt+t"));
    legendAct = new QAction(this);
    legendAct->setText(tr("Legend", 0));
    legendAct->setShortcut(tr("Alt+l"));
    legendAct->setCheckable(true);
    legendAct->setChecked(true);
    pointsAct = new QAction(this);
    pointsAct->setText(tr("Show Points", 0));
    pointsAct->setShortcut(tr("Shift+Alt+l"));
    pointsAct->setCheckable(true);
    pointsAct->setChecked(true);
    logScaleAct = new QAction(this);
    logScaleAct->setText(tr("Log Scale", 0));
    logScaleAct->setShortcut(tr("Alt+s"));
    logScaleAct->setCheckable(true);
    logScaleAct->setChecked(false);
    logScaleAct->setDisabled(true); ///\todo Logscale in 2D plot not working. fix.
}

void milxQtPlot::createConnections()
{
    //axes
    connect(xAxisName, SIGNAL(triggered()), this, SLOT(renameAxes()));
    connect(titleName, SIGNAL(triggered()), this, SLOT(renameTitle()));
    connect(pointsAct, SIGNAL(triggered()), this, SLOT(points()));
    connect(logScaleAct, SIGNAL(triggered()), this, SLOT(logScale()));
    connect(legendAct, SIGNAL(triggered()), this, SLOT(legend()));

    connect(milxQtRenderWindow::refreshAct, SIGNAL(triggered()), this, SLOT(refresh()));
}

void milxQtPlot::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection
    contextMenu->setTitle(tr("MainWindow", "Plotting", 0));

    createMenu(contextMenu);

    contextMenu->exec(currentEvent->globalPos());
}
