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
#ifndef MILXQTPLOT_H
#define MILXQTPLOT_H

//VTK
#include <vtkChartXY.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolume.h>

#include "milxQtModel.h"

/**
    \class milxQtPlot
    \brief The class provides 1D/2D and 3D plotting capability. This includes scatter and surface plots.
    \author Shekhar S. Chandra, 2013

    Usage:
    \code
    //Example to use plot that displays itself
    vtkSmartPointer<vtkTable> table;
    QPointer<milxQtPlot> plot = new milxQtPlot;
    success = milx::File::OpenDelimitedText(filename.toStdString(), table);

    plot->setPlotType2D(xCol, yCol);
    plot->setName(filename);
    plot->setConsole(console);
    plot->SetSource(table);
    plot->generatePlot();
    plot->show();
    \endcode

    Usage with connections:
    \code
    //Exmaple to use plot that triggers plot views
    vtkSmartPointer<vtkTable> table;
    QPointer<milxQtPlot> plot = new milxQtPlot;
    success = milx::File::OpenDelimitedText(filename.toStdString(), table);

    //Connect for display
    connect(plot, SIGNAL(resultAvailable(milxQtModel*)), this, SLOT(display(milxQtModel*)));
    connect(plot, SIGNAL(resultAvailable(milxQtRenderWindow*)), this, SLOT(display(milxQtRenderWindow*)));

    plot->setPlotType2D(xCol, yCol);
    plot->setName(filename);
    plot->setConsole(console);
    plot->SetSource(table);
    plot->generatePlot(); //emits result
    \endcode
*/
class MILXQT_EXPORT milxQtPlot : public milxQtModel
{
    Q_OBJECT

public:
    /**
        \fn milxQtPlot::milxQtPlot(QWidget *theParent = 0, bool contextSystem = true)
        \brief Default constructor
    */
    milxQtPlot(QWidget *theParent = 0, bool contextSystem = true);
    /**
        \fn milxQtPlot::~milxQtPlot()
        \brief Default destructor
    */
    virtual ~milxQtPlot();

    void setPlotType2D(const int indexX = 0, const int indexY = 1)
    {
        plotType2D = true;
        plotType3D = false;
        plotTypeSurface = false;
        plotTypeVolume = false;
        xIndex = indexX;
        yIndex = indexY;
    }
    void setPlotType3D(const int indexX = 0, const int indexY = 1, const int indexZ = 2)
    {
        plotType2D = false;
        plotType3D = true;
        plotTypeSurface = false;
        plotTypeVolume = false;
        xIndex = indexX;
        yIndex = indexY;
        zIndex = indexZ;
    }
    inline void setPlotTypeSurface(const int zAxis = 2)
    {
        plotTypeSurface = true;
        plotTypeVolume = false;
        plotType2D = false;
        plotType3D = false;
        displaceAxis = zAxis;
    }
    inline void setPlotTypeVolume()
    {
        plotTypeVolume = true;
        plotTypeSurface = false;
        plotType2D = false;
        plotType3D = false;
    }

    /*!
        \fn milxQtPlot::SetSource(vtkSmartPointer<vtkImageData> img, const bool eightbit, const bool useGPU = false)
        \brief Assigns the image provided to the class, preparing for plot. Call generatePlot() and then show() to display.
    */
    void SetSource(vtkSmartPointer<vtkImageData> img, const bool eightbit, const bool useGPU = false)
    {
        imageData = vtkSmartPointer<vtkImageData>::New();
            imageData->DeepCopy(img);
        sourceLoaded = true;
        sourceEightbit = eightbit;
        plotTypeGPU = useGPU;
    }
    /*!
        \fn milxQtPlot::SetSource(vtkSmartPointer<vtkTable> tbl)
        \brief Assigns the table provided to the class, preparing for plot. Call generatePlot() and then show() to display.
    */
    void SetSource(vtkSmartPointer<vtkTable> tbl)
    {
        table->DeepCopy(tbl);
        printDebug("Source Table is " + QString::number(table->GetNumberOfRows()) + "x" + QString::number(table->GetNumberOfColumns()));
        sourceLoaded = true;
    }
    /*!
        \fn milxQtPlot::GetSource()
        \brief Returns the current table source for plot held.
    */
    inline vtkSmartPointer<vtkTable> GetSource()
    {
        return table;
    }
    /*!
        \fn milxQtPlot::GetSourceImage()
        \brief Returns the current image source for plot held.
    */
    inline vtkSmartPointer<vtkImageData> GetSourceImage()
    {
        return imageData;
    }
    /*!
      \brief Get the underlying data, return the vtkImageData object if volume plot else polydata, useful for getting scalar range etc.
    */
    inline virtual vtkDataSet* GetDataSet()
    {
        if(plotTypeSurface)
            return model.Result();

        return imageData;
    }

    inline void setOpacityTransferValues(std::vector<float> values)
    {
        opacityTransferValues = values;
    }
    inline std::vector<float> getOpacityTransferValues()
    {
        return opacityTransferValues;
    }
    inline void setColourTransferValues(std::vector<float> values)
    {
        colourTransferValues = values;
    }
    inline std::vector<float> getColourTransferValues()
    {
        return colourTransferValues;
    }

public slots:
    /*!
        \fn milxQtPlot::refresh()
        \brief Refresh the display of the window.
    */
    inline void refresh()
    {
        milxQtRenderWindow::refresh();
    }
    /*!
        \fn milxQtPlot::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);

    /**
        \fn milxQtPlot::generatePlot()
        \brief Generate appropriate plot from the table of data provided by SetInput().
    */
    void generatePlot();

    /**
        \fn milxQtPlot::scatterPlot(vtkSmartPointer<vtkTable> table, const int xColumn, const int yColumn)
        \brief Display a scatter x-y plot from data in table given at column indices x and y.
    */
    void scatterPlot(vtkSmartPointer<vtkTable> table, const int xColumn, const int yColumn);
    /**
        \fn milxQtPlot::scatterPlot(vtkSmartPointer<vtkTable> table, const int xColumn, const int yColumn, const int zColumn)
        \brief Display a 3D scatter x-y plot from data in table given at column indices x, y and z.
    */
    void scatterPlot(vtkSmartPointer<vtkTable> table, const int xColumn, const int yColumn, const int zColumn);
    /**
        \fn milxQtPlot::surfacePlot(vtkSmartPointer<vtkImageData> img, const int zDirection = 2)
        \brief Display a surface plot from image data given.

        zDirection is the axis number of the proposed height/displacement axis of the surface (0 - x, 1 - y, 2 - z). The default is z (2).
    */
    void surfacePlot(vtkSmartPointer<vtkImageData> img, const int zDirection = 2);
    /**
        \fn milxQtPlot::surfacePlot(vtkSmartPointer<vtkTable> table)
        \brief Display a surface plot from data in table given. Table is converted to image data and then plotted.
    */
    void surfacePlot(vtkSmartPointer<vtkTable> table);
    /**
        \fn milxQtPlot::volumePlot(vtkSmartPointer<vtkImageData> img, const bool eightbit, const bool quiet)
        \brief Display a volume plot from data in image given.

        Hardware and datatype should be automatically handled.
    */
    void volumePlot(vtkSmartPointer<vtkImageData> img, const bool eightbit, const bool quiet = false);

    /**
        \fn milxQtPlot::renameAxes()
        \brief Rename the axes
    */
    void renameAxes(QString xLabel = "", QString yLabel = "", QString zLabel = "");
    /**
        \fn milxQtPlot::renameTitle()
        \brief Rename the title of the plot
    */
    void renameTitle(QString newTitle = "");
    /**
        \fn milxQtPlot::points(const bool showIt = false)
        \brief Show/hide points for the line plot
    */
    void points(const bool showIt = true);
    /**
        \fn milxQtPlot::logScale(const bool showIt = false)
        \brief Show/hide log scale (of the y-axis) for the line plot
    */
    void logScale(const bool showIt = true);
    /**
        \fn milxQtPlot::legend(const bool showIt = false)
        \brief Show/hide legend for plot
    */
    void legend(const bool showIt = true);
    /*!
        \fn milxQtPlot::enableScale(QString title = "")
        \brief Enable scale bar display with the title provided.
    */
    void enableScale(QString title = "");
    /*!
        \fn milxQtPlot::scaleDisplay()
        \brief Toggles the scale bar display.
    */
    void scaleDisplay();
    /*!
        \fn milxQtPlot::updateCoords(vtkObject *obj)
        \brief Picks the coordinates and pixel value from the current mouse position in the window.
    */
    virtual void updateCoords(vtkObject *obj);
    /*!
        \fn milxQtPlot::updateLookupTable()
        \brief Sets the necessary LUTs to model view.
    */
    virtual void updateLookupTable();
    /*!
        \fn milxQtPlot::updateVolumePlot(int value)
        \brief Updates the volume rendering display. Used when slider values are changed.
    */
    virtual void updateVolumePlot(int value);

signals:

protected:
    bool sourceLoaded; //!< Table/data loaded in class?
    bool sourceEightbit; //!< Table/data loaded is eightbit?
    bool plotType2D; //!< 2D plot type?
    bool plotType3D; //!< 3D plot type?
    bool plotTypeSurface; //!< 2D plot type?
    bool plotTypeVolume; //!< 2D plot type?
    bool plotTypeGPU; //!< Use GPU whenever possible?
    bool plotted; //!< plot done?
    int xIndex; //!< x axis index of table
    int yIndex; //!< y axis index of table
    int zIndex; //!< z axis index of table
    int displaceAxis; //!< warp axis for surface plots

    vtkSmartPointer<vtkTable> table; //!< data presented as a table
    vtkSmartPointer<vtkImageData> imageData; //!< data presented as a image
    vtkSmartPointer<vtkChartXY> templateChart; //!< chart (used only for xy scatter plot)
    vtkSmartPointer<vtkVolume> volume; //! Maintains volume rendering when needed
    vtkSmartPointer<vtkSmartVolumeMapper> volumeFixedMapper; //!< Volume rendering member

    //------------------
    //Context Menu
    QAction* xAxisName; //!< Action for x-axis name
    QAction* titleName; //!< Action for title name
    QAction* pointsAct; //!< Action for showing line points
    QAction* logScaleAct; //!< Action for showing a log scale
    QAction* legendAct; //!< Action for showing legend

    std::vector<float> opacityTransferValues; //!< values for volume rendering opacity transfer function
    std::vector<float> colourTransferValues; //!< values for volume rendering colour transfer function

    /*!
        \fn milxQtPlot::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
        \fn milxQtPlot::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
    	\fn milxQtPlot::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);

private:
    size_t tickNumber; //!< Number of ticks in the sliders
    float slideInterval; //!< interval for any sliders used, as sliders are integer-only

    QSlider *lowSldr;
    QSlider *highSldr;
};

#endif // MILXQTPLOT_H
