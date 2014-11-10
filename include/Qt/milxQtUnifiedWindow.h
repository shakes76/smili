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
#ifndef MILXQTUNIFIEDWINDOW_H
#define MILXQTUNIFIEDWINDOW_H

#include <QList>
//VTK
#include <vtkCheckerboardWidget.h>
#include <vtkImageCheckerboard.h>

#include "milxQtImage.h"
#include "milxQtModel.h"

/**
    \class milxQtUnifiedWindow
    \brief The class maintains a state and render of multiple display objects (such as a milxQtModel or milxQtImage) in a unifying fashion.
    \author Shekhar S. Chandra, 2013

    This window supports operations on and displays of multiple types of data. Combining, differencing etc. of models and/or images is possible.
*/
class MILXQT_EXPORT milxQtUnifiedWindow : public milxQtRenderWindow
{
    Q_OBJECT

public:
    /**
        \fn milxQtUnifiedWindow::milxQtUnifiedWindow(QWidget *theParent = 0)
        \brief Default constructor
    */
    milxQtUnifiedWindow(QWidget *theParent = 0);
    /**
        \fn milxQtUnifiedWindow::~milxQtUnifiedWindow()
        \brief Default destructor
    */
    virtual ~milxQtUnifiedWindow();

public slots:
    /*!
        \fn milxQtUnifiedWindow::refresh()
        \brief Refresh the display of the window.
    */
    void refresh();

    /**
        \fn milxQtUnifiedWindow::addToWindow(milxQtModel *model)
        \brief Add this model window to the unified display.
    */
    void addToWindow(milxQtModel *model);
    /**
        \fn milxQtUnifiedWindow::addToWindow(milxQtImage *image)
        \brief Add this image window to the unified display.
    */
    void addToWindow(milxQtImage *image);

    /**
        \fn milxQtUnifiedWindow::removeFromWindow(QWidget *passedWindow)
        \brief Remove this window from the unified display.
    */
    void removeFromWindow(QWidget *passedWindow);

    /*!
        \fn milxQtUnifiedWindow::generateUnion()
        \brief Generate a union of the data to display, i.e. show all data at once unaltered.
    */
    void generateUnion();
    /*!
        \fn milxQtUnifiedWindow::generateDifference()
        \brief Generate a difference of the geometric data to display, i.e. show the net displacement vectors of the data.
    */
    void generateDifference(double pseudoInfinityFactor = -1.0);
    /*!
        \fn milxQtUnifiedWindow::generateScalarDifference()
        \brief Generate a difference of the scalar data to display, i.e. show the net displacement vectors of the data.
    */
    void generateScalarDifference();
    /*!
        \fn milxQtUnifiedWindow::generateCheckerBoard()
        \brief Generate a checker board widget to compare images, i.e. show checkboards with the compared images.
    */
    void generateCheckerBoard();

    /**
      \brief Function for extracting pixel values from image given surface points.
    */
    vtkSmartPointer<vtkFloatArray> surfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<floatImageType> img, const bool absoluteValues);

    /*!
        \fn milxQtUnifiedWindow::customOperation()
        \brief Custom operation, data dependent for viewing in unified environment.
    */
    virtual void customOperation() {};

signals:
    void imageAvailable(milxQtImage *);
    void modelAvailable(milxQtModel *);
    void modelAvailable(QWidget *);

protected:
    QList< milxQtModel* > unifyModels; //!< Model to maintain and operate on
    QList< milxQtImage* > unifyImages; //!< Images to maintain and operate on
    QList< vtkSmartPointer<vtkImageReslice> > slices; //!< Images to maintain and operate on

//    QPointer<milxQtModel> geoDiffModel; //!< Model maintaining the geometric difference of a vector field
    vtkSmartPointer<vtkCheckerboardWidget> checkerWidget; //!< Checker board widget for comparing images.
    vtkSmartPointer<vtkImageCheckerboard> checker;

    //States (as actions)
    QActionGroup* modeGroup; //!< Grouping for check boxes
    QAction* unionAct; //!< Union of data
    QAction* geoDifferenceAct; //!< Geometric Difference of data
    QAction* scalarDifferenceAct; //!< Scalar Difference of data
    QAction* checkerBoardAct; //!< Checkerboard of image data

    /*!
        \fn milxQtUnifiedWindow::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
        \fn milxQtUnifiedWindow::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
    	\fn milxQtUnifiedWindow::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);
    void setCommonProperties(QWidget *passedWindow);

private:

};

#endif // MILXQTUNIFIEDWINDOW_H
