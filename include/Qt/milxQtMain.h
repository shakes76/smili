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
#ifndef MILXQTMAIN_H
#define MILXQTMAIN_H
//Qt
#include <QPointer>
#include <QMainWindow>
#include <QTabWidget>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QList>
#include <QSlider>
#include <QPushButton>
#include <QProgressBar>
#include <QComboBox>
#include <QCheckBox>
#include <QWebEngineView>

//VTK
#include <vtkEventQtSlotConnect.h>
//Displays
#include "milxQtImage.h"
#include "milxQtModel.h"
#include "milxQtUnifiedWindow.h"
#include "milxQtPluginInterface.h"

/*!
    \class milxQtMain
    \brief This class represents the MILX Qt Main Window object using Qt.
    \author Shekhar S. Chandra, 2013

    The class is able to display milxQtImage, milxQtModel and other classes into tabbed workspaces. This forms the main user interface of the library and
    wraps the overall usage of the display classes.

    The rendering is encapsulated within a Qt QMainWindow widget.

    Usage Example:
    \code
    QPointer<milxQtMain> MainWindow = new milxQtMain;

    ///Create model from a sphere source thats in VTK
    QPointer<milxQtModel> model = new milxQtModel;
        model->setName("Original Model");
        model->SetInput(sphere->GetOutput());
        model->generateModel();
        model->generateLabels();

    ///Create matrix
    vnl_matrix<double> rotationMatrix = doSomething();

    ///Image viewer
    QPointer<milxQtImage> imageViewer = new milxQtImage;
        imageViewer->setName("Rotation Matrix");
        imageViewer->setData(rotationMatrix);
        imageViewer->generateImage();

    ///Add displays to main window
    MainWindow->display(model); //Model
        MainWindow->display(imageViewer); //Image
        MainWindow->show();
    \endcode

    One can also iterate through the all the windows in the current tab as the following example:
    \code
    MainWindow.loadFiles(files); ///Load all the models

    MainWindow.initialiseWindowTraversal();
    for(int j = 0; j < MainWindow.getNumberOfWindows(); j ++)
    {
        QPointer<milxQtModel> mdl = MainWindow.nextModel();

        if(mdl)
            procrustes->SetInput(j, mdl->GetOutput()); //Do something
    }
    \endcode
*/
class MILXQT_EXPORT milxQtMain : public QMainWindow
{
    Q_OBJECT

public:
    /*!
        \fn milxQtMain::milxQtMain(QWidget *parent = 0)
        \brief The standard constructor
    */
    milxQtMain(QWidget *theParent = 0);
    /*!
        \fn milxQtMain::~milxQtMain()
        \brief The standard destructor
    */
    virtual ~milxQtMain();

    /*!
        \fn milxQtMain::activeName()
        \brief Returns the name of the active window regardless of being an image, table or surface plot.
    */
    QString activeName();
    /*!
        \fn milxQtMain::activeNamePrefix()
        \brief Returns the name with prefix of the active window regardless of being an image, table or surface plot.
    */
    QString activeNamePrefix();

    //Options
    /*!
        \fn milxQtMain::preferWhiteBackground(const bool whiteBack)
        \brief Sets whether white background is to be used whenever possible. Default: false
    */
    inline void preferWhiteBackground(const bool whiteBack)
    {
        whiteBackground = whiteBack;
    }
    inline bool isWhiteBackgroundPreferred()
    {
        return whiteBackground;
    }
    /*!
        \fn milxQtMain::preferHumanGlyph(const bool human)
        \brief Sets whether human orientation glyph is to be shown whenever possible. Default: true
    */
    inline void preferHumanGlyph(const bool human)
    {
        humanGlyph = human;
    }
    inline bool isHumanGlyphPreferred()
    {
        return humanGlyph;
    }
    /*!
        \fn milxQtMain::preferSubWindowSize(const int winSize)
        \brief Sets number of sub-window size is to be used whenever possible. Default: double of minWindowSize (around 512).
    */
    inline void preferSubWindowSize(const int winSize)
    {
        subWindowSize = winSize;
    }
    inline int hasPreferredSubWindowSize()
    {
        return subWindowSize;
    }
    /*!
        \fn milxQtMain::preferMaximumProcessors(const int procs)
        \brief Sets number of max processors is to be used whenever possible. Default: Half of max available on CPU to a minimum of one.
    */
    inline void preferMaximumProcessors(const int procs)
    {
        maxProcessors = procs;
    }
    inline int hasMaximumProcessors()
    {
        return maxProcessors;
    }
    /*!
        \fn milxQtMain::preferScreenshotMagnifyFactor(const int factor)
        \brief Sets the magnify factor when saving screenshots. Default: 2.
    */
    inline void preferScreenshotMagnifyFactor(const int factor)
    {
        magnifyFactor = factor;
    }
    inline int hasScreenshotMagnifyFactor()
    {
        return magnifyFactor;
    }
    /*!
        \fn milxQtMain::preferTimestamps(const bool timestamp)
        \brief Sets whether timestamps are to be shown whenever possible. Default: true
    */
    inline void preferTimestamps(const bool timestamps)
    {
        timestamping = timestamps;
    }
    inline bool isTimestampsPreferred()
    {
        return timestamping;
    }
	/*!
	\fn milxQtMain::setCurrentTheme(const QString theme)
	\brief Sets whether timestamps are to be shown whenever possible. Default: true
	*/
	inline void setCurrentTheme(const QString theme)
	{
		appTheme = theme;
	}
	inline QString currentTheme()
	{
		return appTheme;
	}
	/*!
        \fn milxQtMain::preferImageInterpolation(const bool interp)
        \brief Sets whether interpolation is to be shown whenever possible for images. Default: true
    */
    inline void preferImageInterpolation(const bool interp)
    {
        interpolationImages = interp;
    }
    inline bool isImageInterpolationPreferred()
    {
        return interpolationImages;
    }
    /*!
        \fn milxQtMain::preferOrientation(const bool orient)
        \brief Sets whether image orientation is to be applied whenever possible for images. Default: true
    */
    inline void preferOrientation(const bool orient)
    {
        orientationImages = orient;
    }
    inline bool isOrientationPreferred()
    {
        return orientationImages;
    }
    /*!
        \fn milxQtMain::preferModelInterpolation(const bool interp)
        \brief Sets whether interpolation is to be shown whenever possible for models. Default: false
    */
    inline void preferModelInterpolation(const bool interp)
    {
        interpolationModels = interp;
    }
    inline bool isModelInterpolationPreferred()
    {
        return interpolationModels;
    }
    /*!
        \fn milxQtMain::preferScalarBar(const bool bar)
        \brief Sets whether scalar bar is to be show whenever possible for models. Default: false
    */
    inline void preferScalarBar(const bool bar)
    {
        scalarBarModels = bar;
    }
    inline bool isScalarBarPreferred()
    {
        return scalarBarModels;
    }
	/*!
	\fn milxQtMain::enableColourMaps(const bool enable)
	\brief Sets whether to use custom colour maps. Default: true
	*/
	inline void activateColourMaps(const bool enable)
	{
		enableColourMaps = enable;
	}
	inline bool isColourMapEnabled()
	{
		return enableColourMaps;
	}

public slots:
    /*!
        \fn milxQtMain::isRender(QWidget *win)
        \brief Returns true if window is a generic render (milxQtRenderWindow object)
    */
    bool isRender(QWidget *win);
    /*!
        \fn milxQtMain::isActiveRender()
        \brief Returns true if active window is a generic render (milxQtRenderWindow object)
    */
    bool isActiveRender();
    /*!
        \fn milxQtMain::activeRender()
        \brief Returns the milxQtRenderWindow object, returns 0 if active window is not a milxQtRenderWindow object
    */
    milxQtRenderWindow* activeRender();
    /*!
        \fn milxQtMain::isImage(QWidget *win)
        \brief Returns true if window is an Image (milxQtImage object)
    */
    bool isImage(QWidget *win);
    /*!
        \fn milxQtMain::isActiveImage()
        \brief Returns true if active window is an Image (milxQtImage object)
    */
    bool isActiveImage();
    /*!
        \fn milxQtMain::activeImage()
        \brief Returns the milxQtImage object, returns 0 if active window is not a milxQtImage object
    */
    milxQtImage* activeImage();
    /*!
        \fn milxQtMain::isModel(QWidget *win)
        \brief Returns true if window is an Model (milxQtModel object)
    */
    bool isModel(QWidget *win);
    /*!
        \fn milxQtMain::isActiveModel()
        \brief Returns true if active window is an Model (milxQtModel object)
    */
    bool isActiveModel();
    /*!
        \fn milxQtMain::activeModel()
        \brief Returns the milxQtModel object, returns 0 if active window is not a milxQtModel object
    */
    milxQtModel* activeModel();
    /*!
        \fn milxQtMain::isPlot(QWidget *win)
        \brief Returns true if window is an Plot (milxQtPlot object)
    */
    bool isPlot(QWidget *win);
    /*!
        \fn milxQtMain::isActivePlot()
        \brief Returns true if active window is an Plot (milxQtPlot object)
    */
    bool isActivePlot();
    /*!
        \fn milxQtMain::activePlot()
        \brief Returns the milxQtPlot object, returns 0 if active window is not a milxQtPlot object
    */
    milxQtModel* activePlot();
    /*!
        \fn milxQtMain::isUnifiedWindow(QWidget *win)
        \brief Returns true if window is an UnifiedWindow (milxQtUnifiedWindow object)
    */
    bool isUnifiedWindow(QWidget *win);
    /*!
        \fn milxQtMain::isActiveUnifiedWindow()
        \brief Returns true if active window is an UnifiedWindow (milxQtUnifiedWindow object)
    */
    bool isActiveUnifiedWindow();
    /*!
        \fn milxQtMain::activeUnifiedWindow()
        \brief Returns the milxQtUnifiedWindow object, returns 0 if active window is not a milxQtUnifiedWindow object
    */
    milxQtUnifiedWindow* activeUnifiedWindow();
    /*!
        \fn milxQtMain::isActiveWebView()
        \brief Returns true if active window is a web viewer
    */
    bool isActiveWebView();
    /*!
        \fn milxQtMain::activeWebView()
        \brief Returns the QWebView object, returns 0 if active window is not a QWebView object
    */
    QWebEngineView* activeWebView();
    /*!
        \fn milxQtMain::setActiveWindow(QWidget *window)
        \brief Makes the window the active window.
    */
    void setActiveWindow(QWidget *currentWindow);

    /**
        \fn milxQtMain::newTab()
        \brief Creates a new tab for the workspace.
    */
    void newTab();
    /*!
        \fn milxQtMain::open()
        \brief Opens data for viewing in the current tab.
    */
    bool open();
    /*!
        \fn milxQtMain::openCollection()
        \brief Opens a collection of models for viewing in the same window.
    */
    void openCollection();
    /*!
        \fn milxQtMain::openSeries()
        \brief Opens a DICOM series for viewing in the same window.
    */
    void openSeries();
    /**
        \fn milxQtMain::openRecentFile()
        \brief Opens recent file (from recent file menu) for analysis and display.
    */
    void openRecentFile();
    /*!
        \fn milxQtMain::loadFiles(const QStringList &filenames)
        \brief Opens multiple files for viewing in the current tab.
    */
    void loadFiles(const QStringList &filenames);
    /*!
        \fn milxQtMain::loadFile(const QString &filename)
        \brief Opens a file for viewing in the current tab.
    */
    bool loadFile(const QString &filename);
    /*!
        \fn milxQtMain::setCurrentFile(const QString &fileName)
        \brief Saves the recent opened/saved file into the recent files list.
    */
    void setCurrentFile(const QString &fileName);
    /*!
        \fn milxQtMain::updateRecentFileActions()
        \brief Updates the recent files list in the file menu.
    */
    void updateRecentFileActions();
    /*!
        \fn milxQtMain::updateWindowMenu()
        \brief Updates the windows list in the windows menu.
    */
    QActionGroup* updateWindowMenu();
    /*!
        \fn milxQtMain::updateWindowListMenu(bool applyMapper = true)
        \brief Updates the windows list in the windows list menu. This is a raw list for children to use.
    */
    QActionGroup* updateWindowListMenu(bool applyMapper = true);
    /*!
        \fn milxQtMain::updateImportFromMenu(bool applyMapper = false)
        \brief Updates the Import form menu in the windows list menu. This is a raw list for children to use.
    */
    QActionGroup* updateImportFromMenu(bool applyMapper = false);
    /*!
        \fn milxQtMain::updateWindowsWithValue(int value)
        \brief Updates the windows based on single value (in percent), such as from a slider etc.
    */
    void updateWindowsWithValue(int value);
    /*!
    \fn milxQtMain::updateWindowsWithAutoLevel()
    \brief Updates the windows with auto levelling.
    */
    void updateWindowsWithAutoLevel();
    /*!
    \fn milxQtMain::updateWindowsWithRefresh()
    \brief Updates the windows with default views and window levels.
    */
    void updateWindowsWithRefresh();
    /*!
    \fn milxQtMain::updateWindowsWithCursors()
    \brief Updates the windows with cursors.
    */
    void updateWindowsWithCursors();
    /*!
        \fn milxQtMain::updateWindowsWithView(int value)
        \brief Updates the windows view, such as from a combo box etc.
    */
    void updateWindowsWithView(int value);
    /*!
        \fn milxQtMain::updateWindowsWithViewType(int value)
        \brief Updates the windows view type based on type set, such as from a combo box etc.
    */
    void updateWindowsWithViewType(int value);
    /*!
        \fn milxQtMain::updateWindowsWithViewOrientation(int value)
        \brief Updates the windows view orientation based on orientation set, such as from a combo box etc.
    */
    void updateWindowsWithViewOrientation(int value);
    /*!
        \fn milxQtMain::windowActionList(QMenu *menuForList, bool groupTogether = true, bool applyMapper = false)
        \brief Creates an action list of windows currently opened to the menu provided.

        Action Group is applied only if Boolean is set to true.
    */
    QActionGroup* windowActionList(QMenu *menuForList, bool groupTogether = true, bool applyMapper = false);

    /*!
        \fn milxQtMain::save(QString filename = "")
        \brief Saves the data in the current tab to file. The saving is content sensitive and depends on the extension provided.
    */
    void save(QString filename = "");
    /*!
        \fn milxQtMain::saveScreen(QString filename = "")
        \brief Saves a screenshot of the active window in the current tab to file. The saving is content sensitive and depends on the extension provided.
    */
    void saveScreen(QString filename = "");

    /*!
        \fn milxQtMain::close()
        \brief Closes window but ensures all windows and tabs are porperly disposed of
    */
    void close();

    /*!
        \fn milxQtMain::setTabName(QMdiSubWindow *window)
        \brief Set the tab name.
    */
    void setTabName(QMdiSubWindow *fromWindow);
    /*!
        \fn milxQtMain::setTabName(QWidget *window)
        \brief Set the tab name.
    */
    void setTabName(QWidget *fromWindow);
    /*!
        \fn milxQtMain::setTabName(const QString newName)
        \brief Set the tab name of the current tab.
    */
    void setTabName(const QString newName);
    /*!
        \fn milxQtMain::closeTab(int index)
        \brief Close the tab.
    */
    void closeTab(int index);
    /*!
        \fn milxQtMain::closeTabActiveWindow()
        \brief Close active window in the tab.
    */
    inline void closeTabActiveWindow()
    {
        qobject_cast<QMdiArea *>(workspaces->currentWidget())->closeActiveSubWindow();
    }
    /*!
        \fn milxQtMain::closeTabAllWindows()
        \brief Close all the windows in the current tab.
    */
    inline void closeTabAllWindows()
    {
        qobject_cast<QMdiArea *>(workspaces->currentWidget())->closeAllSubWindows();
    }
    /*!
        \fn milxQtMain::cascadeTab()
        \brief Cascade all the windows in the current tab.
    */
    inline void cascadeTab()
    {
        qobject_cast<QMdiArea *>(workspaces->currentWidget())->cascadeSubWindows();
    }
    /*!
        \fn milxQtMain::tileTab()
        \brief Tile all the windows in the current tab.
    */
    inline void tileTab()
    {
        qobject_cast<QMdiArea *>(workspaces->currentWidget())->tileSubWindows();
    }
    /*!
        \fn milxQtMain::tileTabVertically()
        \brief Tile all the windows vertically in the current tab.
    */
    void tileTabVertically();
    /*!
        \fn milxQtMain::tileTabHorizontally()
        \brief Tile all the windows horizontally in the current tab.
    */
    void tileTabHorizontally();
    /*!
        \fn milxQtMain::helpContents()
        \brief Show the help contents browser.
    */
    void helpContents();
	/*!
	\fn milxQtMain::about()
	\brief Show the about information.
	*/
	void about();
	/*!
	\fn milxQtMain::controls()
	\brief Show the controls information.
	*/
	void controls();
	/*!
        \fn milxQtMain::preferences()
        \brief Show the preferences for the program.
    */
    void preferences();
    /*!
        \fn milxQtMain::working(int value)
        \brief Show progress bar when computing work. NEgative value shows a busy bar.
    */
    void working(int value);
    /*!
        \fn milxQtMain::done(int value)
        \brief Hide progress bar when done computing.
    */
    void done(int value);

    /*!
        \fn milxQtMain::display(milxQtRenderWindow* rnd)
        \brief Handles the displaying of renders as they are produced.
    */
    void display(milxQtRenderWindow* rnd);
    /*!
        \fn milxQtMain::predisplay(milxQtImage* img)
        \brief Handles the pre-displaying tasks for images, such as checking view type and if 3D viewing is needed etc.

        This has to be used in a number of places so its made modular.
    */
    void predisplay(milxQtImage* img);
    /*!
        \fn milxQtMain::display(milxQtImage* img)
        \brief Handles the displaying of images as they are produced.
    */
    void display(milxQtImage* img);
    /*!
        \fn milxQtMain::display(vtkImageData* img, QString nameOfImage)
        \brief Handles the displaying of images as they are produced.
    */
    void display(vtkImageData* img, QString nameOfImage);
    /*!
        \fn milxQtMain::display(milxQtModel* mdl)
        \brief Handles the displaying of models as they are produced.
    */
    void display(milxQtModel* mdl);
    /*!
        \fn milxQtMain::display(vtkPolyData* newModel, QString nameOfModel)
        \brief Handles the displaying of polydata as they are produced.
    */
    void display(vtkPolyData* newModel, QString nameOfModel);
    /*!
        \fn milxQtMain::display(milxQtUnifiedWindow* uni)
        \brief Handles the displaying of unified (multi) display windows as they are produced.
    */
    void display(milxQtUnifiedWindow* uni);
    /*!
        \fn milxQtMain::display(vtkPolyDataCollection *modelCollection, QStringList &filenames)
        \brief Handles the display of collections
    */
    void display(vtkPolyDataCollection *modelCollection, QStringList &filenames);

    /**
        \fn milxQtMain::addRender(milxQtRenderWindow *rnd)
        \brief Adds a generic render widget and links its results.
    */
    void addRender(milxQtRenderWindow *rnd);
    /**
        \fn milxQtMain::addImage(milxQtImage *img)
        \brief Adds an image widget and links its results.
    */
    void addImage(milxQtImage *img);
    /**
        \fn milxQtMain::addModel(milxQtModel *mdl)
        \brief Adds a model widget and links its results.
    */
    void addModel(milxQtModel *mdl);
    /**
        \fn milxQtMain::addUnifiedWindow(milxQtUnifiedWindow *uni)
        \brief Adds a unified window and links its results.
    */
    void addUnifiedWindow(milxQtUnifiedWindow *uni);
    /**
        \fn milxQtMain::addToUnifiedWindow(milxQtImage *img)
        \brief Adds an image to the unified window and links its results.
    */
    inline void addToUnifiedWindow(milxQtImage *img)
    {
        currentUnifiedWindow->addToWindow(img);
    }
    /**
        \fn milxQtMain::addToUnifiedWindow(milxQtModel *mdl)
        \brief Adds a model to the unified window and links its results.
    */
    inline void addToUnifiedWindow(milxQtModel *mdl)
    {
        currentUnifiedWindow->addToWindow(mdl);
    }
    /**
        \fn milxQtMain::showUnifiedWindow()
        \brief Shows the unified window.
    */
    inline void showUnifiedWindow()
    {
        display(currentUnifiedWindow);
    }
    /**
        \fn milxQtMain::getUnifiedWindow()
        \brief Returns the unified window.
    */
    inline milxQtUnifiedWindow* getUnifiedWindow()
    {
        return currentUnifiedWindow;
    }
    /**
        \fn milxQtMain::cleanUpOnClose(QWidget *win)
        \brief Removes references to windows that are closed so they are no longer "open"
    */
    void cleanUpOnClose(QWidget *win);

    //Window traversal
    /**
        \fn milxQtMain::getListOfWindows()
        \brief Get a list of widgets/windows that are in the current tab
    */
    inline QList<QMdiSubWindow*> getListOfWindows()
    {
        return qobject_cast<QMdiArea *>(workspaces->currentWidget())->subWindowList();
    }
    /**
        \fn milxQtMain::getNumberOfWindows()
        \brief Return the number of windows in the current tab.
    */
    inline int getNumberOfWindows()
    {
        return qobject_cast<QMdiArea *>(workspaces->currentWidget())->subWindowList().size();
    }
    /**
        \fn milxQtMain::getNumberOfImageWindows()
        \brief Return the number of image windows in the current tab.
    */
    int getNumberOfImageWindows();
    /**
        \fn milxQtMain::getNumberOfModelWindows()
        \brief Return the number of model windows in the current tab.
    */
    int getNumberOfModelWindows();
    /**
        \fn milxQtMain::getNumberOfTabs()
        \brief Return the number of tabs in the main window.
    */
    inline int getNumberOfTabs()
    {
        return workspaces->count();
    }
    /**
        \fn milxQtMain::initialiseWindowTraversal()
        \brief Initialise the window iterator to the first window opened.
    */
    inline void initialiseWindowTraversal()
    {
        windowIterator = 0;
    }
    /**
        \fn milxQtMain::currentWindow()
        \brief Get the current window opened in the current tab.

        Use initialiseWindowTraversal() is reset to the first window.
    */
    milxQtWindow* currentWindow();
    /**
        \fn milxQtMain::nextWindow()
        \brief Get the next window opened in the current tab.

        Use initialiseWindowTraversal() is reset to the first window.
    */
    milxQtWindow* nextWindow();
    /**
        \fn milxQtMain::nextRenderWindow()
        \brief Get the next render window opened in the current tab

        Use initialiseWindowTraversal() is reset to the first window.
    */
    milxQtRenderWindow* nextRenderWindow();
    /**
        \fn milxQtMain::nextModel()
        \brief Get the next model opened in the current tab

        Use initialiseWindowTraversal() is reset to the first window.
    */
    milxQtModel* nextModel();
    /**
        \fn milxQtMain::nextImage()
        \brief Get the next image opened in the current tab

        Use initialiseWindowTraversal() is reset to the first window.
    */
    milxQtImage* nextImage();

    //Image-Model Inter-Members
    /*!
        \fn milxQtMain::imageToSurface(vtkSmartPointer<vtkImageData> img, const float value = numeric_limits<float>::max())
        \brief Compute image to surface process.
    */
    void imageToSurface(vtkSmartPointer<vtkImageData> img, const float value = numeric_limits<float>::max());
    /*!
        \fn milxQtMain::imageToPolyData(vtkSmartPointer<vtkImageData> img)
        \brief Compute image to poly data process.
    */
    void imageToPolyData(vtkSmartPointer<vtkImageData> img);
    /*!
        \fn milxQtMain::imageToPseudoImage(vectorImageType::Pointer img)
        \brief Compute pseudo image from vector image.
    */
    void imageToPseudoImage(vectorImageType::Pointer img);
    /*!
        \fn milxQtMain::imageToVectorField(vectorImageType::Pointer img, floatImageType::Pointer magImg, int subsampleFactor = 0, float scaling = 0.0)
        \brief Compute vector field from image.
    */
    void imageToVectorField(vectorImageType::Pointer img, floatImageType::Pointer magImg, int subsampleFactor = 0, float scaling = 0.0);
    /*!
        \fn milxQtMain::imageToTensorField(vectorImageType::Pointer img, floatImageType::Pointer magImg, int subsampleFactor = 0, float scaling = 0.0)
        \brief Compute tensor field from image.
    */
    void imageToTensorField(vectorImageType::Pointer img, floatImageType::Pointer magImg, int subsampleFactor = 0, float scaling = 0.0);
    /*!
        \fn milxQtMain::imageToStreamLines(vectorImageType::Pointer img, floatImageType::Pointer magImg, size_t subsampleFactor = 0)
        \brief Compute streamlines from vector image from the extents supplied (pixels of which are used as start points for the streams).
    */
    void imageToStreamLines(vectorImageType::Pointer img, floatImageType::Pointer magImg, size_t subsampleFactor = 0);
    /*!
        \fn milxQtMain::imageToVolume(vtkSmartPointer<vtkImageData> img, bool eightbit)
        \brief Display image as volume rendering
    */
    void imageToVolume(vtkSmartPointer<vtkImageData> img, bool eightbit);
    /*!
        \fn milxQtMain::imageToPlot(vtkSmartPointer<vtkImageData> img, int displaceAxis = 2)
        \brief Display image as a surface plot.

        displaceAxis is the axis to displace the surface along.
    */
    void imageToPlot(vtkSmartPointer<vtkImageData> img, int displaceAxis = 2);
    /*!
        \fn milxQtMain::tableToPlot(vtkSmartPointer<vtkTable> tbl, QString title)
        \brief Display table as a plot
    */
    void tableToPlot(vtkSmartPointer<vtkTable> tbl, QString title);
    /*!
        \fn milxQtMain::voxeliseSurface(vtkSmartPointer<vtkPolyData> surface)
        \brief Convert poly data to an image volume.
    */
    void voxeliseSurface(vtkSmartPointer<vtkPolyData> surface);

    //Batch
    /*!
        \fn milxQtMain::imagesMix()
        \brief Blend all open images together to create new image usaing a mixer dialog.

        All data from open image windows are blended to the first image.
    */
    void imagesMix();
    /*!
        \fn milxQtMain::imagesBlend(QVector<float> opacities)
        \brief Blend all open images together to create new image. 'opacities' would hold a list of floats for blending each image

        All data from open image windows are blended to the first image.
    */
    void imagesBlend(QVector<float> opacities);
    /*!
        \fn milxQtMain::imagesAdd()
        \brief Add all open images together to create new image.

        All data from open image windows are added to the first image.
    */
    void imagesAdd();
    /*!
        \fn milxQtMain::imagesAverage()
        \brief Average all open images together to create new image.
    */
    void imagesAverage();
    /*!
        \fn milxQtMain::imagesSubtract()
        \brief Subtract all open images together to create new image.

        All data from open image windows are subtracted from the first image.
    */
    void imagesSubtract();
    /*!
    \fn milxQtMain::imagesMultiply()
    \brief Multiply all open images together to create new image.

    All data from open image windows are multiplied to the first image.
    */
    void imagesMultiply();
    /*!
        \fn milxQtMain::imagesConvolve()
        \brief Convolve all open images together to create new image.

        All data from open image windows are convolved into the first image.
    */
    void imagesConvolve();
    /*!
        \fn milxQtMain::imagesMergeLabels()
        \brief Merge all open labelled images together to create new labelled image.

        All data from open image windows are merged into a single image.
    */
    void imagesMergeLabels();

    /**
        \fn milxQtMain::unify()
        \brief Add the active window to the current multi (unified) display.
    */
    void unify();
    /**
        \fn milxQtMain::link()
        \brief Link all windows in the current tab. Changing one window camera updates all others.
    */
    inline void link()
    {
        actionLinkWindows->setChecked(true);
    }
    /**
        \fn milxQtMain::unlink()
        \brief Unlink all windows in the current tab.
    */
    inline void unlink()
    {
        actionLinkWindows->setChecked(false);
    }
    /**
        \fn milxQtMain::update()
        \brief Update the GUI elements, such as menus etc. to most up-to-date status.
    */
    void update();
    /*!
        \fn milxQtMain::writeSettings()
        \brief Write the necessary GUI settings/state for the main class. If the settings have been reset, no settings are written until re-read.
    */
    void writeSettings();
	/*!
	\fn milxQtMain::resetSettings()
	\brief Reset the necessary GUI settings/state for the main class.
	*/
	void resetSettings();
    /*!
        \fn milxQtMain::readSettings()
        \brief Load the necessary GUI settings/state for the main class.
    */
    void readSettings();

    /*!
        \fn milxQtMain::loadPlugins()
        \brief Loads plugins (which are assumed to be DLLs) according to the milxQtPluginInterface.
    */
    bool loadPlugins();
    /*!
        \fn milxQtMain::getPlugins()
        \brief Returns a list of the loaded plugins (which are assumed to be DLLs) according to the milxQtPluginInterface.
    */
    inline QList< QPointer<milxQtPluginInterface> > getPlugins()
    {
        return plugins;
    }

    //Print Members
    /*!
        \fn milxQtMain::printError(QString msg)
        \brief Error message wrapper for console.
    */
    inline void printError(QString msg)
    {
        console->printError(msg);
    }
//    {   cerr << msg.toStdString() << endl;   }
    /*!
        \fn milxQtMain::printWarning(QString msg)
        \brief Warning message wrapper for console.
    */
    inline void printWarning(QString msg)
    {
        console->printWarning(msg);
    }
//    {   cerr << msg.toStdString() << endl;   }
    /*!
        \fn milxQtMain::printDebug(QString msg)
        \brief Debug message wrapper for console.
    */
    inline void printDebug(QString msg)
    {
        console->printDebug(msg);
    }
//    {   cerr << msg.toStdString() << endl;   }
    /*!
        \fn milxQtMain::printInfo(QString msg)
        \brief Info message wrapper for console.
    */
    inline void printInfo(QString msg)
    {
        console->printInfo(msg);
    }
//    {   cerr << msg.toStdString() << endl;   }

protected slots:
    /*!
        \fn milxQtMain::redirectWindowActivated(QMdiSubWindow *win)
        \brief Redirect the workspace signal to milxQtMain object level.
    */
    inline void redirectWindowActivated(QMdiSubWindow *win)
    {
		if (win) {
			emit windowActivated(win->widget());
		}
    }
    /*!
        \fn milxQtMain::transferViewToWindows(vtkObject *obj, unsigned long, void *client_data, void *, vtkCommand *command)
        \brief If the windows are linked, copy viewing options from one window to all others.
    */
    void transferViewToWindows(vtkObject *obj, unsigned long, void *client_data, void *, vtkCommand *command);
    /*!
        \fn milxQtMain::dataMenu()
        \brief Create the data menu based on Window selected.
    */
    void dataMenu();
    /*!
        \brief Update the Qt events, used to keep UI responsive
    */
    inline void updateQtEvents()
    {
        qApp->processEvents();
    }

signals:
    /**
        \fn milxQtMain::windowActivated(QWidget * win)
        \brief Signal for when a window has been activated.
    */
    void windowActivated(QWidget * win);
    /**
        \fn milxQtMain::displayed(milxQtRenderWindow*)
        \brief Signal for when a render window has been displayed.
    */
    void displayed(milxQtRenderWindow*);
    /**
        \fn milxQtMain::displayed(milxQtImage*)
        \brief Signal for when an image has been displayed.
    */
    void displayed(milxQtImage*);
    /**
        \fn milxQtMain::displayed(milxQtModel*)
        \brief Signal for when a model has been displayed.
    */
    void displayed(milxQtModel*);
    /**
        \fn milxQtMain::displayed(milxQtUnifiedWindow*)
        \brief Signal for when a unified window has been displayed.
    */
    void displayed(milxQtUnifiedWindow*);
    /**
        \fn milxQtMain::updatedWindowListMenu(QActionGroup *)
        \brief Signal for when the window list menu is updated
    */
    void updatedWindowListMenu(QActionGroup *);
    /**
        \fn milxQtMain::updatedWindowListMenu(QMenu *)
        \brief Signal for when the window list menu is updated
    */
    void updatedWindowListMenu(QMenu *);
    /**
        \fn milxQtMain::updatedImportFromMenu(QActionGroup *)
        \brief Signal for when the import from menu is updated
    */
    void updatedImportFromMenu(QActionGroup *);
    /**
        \fn milxQtMain::updatedImportFromMenu(QMenu *)
        \brief Signal for when the import from menu is updated
    */
    void updatedImportFromMenu(QMenu *);

protected:
    bool debugMode; //!< Debug output mode?

    //Settings
    bool whiteBackground; //!< Prefer white backgrounds?
    bool humanGlyph; //!< Prefer showing human glyph?
    int subWindowSize; //!< Window size of child windows to use
    int maxProcessors; //!< Max processors to use
    int magnifyFactor; //!< Screenshot magnify factor
    bool timestamping; //!< Prefer showing timestamp?
	QString appTheme; //!< The application style theme
    bool interpolationImages; //!< Prefer showing interpolation for images?
    bool orientationImages; //!< Prefer applying orientation to images?
    bool interpolationModels; //!< Prefer applying interpolation to models?
    bool scalarBarModels; //!< Show scalar bar for models?
	bool resettingInterface; //!< Flag if reseting interface, ignore write settings on close
	bool enableColourMaps = true; //!< Enable custom colour maps

    enum { MaxRecentFiles = 10 };
    //Menus (hierarchical deletion)
    QMenuBar* menuBar; //!< Menu bar for the window
    QMenu* menuFile; //!< File menu
    QMenu* menuImages; //!< File menu
    QMenu* menuData; //!< Data menu
    //QMenu* menuNew; //!< File menu
    QMenu* menuWindows; //!< Windows menu
    QMenu* menuHelp; //!< Help menu

    size_t progressCallCount; //!< Takes account of how many calls of working have been made.
    int windowIterator; //!< Keep track of current window being traversed

    //File dialog attributes
    QString saveSupport; //!< Save file extension support list, cats all known extensions.
    QString openSupport; //!< Load file extension support list, cats all known extensions.

    //----File---- (hierarchical deletion)
    QAction* actionOpen; //!< open action
    QAction* actionOpenSeries; //!< open series action
    QAction* actionOpenCollect; //!< open collection action
    QAction* actionSave; //!< save action
    QAction* actionSaveScreen; //!< save action
    QAction* actionCloseActive; //!< Close active window in tab action
    QAction* actionCloseAll; //!< Close all windows in tab action
    QList<QAction*> actionsRecentFile; //!< Array of recent file actions
    QAction* actionRecentFileSeparator; //!< Pointer to separator so it can be turned on and off
    QAction* actionExit; //!< Exit action
    //----New---- (hierarchical deletion)
    QAction* actionNewTab; //!< New Tab action
    //----Images---- (hierarchical deletion)
    QAction* actionBlendImages; //!< blend images batch operation action
    QAction* actionAddImages; //!< add images batch operation action
    QAction* actionAverageImages; //!< average images batch operation action
    QAction* actionSubtractImages; //!< subtract images batch operation action
    QAction* actionMultiplyImages; //!< multiply images batch operation action
    QAction* actionConvolveImages; //!< Convolve images batch operation action
    QAction* actionMergeLabels; //!< Merge labelled images batch operation action
    //----Windows---- (hierarchical deletion)
    QAction* actionLinkWindows; //!< Using linked cameras for all windows in a tab?
    QAction* actionCascade; //!< Cascade windows in workspace action
    QAction* actionTile; //!< Tile windows in workspace action
    QAction* actionTileVertically; //!< Tile vertically windows in workspace action
    QAction* actionTileHorizontally; //!< Tile horizontally windows in workspace action
    QAction* actionConsole; //!< toggle action for console
    //----Help---- (hierarchical deletion)
    QAction* actionContents; //!< Action for showing contents of help
    QAction* actionPreferences; //!< Action for showing program preferences
    QAction* actionControls; //!< Action for showing controls information
    QAction* actionAbout; //!< Action for showing about information

    //Common actions/menus
    QAction* actionCompare; //!< Action for comparing data
    QMenu* menuWindowList; //!< Menu for list of windows
    QMenu* importFromMenu; //!< Menu import from data menu

    //Toolbars
    QToolBar* fileToolBar; //!< Some actions from file menu
//    QToolBar* editToolBar; //!< Some actions from edit menu, like copy, paste etc.
    QToolBar* windowToolBar; //!< Some actions from window menu, like tile, cascade etc.
    QToolBar* defaultToolBar; //!< Some actions from default stuff menu, like view etc.
    QToolBar* imageToolBar; //!< Some actions for images, like view etc.

    //Image toolbar actions
    QAction* actionImageText; //!< toggle text annotate mode
    QSlider* imageLevelSlider; //!< adjust window level of image display
    QPushButton* imageLevelButton; //!< auto level button
    //QDial* imageLevelDial; //!< adjust window level of image display via a dial
    QPushButton* refreshButton; //!< window refresh button
    QPushButton* cursorButton; //!< crosshairs button

    //Workspaces (hierarchical deletion)
    QTabWidget* workspaces; //!< Pointer to the Workspace environment for the user.
    QSignalMapper *windowMapper; //!< Mapper of mulit-connections
    //Plugins (Smart Pointer deletion)
    QList< QPointer<milxQtPluginInterface> > plugins; //!< List of plugins loaded succesfully.
    QList< QAction* > renderExtsActions; //!< List of render window extenstion actions loaded succesfully.
    QList< QAction* > modelExtsActions; //!< List of model extenstion actions loaded succesfully.
    QList< QAction* > imageExtsActions; //!< List of image extenstion actions loaded succesfully.
    QList< QAction* > dockActions; //!< List of dock actions of dock widgets loaded succesfully.
    //Bar (hierarchical deletion)
    QProgressBar *progressBar; //!< Progress bar for computation
    //Console
    QPointer<milxQtConsole> console; //!< console docked window

    //Default Comboboxes
    QComboBox *defaultViewBox; //!< Box for default view
    QComboBox *defaultViewTypeBox; //!< Box for default view type (1 view or scanner type/4 view)
    QComboBox *defaultOrientationTypeBox; //!< Box for default orientation type (radiological or neurological orientation)

    QPointer<milxQtUnifiedWindow> currentUnifiedWindow; //!< Current window for multi-display
    QList< QPointer<milxQtModel> > modelWindows; //!< List of model windows opened.
    QList< QPointer<milxQtImage> > imageWindows; //!< List of model windows opened.

    vtkSmartPointer<vtkEventQtSlotConnect> Connector; //!< VTK Events to slots convertor

    /*!
        \fn milxQtMain::linkProgressEventOf(vtkObject * obj)
        \brief Link the progress of filters etc to keep the UI responsive.
    */
    void linkProgressEventOf(vtkObject * obj);
    /*!
        \fn milxQtMain::commonChildProperties(QWidget *widget)
        \brief Sets all the common child properties of the data windows such as images, plots etc.
    */
    void commonChildProperties(QWidget *widget);
    /*!
        \fn milxQtMain::createMenu()
        \brief Creates the menu actions.
    */
    void createMenu();
    /*!
        \fn milxQtMain::createComboBoxes()
        \brief Creates the combo boxes for main window.
    */
    void createComboBoxes();
    /*!
        \fn milxQtMain::createToolBars()
        \brief Creates the toolbar for the main window.
    */
    void createToolBars();
    /*!
        \fn milxQtMain::createConnections()
        \brief Creates the signals and slots connections within the main window.
    */
    void createConnections();
    /*!
        \fn milxQtMain::createProgressBar()
        \brief Creates the progress bar in the status bar within the main window.
    */
    void createProgressBar();
    /*!
        \fn milxQtMain::setupTooltips()
        \brief Assign tooltips to each action for users benefit. The tooltips explains the function of each action
    */
    void setupTooltips();
    /*!
    	\fn milxQtMain::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);
    /*!
        \fn milxQtMain::strippedName(const QString &fullFileName)
        \brief Returns the filename with the path stripped.
    */
    inline QString strippedName(const QString &fullFileName)
    {
        return QFileInfo(fullFileName).fileName();
    }
    /*!
        \fn milxQtMain::dragEnterEvent(QDragEnterEvent *event)
        \brief Part of the Drag and Drop feature members. Tells what drags to accept.
    */
    void dragEnterEvent(QDragEnterEvent *event);
    /*!
        \fn milxQtMain::dropEvent(QDropEvent *event)
        \brief Part of the Drag and Drop feature members. Opens the dropped files.
    */
    void dropEvent(QDropEvent *event);
    /*!
        \fn milxQtMain::closeEvent(QCloseEvent *event)
        \brief Tasks to complete when closing.
    */
    void closeEvent(QCloseEvent *event);

private:

};

#endif // MILXQTMAIN_H
