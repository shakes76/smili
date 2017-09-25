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
#include "milxQtMain.h"
//Qt

//VTK
#include <vtkWindowToImageFilter.h>
#include <vtkRendererCollection.h>
#include <vtkCamera.h>
#include <vtkLookupTable.h>
#include <vtkImageBlend.h>
#include <vtkTensor.h>
#include <vtkMultiThreader.h>
#include <vtkMath.h>
//milxQt
#include "milxQtFile.h"
#include "milxQtPlot.h"

//Forms
#include "milxQtAboutForm.h"
#include "milxQtPreferencesForm.h"

milxQtMain::milxQtMain(QWidget *theParent) : QMainWindow(theParent)
{
    debugMode = true;

    //Settings
    whiteBackground = false;
    humanGlyph = true;
    subWindowSize = minWindowSize*2;
    maxProcessors = milx::NumberOfProcessors();
    if(maxProcessors > 1)
      maxProcessors = milx::NumberOfProcessors()/2;
    magnifyFactor = 1;
    timestamping = true;
    interpolationImages = true;
    orientationImages = true;
    interpolationModels = false;
    scalarBarModels = false;
	resettingInterface = false;

    progressCallCount = 0;
    windowIterator = 0;

    saveSupport = extensionsSave.c_str();
    openSupport = extensionsOpen.c_str();

    currentUnifiedWindow = new milxQtUnifiedWindow;
    currentUnifiedWindow->setName("Multi-Display");
    currentUnifiedWindow->setDeletableOnClose(false);

    Connector = vtkSmartPointer<vtkEventQtSlotConnect>::New();

    ///Connect Progress Object for progress updates in file and imaging
    linkProgressEventOf(milx::ProgressUpdates->GetUpdateObject()); //link itk observer propagator
    linkProgressEventOf(milx::VTKProgressUpdates->GetUpdateObject()); //link itk observer propagator

    ///Establish a workspace for user
    workspaces = new QTabWidget;
    workspaces->setTabsClosable(true);
    newTab();
    QMainWindow::setCentralWidget(workspaces); //Hierachy deletion
    windowMapper = new QSignalMapper(this);

    ///Setup Console
    console = new milxQtConsole;
    actionConsole = console->dockWidget()->toggleViewAction();
    actionConsole->setIcon(QIcon(":/resources/toolbar/console.png"));
    dockActions.append(actionConsole);
    QObject::connect(console->dockWidget(), SIGNAL(dockLocationChanged(Qt::DockWidgetArea)), console, SLOT(setDockDefaultArea(Qt::DockWidgetArea)));
    addDockWidget(console->dockDefaultArea(), console->dockWidget());
    console->show();

    ///Setup the UI Menus and boxes
    createMenu();
    createComboBoxes();

    ///Setup the UI Toolbar
    createToolBars();

    ///Setup Connections
    createConnections();

    ///Setup tooltips
    setupTooltips();

    ///Load plugins
    loadPlugins();

    ///Setup progress bar
    createProgressBar();

    setAcceptDrops(true);

    setToolTip("<p style='white-space:pre'>Welcome to SMILX. Use the <b>context menu</b> (right click) for all operations.\n Load data from the <b>File Menu</b>.</p>");

    ///Read the application settings/preferences
    readSettings();

    ///Setup ITK Threads
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(maxProcessors);
    vtkMultiThreader::SetGlobalDefaultNumberOfThreads(maxProcessors);

    ///Program Info
    printInfo("--------------------------------------------------------");
    printInfo("sMILX Visualisation Tool for Medical Imaging");
    printInfo("(c) Copyright CSIRO, 2015.");
    printInfo("University of Queensland, Australia.");
    printInfo("Australian e-Health Research Centre, CSIRO.");
    printInfo("SMILI Version: " + QString::number(milx::Version));
    printInfo("Application Version: " + QString::number(milxQtVersion));
    printInfo("Processors to be used: " + QString::number(maxProcessors));
    printInfo("--------------------------------------------------------\n");

    ///Style

    update();
}

milxQtMain::~milxQtMain()
{
    foreach (QPointer<milxQtPluginInterface> loadedPlugin, plugins)
    {
        if(loadedPlugin->isThreaded() && loadedPlugin->isRunning())
            loadedPlugin->quit();
    }
    plugins.clear();
    renderExtsActions.clear();
    modelExtsActions.clear();
    imageExtsActions.clear();
    dockActions.clear();

    modelWindows.clear();
    imageWindows.clear();
    printDebug("Main Destroyed.");
}

QString milxQtMain::activeName()
{
    if(qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow() == 0)
    {
        ///\todo Possible Bug: Handling when no window active.
        return 0;
    }
    else
    {
        if(isActiveImage())
        {
            milxQtImage *childImg = activeImage(); ///Check for image validity already done so safely do...
            return childImg->strippedName();
        }
        else if(isActiveModel())
        {
            milxQtModel *childMdl = activeModel(); ///Check for image validity already done so safely do...
            return childMdl->strippedName();
        }
        else if(isActiveRender())
        {
            milxQtRenderWindow *childRnd = activeRender(); ///Check for image validity already done so safely do...
            return childRnd->strippedName();
        }
        else
            return "Unrecognised Display";
    }
}

QString milxQtMain::activeNamePrefix()
{
    if(qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow() == 0)
    {
        ///\todo Possible Bug: Handling when no window active.
        return 0;
    }
    else
    {
        if(isActiveImage())
        {
            milxQtImage *childImg = activeImage(); ///Check for image validity already done so safely do...
            return childImg->strippedNamePrefix();
        }
        else if(isActiveModel())
        {
            milxQtModel *childMdl = activeModel(); ///Check for image validity already done so safely do...
            return childMdl->strippedNamePrefix();
        }
        else if(isActiveRender())
        {
            milxQtRenderWindow *childRnd = activeRender(); ///Check for image validity already done so safely do...
            return childRnd->strippedNamePrefix();
        }
        else if(isActiveWebView())
        {
            QWebView *webViewer = activeWebView(); ///Check for image validity already done so safely do...
            return webViewer->windowTitle();
        }
        else
            return "Unrecognised Display";
    }
}

void milxQtMain::newTab()
{
    QWorkspace *tmpPtr = new QWorkspace; ///\todo Workspace class is deprecated. Update to MdiArea.
    workspaces->addTab(tmpPtr, tr("Empty"));
    workspaces->setCurrentWidget(tmpPtr); //Hierachy deletion
    tmpPtr->setAttribute(Qt::WA_DeleteOnClose);

    connect(tmpPtr, SIGNAL(windowActivated(QWidget *)), this, SLOT(setTabName(QWidget *)));
    connect(tmpPtr, SIGNAL(windowActivated(QWidget *)), this, SLOT(redirectWindowActivated(QWidget *)));
}

void milxQtMain::addRender(milxQtRenderWindow *rnd)
{
    qobject_cast<QWorkspace *>(workspaces->currentWidget())->addWindow(rnd);

    commonChildProperties(rnd);

    rnd->enableUpdates(statusBar());

    connect(rnd, SIGNAL(nameChanged(const QString &)), this, SLOT(setTabName(const QString &)));
    connect(rnd, SIGNAL(working(int)), this, SLOT(working(int)));
    connect(rnd, SIGNAL(done(int)), this, SLOT(done(int)));
}

void milxQtMain::addImage(milxQtImage *img)
{
    addRender(img);
}

void milxQtMain::addModel(milxQtModel *mdl)
{
    mdl->modelInfo();

    addRender(mdl);
}

void milxQtMain::addUnifiedWindow(milxQtUnifiedWindow *uni)
{
    addRender(uni);
}

void milxQtMain::cleanUpOnClose(QWidget *win)
{
    if(isImage(win))
    {
        milxQtImage *img = qobject_cast<milxQtImage *>(win);
        printInfo("Removing Image " + img->getName() + " from memory");
        imageWindows.removeAll(img);
        img->disconnect(); //disconnect all
        if(currentUnifiedWindow)
            currentUnifiedWindow->removeFromWindow(img);
    }
    else if(isModel(win))
    {
        milxQtModel *mdl = qobject_cast<milxQtModel *>(win);
        printInfo("Removing Model " + mdl->getName() + " from memory");
        modelWindows.removeAll(mdl);
        mdl->disconnect(); //disconnect all
        if(currentUnifiedWindow)
            currentUnifiedWindow->removeFromWindow(mdl);
    }
}

//Slots
bool milxQtMain::isRender(QWidget *win)
{
    if(qobject_cast<milxQtRenderWindow *>(win) == 0)
        return false;
    else
        return true;
}

bool milxQtMain::isActiveRender()
{
    if(activeRender() == 0)
        return false;
    else
        return true;
}

milxQtRenderWindow* milxQtMain::activeRender()
{
    if(QWidget *activeWin = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow())
        return qobject_cast<milxQtRenderWindow *>(activeWin);
    return 0;
}

bool milxQtMain::isImage(QWidget *win)
{
    if(qobject_cast<milxQtImage *>(win) == 0)
        return false;
    else
        return true;
}

bool milxQtMain::isActiveImage()
{
    if(activeImage() == 0)
        return false;
    else
        return true;
}

milxQtImage* milxQtMain::activeImage()
{
    if(QWidget *activeWin = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow())
        return qobject_cast<milxQtImage *>(activeWin);
    return 0;
}

bool milxQtMain::isModel(QWidget *win)
{
    if(qobject_cast<milxQtModel *>(win) == 0)
        return false;
    else
        return true;
}

bool milxQtMain::isActiveModel()
{
    if(activeModel() == 0)
        return false;
    else
        return true;
}

milxQtModel* milxQtMain::activeModel()
{
    if(QWidget *activeWin = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow())
        return qobject_cast<milxQtModel *>(activeWin);
    return 0;
}

bool milxQtMain::isPlot(QWidget *win)
{
    if(qobject_cast<milxQtPlot *>(win) == 0)
        return false;
    else
        return true;
}

bool milxQtMain::isActivePlot()
{
    if(activePlot() == 0)
        return false;
    else
        return true;
}

milxQtModel* milxQtMain::activePlot()
{
    if(QWidget *activeWin = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow())
        return qobject_cast<milxQtPlot *>(activeWin);
    return 0;
}

bool milxQtMain::isUnifiedWindow(QWidget *win)
{
    if(qobject_cast<milxQtUnifiedWindow *>(win) == 0)
        return false;
    else
        return true;
}

bool milxQtMain::isActiveUnifiedWindow()
{
    if(activeUnifiedWindow() == 0)
        return false;
    else
        return true;
}

milxQtUnifiedWindow* milxQtMain::activeUnifiedWindow()
{
    if(QWidget *activeWin = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow())
        return qobject_cast<milxQtUnifiedWindow *>(activeWin);
    return 0;
}

bool milxQtMain::isActiveWebView()
{
    if(activeWebView() == 0)
        return false;
    else
        return true;
}

QWebView* milxQtMain::activeWebView()
{
    if(QWidget *activeWin = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow())
        return qobject_cast<QWebView *>(activeWin);
    return 0;
}

void milxQtMain::setActiveWindow(QWidget *currentWindow)
{
    if(!currentWindow)
        return;
    else
        qobject_cast<QWorkspace *>(workspaces->currentWidget())->setActiveWindow(currentWindow);
}

bool milxQtMain::open()
{
    QFileDialog *fileOpener = new QFileDialog(this);
    QSettings settings("Shekhar Chandra", "milxQt");

    //Add supported file entry
    QString exts = "Supported Files (" + openSupport + ");;";
    QString path = settings.value("recentPath").toString();

    exts += openExts.c_str();

    foreach (QPointer<milxQtPluginInterface> loadedPlugin, plugins)
    {
        if(loadedPlugin->openFileSupport() != "")
        {
            exts += ";;";
            exts += loadedPlugin->openFileSupport();
        }
    }

    QStringList filenames = fileOpener->getOpenFileNames(this,
                            tr("Select File(s) to Open"),
                            path,
                            tr(exts.toStdString().c_str()) ); //!< \todo Check and validate extensions support at Open in Main class

    if(filenames.isEmpty())
        return false;

    loadFiles(filenames);

    return true;
}

void milxQtMain::openCollection()
{
    QPointer<milxQtFile> file = new milxQtFile;
    QStringList filenames;

    /*QMessageBox msgBox;
    msgBox.setText("Open an image (DICOM) series instead?");
    msgBox.setInformativeText("Would you like to open an image (DICOM) series instead of a model collection?.");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();

    if(ret == QMessageBox::Yes)
    {
        printDebug("Opening Image Series.");
        QPointer<milxQtFile> reader = new milxQtFile;
        QPointer<milxQtImage> img = new milxQtImage;  //list deletion

        printDebug("Supported Image formats: " + reader->supportedImageFormats());
        bool success = reader->openImageCollection(img, filenames);

        if(success)
        {
            //name saved internally in openImageCollection
            img->setConsole(console);
            img->generateImage();

            predisplay(img); //do the multi-view stuff if requested
        }
        else
            printError("File format didn't appear to be supported. Check the file or add plugins for support.");
    }
    else
    {*/
        vtkPolyDataCollection* modelCollection = vtkPolyDataCollection::New();

        if(!file->openModelCollection(modelCollection, filenames))
        {
            printError("Unable to read selected files.");
            return;
        }

        display(modelCollection, filenames);
//    }
}

void milxQtMain::openSeries()
{
  QPointer<milxQtFile> reader = new milxQtFile;
  QPointer<milxQtImage> img = new milxQtImage;  //list deletion

  printDebug("Supported Image formats: " + reader->supportedImageFormats());
  bool success = reader->openImageSeries(img);

  if(success)
  {
      //name saved internally in openImageCollection
      img->setConsole(console);
      img->generateImage();

      predisplay(img); //do the multi-view stuff if requested
  }
  else
      printError("An error was encountered in reading the image series");
}

void milxQtMain::openRecentFile()
{
    QAction *action = qobject_cast<QAction *>(sender());
    if(action)
        loadFile(action->data().toString());
}

void milxQtMain::loadFiles(const QStringList &filenames)
{
    for(int j = 0; j < filenames.size(); j ++)
        loadFile(filenames[j]); //!< \todo Add error handling
}

bool milxQtMain::loadFile(const QString &filename)
{
    QPointer<milxQtFile> reader = new milxQtFile;
    bool success = false;

    if (filename.isEmpty())
        return success;

    //Convert string to native paths
    QString nativeFilename = QDir::toNativeSeparators(filename);

    ///Check filename with plugins
    printDebug("Check Plugins if they can open the file.");
    foreach (QPointer<milxQtPluginInterface> loadedPlugin, plugins)
    {
        if(!loadedPlugin->hasOpenSupport())
            continue; //Skip trying to check

        QStringList extensionsList = loadedPlugin->openExtensions();

        foreach (QString extension, extensionsList)
        {
            printDebug("Exts: " + extension);
            if(filename.contains(extension, Qt::CaseInsensitive))
            {
                loadedPlugin->setFileName(filename);
                if(loadedPlugin->isThreaded()) ///Check if threaded
                {
                    printInfo("Using threaded opening for plugin.");
                    loadedPlugin->preStartTasks();
                    loadedPlugin->start();
                    //loadedPlugin->wait();
                    success = true;
                }
                else ///if not then use serial methods
                {
                    printInfo("Using standard opening for plugin.");
                    loadedPlugin->open(nativeFilename);
                }

                QPointer<milxQtRenderWindow> renWin = loadedPlugin->genericResult();
                if(renWin)
                {
                    printInfo("Loading Plugin Generic Result");
                    if(renWin->getName() == "")
                        renWin->setName(filename);
                    renWin->setConsole(console);
                    renWin->generateRender();
                    display(renWin);
                    success = true;
                }

                QPointer<milxQtModel> model = loadedPlugin->modelResult();
                if(model)
                {
                    printInfo("Loading Plugin Model Result");
                    if(model->getName() == "")
                        model->setName(filename);
                    model->setConsole(console);
                    model->generateModel();
                    display(model);
                    success = true;
                }

                QPointer<milxQtImage> image = loadedPlugin->imageResult();
                if(image)
                {
                    printInfo("Loading Plugin Image Result");
                    if(image->getName() == "")
                        image->setName(filename);
                    image->setConsole(console);
                    image->generateImage();
                    display(image);
                    success = true;
                }
            }

            if(success)
                break;
        }

        if(success && !loadedPlugin->isThreaded())//!< If success, Plugin loaded it so return, else continue
        {
            printInfo("Loaded File via Plugin, now Updating");
            loadedPlugin->update();
            setCurrentFile(nativeFilename);
            return success;
        }
    }

    ///Check filename with Model and Image
    printDebug("File format not supported by plugins, load using standard methods.");
    if(filename.contains(".vtp", Qt::CaseInsensitive) || filename.contains(".vtk", Qt::CaseInsensitive)
            || filename.contains(".ply", Qt::CaseInsensitive) || filename.contains(".obj", Qt::CaseInsensitive) || filename.contains(".stl", Qt::CaseInsensitive))
    {
        printDebug("Opening Model.");
        QPointer<milxQtModel> model = new milxQtModel; //list deletion
        success = reader->openModel(nativeFilename, model);

        if(success)
        {
            model->setName(filename);
            model->setConsole(console);
            model->generateModel();
            display(model);
        }
    }
    else if(filename.contains(".csv", Qt::CaseInsensitive)  || filename.contains(".dat", Qt::CaseInsensitive) || filename.contains(".txt", Qt::CaseInsensitive))
    {
        printDebug("Opening Text Delimited File.");
        //Load the vertex table from CSV file
        vtkSmartPointer<vtkTable> table;
        QPointer<milxQtPlot> plot = new milxQtPlot;
        success = milx::File::OpenDelimitedText(nativeFilename.toStdString(), table);

        int surfaceRet = QMessageBox::No, dimensionRet = QMessageBox::No;
        int xCol = 0, yCol = 1, zCol = 2;
        if(table->GetNumberOfColumns() > 2) ///Ambiguous, choosing 2 columns or 3D scatter or surface plot?
        {
            QMessageBox msgBox;
            msgBox.setText("Choose Plot Type");
            msgBox.setInformativeText("Do you want a surface plot instead of a scatter plot?");
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
            msgBox.setDefaultButton(QMessageBox::No);
            surfaceRet = msgBox.exec();

            if(surfaceRet == QMessageBox::No)
            {
                msgBox.setText("Choose Plot Dimension");
                msgBox.setInformativeText("Do you want a 3D scatter plot?");
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                msgBox.setDefaultButton(QMessageBox::No);
                dimensionRet = msgBox.exec();

                bool ok1 = true, ok2 = true, ok3 = true;
                xCol = QInputDialog::getInt(this, tr("Please Provide column number for x axis"),
                                      tr("Column:"), 0, 0, numeric_limits<int>::max(), 1, &ok1);
                yCol = QInputDialog::getInt(this, tr("Please Provide column number for y axis"),
                                      tr("Column:"), 1, 0, numeric_limits<int>::max(), 1, &ok2);

                if(dimensionRet == QMessageBox::Yes) //3D
                {
                    zCol = QInputDialog::getInt(this, tr("Please Provide column number for z axis"),
                                      tr("Column:"), 2, 0, numeric_limits<int>::max(), 1, &ok3);
                }

                if(!ok1 || !ok2 || !ok3)
                    return success;

                if(dimensionRet == QMessageBox::No) //2D
                    plot->setPlotType2D(xCol, yCol);
                else
                    plot->setPlotType3D(xCol, yCol, zCol);
            }
            else
            {
                plot->setPlotTypeSurface();
            }
        }

        //Connect for display
        connect(plot, SIGNAL(resultAvailable(milxQtModel*)), this, SLOT(display(milxQtModel*)));
        connect(plot, SIGNAL(resultAvailable(milxQtRenderWindow*)), this, SLOT(display(milxQtRenderWindow*)));

        if(success)
        {
            plot->setName(filename);
            plot->setConsole(console);
            plot->SetSource(table);
            plot->generatePlot(); //emits result
        }
    }
    else
    {
        printDebug("Opening Image.");
        QPointer<milxQtImage> img = new milxQtImage;  //list deletion

        printDebug("Supported Image formats: " + reader->supportedImageFormats());
        success = reader->openImage(nativeFilename, img);

        printInfo("Image Pixel Type: " + reader->getPixelType());
        printInfo("Image Component Type: " + reader->getComponentType());
        printInfo("Image Number of Components: " + QString::number(reader->getNumberOfComponents()));
        printInfo("Image Dimensions: " + QString::number(reader->getNumberOfDimensions()));

        if(success)
        {
            img->setName(filename);
            img->setConsole(console);
            img->generateImage();

            predisplay(img);
        }
        else
          printError("File format didn't appear to be supported. Check the file or add plugins for support.");
    }

    if(success)
        setCurrentFile(nativeFilename);

    return success;
}

void milxQtMain::save(QString filename)
{
    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString();
    QWidget *activeWindow = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow();
    milxQtWindow *currentWindowOpen = qobject_cast<milxQtWindow *>(activeWindow);
    bool pluginSave = false, success = false;

    if(!activeWindow)
        return;

    //Convert string to native paths
    filename = QDir::toNativeSeparators(filename);

    QFileDialog *fileSaver = new QFileDialog(this);

    ///Check if plugin's window
    foreach (QPointer<milxQtPluginInterface> loadedPlugin, plugins)
    {
        if(!loadedPlugin->hasSaveSupport())
            continue;

        if(loadedPlugin->isPluginWindow(activeWindow))
        {
            printInfo("Window is a " + loadedPlugin->name() + " plugin window");
            if(filename.isEmpty())
            {
                filename = fileSaver->getSaveFileName(this,
                                                      tr("Select File Name to Save"),
                                                      currentWindowOpen->getName(),
                                                      tr(loadedPlugin->saveFileSupport().toStdString().c_str()));
            }

            if(!filename.isEmpty())
            {
                loadedPlugin->save(filename);
                success = true;
            }
            pluginSave = true;
        }
    }

    if(!pluginSave)
    {
        if(isActiveImage())
        {
            printDebug("Saving Image with Built-in Capability");
            milxQtImage* image = qobject_cast<milxQtImage *>(activeWindow);

            if(filename.isEmpty())
            {
                filename = fileSaver->getSaveFileName(this,
                                                      tr("Select File Name to Save"),
                                                      image->getName(),
                                                      tr(saveExtsForImages.c_str()));
            }

            if(!filename.isEmpty())
            {
                image->setName(filename);
                QPointer<milxQtFile> writer = new milxQtFile; //Smart deletion
                writer->saveImage(filename, image);
                success = true;
            }
        }
        else if(isActivePlot())
        {
            printDebug("Saving Plot with Built-in Capability");
            milxQtPlot* plot = qobject_cast<milxQtPlot *>(activeWindow);

            if(filename.isEmpty())
            {
                filename = fileSaver->getSaveFileName(this,
                                                      tr("Select File Name to Save"),
                                                      plot->getName(),
                                                      tr(saveOtherExts.c_str()));
            }

            if(!filename.isEmpty())
            {
                plot->setName(filename);
                success = milx::File::SaveDelimitedText(filename.toStdString(), plot->GetSource());
            }
        }
        else if(isActiveModel())
        {
            printDebug("Saving Model with Built-in Capability");
            milxQtModel* model = qobject_cast<milxQtModel *>(activeWindow);

            if(filename.isEmpty())
            {
                filename = fileSaver->getSaveFileName(this,
                                                      tr("Select File Name to Save"),
                                                      model->getName(),
                                                      tr(saveExtsForModels.c_str()));
            }

            if(!filename.isEmpty())
            {
                model->setName(filename);
                QPointer<milxQtFile> writer = new milxQtFile; //Smart deletion
                writer->saveModel(filename, model, false);
                success = true;
            }
        }
        else
        {
            printError("Window not supported.");
            return;
        }
    }

    if(success)
    {
        setCurrentFile(filename);
        printInfo("Write Complete.");
    }
}

void milxQtMain::saveScreen(QString filename)
{
    QSettings settings("Shekhar Chandra", "milxQt");
    QString path = settings.value("recentPath").toString();
    QWidget *activeWindow = qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow();

    //Convert string to native paths
    filename = QDir::toNativeSeparators(filename);

    if(!activeWindow)
        return;
    else
    {
        QFileDialog *fileSaver = new QFileDialog(this);
        vtkSmartPointer<vtkWindowToImageFilter> windowToImage = vtkSmartPointer<vtkWindowToImageFilter>::New();
        QVTKWidget* windowVTK = qobject_cast<QVTKWidget *>(activeWindow);

        if(windowVTK == 0)
            return;

        //Ensure GUI doesnt interfere
//        windowVTK->GetRenderWindow()->StereoRenderOff();
//        windowVTK->GetRenderWindow()->OffScreenRenderingOn();
//        windowVTK->GetRenderWindow()->SetAlphaBitPlanes(1); //Use alpha channel
        windowVTK->GetRenderWindow()->Render();

        windowToImage->SetInput(windowVTK->GetRenderWindow());
        windowToImage->SetMagnification(magnifyFactor);
//        windowToImage->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
        windowToImage->ReadFrontBufferOff();
        windowToImage->Update();

        if(filename.isEmpty())
        {
            filename = fileSaver->getSaveFileName(this,
                               tr("Select File Name to Save"),
                               path,
                               tr(saveExtsForScreens.c_str()));
        }

        QPointer<milxQtFile> writer = new milxQtFile; //Smart deletion
        windowVTK->GetRenderWindow()->Render();

        //Save screenshot
        int extent[6];
        windowToImage->GetOutput()->GetExtent(extent);
        printDebug("Screenshot Size: " + QString::number(extent[1]-extent[0]) + ", " + QString::number(extent[3]-extent[2]) + ", " + QString::number(extent[5]-extent[4]));
        bool success = writer->saveImage(filename, windowToImage->GetOutput());

//        windowVTK->GetRenderWindow()->OffScreenRenderingOff();
        windowVTK->GetRenderWindow()->Render(); //Restore rendering

        if(!success)
        {
            printError("Unable to save screenshot. Ignoring.");
            return;
        }

        printInfo("Write Complete.");
    }
}

void milxQtMain::setTabName(QWidget *fromWindow)
{
    QString tabTitle = activeNamePrefix();
    int index = workspaces->currentIndex();

    if(fromWindow == NULL)
        workspaces->setTabText(index, "Empty");
    else
        workspaces->setTabText(index, tabTitle);
}

void milxQtMain::setTabName(const QString newName)
{
    int index = workspaces->currentIndex();

    workspaces->setTabText(index, newName);
}

void milxQtMain::closeTab(int index)
{
    int newIndex = 0;

    if(workspaces->count() > 1 || index > 0)
    {
        QWorkspace *tmpWorkspace = qobject_cast<QWorkspace *>(workspaces->widget(index));
        disconnect(tmpWorkspace, SIGNAL(windowActivated(QWidget *)), 0, 0);
        tmpWorkspace->closeAllWindows();
        tmpWorkspace->close();

        workspaces->removeTab(index);

        if(index > 0 && workspaces->currentIndex() == index)
            newIndex = workspaces->currentIndex();
        else if(index > 0)
            newIndex = index-1; //Safe to decrement
        workspaces->setCurrentIndex(newIndex); //else use zero
    }
}

void milxQtMain::tileTabVertically()
{
    QWorkspace *wrkSpc = qobject_cast<QWorkspace *>(workspaces->currentWidget());
    QWidgetList windows = wrkSpc->windowList();
    if (windows.count() < 2) {
        tileTab();
        return;
    }
    int wHeight = wrkSpc->height() / windows.count();
    int y = 0;
    foreach(QWidget *widget, windows)
    {
        widget->parentWidget()->resize(wrkSpc->width(), wHeight);
        widget->parentWidget()->move(0, y);
        y += wHeight;
    }
}

void milxQtMain::tileTabHorizontally()
{
    QWorkspace *wrkSpc = qobject_cast<QWorkspace *>(workspaces->currentWidget());
    QWidgetList windows = wrkSpc->windowList();
    if (windows.count() < 2) {
        tileTab();
        return;
    }
    int wWidth = wrkSpc->width() / windows.count();
    int x = 0;
    foreach(QWidget *widget, windows)
    {
        widget->parentWidget()->resize(wWidth, wrkSpc->height());
        widget->parentWidget()->move(x, 0);
        x += wWidth;
    }
}

void milxQtMain::helpContents()
{
    printDebug("Showing Help browser");

    QFile file(":/resources/smilx_doc/home.html");
    QWebView *view = new QWebView(this);
    if(file.open(QIODevice::ReadOnly))
      view->setHtml(file.readAll());
      view->setWindowTitle("sMILX Help");
      qobject_cast<QWorkspace *>(workspaces->currentWidget())->addWindow(view);
      view->show();

    //Quick setup toolbar
    QToolBar *toolBar = addToolBar(QObject::tr("Navigation"));
    toolBar->addAction(view->pageAction(QWebPage::Back));
    toolBar->addAction(view->pageAction(QWebPage::Forward));
    toolBar->addAction(view->pageAction(QWebPage::Reload));
    toolBar->addAction(view->pageAction(QWebPage::Stop));
}

void milxQtMain::preferences()
{
    ///Multi-page dialog is populated with options
    ///Upon acceptance, the settings are directly written
    ///to the MainWindow
    milxQtPreferencesForm prefsForm(this);

    prefsForm.exec();
}

void milxQtMain::controls()
{
    printDebug("Showing controls available...");
    QPixmap pixmap(":resources/controls_splash.png");
    QSplashScreen *controlsSplash = new QSplashScreen(this);
        controlsSplash->setPixmap(pixmap);
        controlsSplash->setMask(pixmap.mask());
        controlsSplash->show();

    qApp->processEvents();
}

void milxQtMain::about()
{
  milxQtAboutForm aboutForm(this);

  aboutForm.exec();
}

void milxQtMain::working(int value)
{
    if(value >= 0)
        progressBar->setValue(value);
    else
    {
        progressCallCount ++;
        progressBar->setMinimum(0);
        progressBar->setMaximum(0);
        printDebug("Working ... ");
    }
}

void milxQtMain::done(int value)
{
    if(value < 0)
        progressCallCount --;

    if(progressCallCount == 0)
    {
        progressBar->setMinimum(0);
        progressBar->setMaximum(100);
        progressBar->setValue(100);
//        progressBar->reset();
        printDebug("Done.");
    }
}

void milxQtMain::display(milxQtRenderWindow* newRender)
{
    ///Set the render into the viewer
    addRender(newRender);
    newRender->setToolTip("<p style='white-space:pre'>VTK Image Viewer 2 Keys - <b>f:</b> Move to, <b>r:</b> Reset, <b>Shift+r:</b> Reset Camera,\n<b>Shift+Mouse1:</b> Translate, <b>Shift+Ctrl+Mouse1:</b> Zoom</p>");
    newRender->setDefaultView(defaultViewBox->currentIndex());
    newRender->setView(defaultViewBox->currentIndex());
    newRender->setDefaultOrientation(defaultOrientationTypeBox->currentIndex());
    newRender->show();

    foreach(QAction *currAct, renderExtsActions) ///Add extension actions
    {
        newRender->addExtensionAction(currAct);
    }

    connect(this, SIGNAL(updatedImportFromMenu(QMenu*)), newRender, SLOT(createCustomConnections(QMenu*)));
    connect(this, SIGNAL(updatedImportFromMenu(QActionGroup*)), newRender, SLOT(setCustomActionGroup(QActionGroup*)));
    connect(newRender, SIGNAL(imageAvailable(vtkImageData*, QString )), this, SLOT(display(vtkImageData*, QString )));
    connect(newRender, SIGNAL(modelAvailable(vtkPolyData*, QString )), this, SLOT(display(vtkPolyData*, QString )));
    Connector->Connect(newRender->GetRenderWindow(),
                       vtkCommand::ModifiedEvent,
                       this,
                       SLOT( transferViewToWindows(vtkObject*, unsigned long, void*, void*, vtkCommand*) ),
                       NULL, 1.0); //High Priority

    ///Add common menus
    newRender->addToContextMenu(menuWindowList);
    newRender->addToContextMenu(importFromMenu);
    newRender->addToContextMenu(actionLinkWindows);

    ///Options
    newRender->background(whiteBackground);
    if(!humanGlyph)
        newRender->disableOrient();

    update();
    emit displayed(newRender);
}

void milxQtMain::predisplay(milxQtImage* newImage)
{
    const QString filename = newImage->getName();
    const size_t numberOfWindows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList().size();

    if(defaultViewTypeBox->currentIndex() != SINGLE && numberOfWindows > 0)
        newTab();
    else if(defaultViewTypeBox->currentIndex() == SINGLE)
    {
        newImage->setDefaultOrientation(defaultOrientationTypeBox->currentIndex()); //do not remove, not redundant
        display(newImage);
    }

    if(defaultViewTypeBox->currentIndex() != SINGLE) //Are we displaying scanner like three views + 3D view?
    {
        int ret = QMessageBox::Yes;
        if (actionLinkWindows->isChecked())
        {
            QMessageBox msgBox;
            msgBox.setText("Linked Views mode detected and has to be disabled for multi-view.");
            msgBox.setInformativeText("Do you want to disable link views and continue?");
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
            msgBox.setDefaultButton(QMessageBox::Yes);
            ret = msgBox.exec();
        }
        if (ret == QMessageBox::No)
            return;
        else if (actionLinkWindows->isChecked())
            actionLinkWindows->setChecked(false);

        //Axial view
        newImage->setNamePrefix("Axial: ");
        newImage->setName(filename);
        newImage->disableDefaultView();//ignore default view
        newImage->setDefaultOrientation(defaultOrientationTypeBox->currentIndex()); //do not remove, not redundant
        newImage->viewToAxial();
        display(newImage);

        //Sagittal view
        QPointer<milxQtImage> imgSagittal = new milxQtImage;  //list deletion
        imgSagittal->setDisplayData(newImage);
        imgSagittal->setNamePrefix("Sagittal: ");
        imgSagittal->setName(filename);
        imgSagittal->setConsole(console);
        imgSagittal->disableDefaultView();
        imgSagittal->generateImage();
        imgSagittal->setDefaultOrientation(defaultOrientationTypeBox->currentIndex()); //do not remove, not redundant
        imgSagittal->viewToSagittal();
        display(imgSagittal);

        //Coronal view
        QPointer<milxQtImage> imgCoronal = new milxQtImage;  //list deletion
        imgCoronal->setDisplayData(newImage);
        imgCoronal->setNamePrefix("Coronal: ");
        imgCoronal->setName(filename);
        imgCoronal->setConsole(console);
        imgCoronal->disableDefaultView();
        imgCoronal->generateImage();
        imgCoronal->setDefaultOrientation(defaultOrientationTypeBox->currentIndex()); //do not remove, not redundant
        imgCoronal->viewToCoronal();
        display(imgCoronal);

        //3D view
        QPointer<milxQtRenderWindow> slicesView = new milxQtRenderWindow;  //list deletion
        slicesView->setNamePrefix("3D View: ");
        slicesView->setName(filename);
        slicesView->setConsole(console);
        slicesView->addImageActor(newImage->GetImageActor(), newImage->getTransformMatrix());
        slicesView->addImageActor(imgSagittal->GetImageActor(), imgSagittal->getTransformMatrix());
        slicesView->addImageActor(imgCoronal->GetImageActor(), imgCoronal->getTransformMatrix());
        //slicesView->addActor(newImage->GetCursorActor(), newImage->getTransformMatrix());
        slicesView->generateRender();

        //setup tracking slices and crosshairs
        newImage->trackView(imgSagittal, AXIAL);
        newImage->trackView(imgCoronal, AXIAL);
        imgSagittal->trackView(newImage, SAGITTAL);
        imgSagittal->trackView(imgCoronal, SAGITTAL);
        imgCoronal->trackView(newImage, CORONAL);
        imgCoronal->trackView(imgSagittal, CORONAL);

        ///Here we use the Qt signals and slots directly as it was found that the VTK-Qt connector caused problems
        ///with the image actors.
        connect(newImage, SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), slicesView, SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        connect(imgSagittal, SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), slicesView, SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        connect(imgCoronal, SIGNAL(modified(vtkSmartPointer<vtkImageActor> )), slicesView, SLOT(updateImageActor(vtkSmartPointer<vtkImageActor>)));
        ///link the data properties too, one change updates others
        connect(newImage, SIGNAL(modified(QPointer<milxQtImage> )), imgSagittal, SLOT(updateDisplay(QPointer<milxQtImage>)));
        connect(newImage, SIGNAL(modified(QPointer<milxQtImage> )), imgCoronal, SLOT(updateDisplay(QPointer<milxQtImage>)));
        connect(imgSagittal, SIGNAL(modified(QPointer<milxQtImage> )), newImage, SLOT(updateDisplay(QPointer<milxQtImage>)));
        connect(imgSagittal, SIGNAL(modified(QPointer<milxQtImage> )), imgCoronal, SLOT(updateDisplay(QPointer<milxQtImage>)));
        connect(imgCoronal, SIGNAL(modified(QPointer<milxQtImage> )), newImage, SLOT(updateDisplay(QPointer<milxQtImage>)));
        connect(imgCoronal, SIGNAL(modified(QPointer<milxQtImage> )), imgSagittal, SLOT(updateDisplay(QPointer<milxQtImage>)));

        display(slicesView);

        updateQtEvents(); //ensure all complete before tiling
        tileTab();
  }
}

void milxQtMain::display(milxQtImage* newImage)
{
    if(!newImage->isLoaded())
    {
        printError("Image not loaded properly. Aborting Display.");
        return; //Check if model is valid
    }

    ///Set the image into the viewer
    addImage(newImage);
    newImage->setToolTip("<p style='white-space:pre'>VTK Image Viewer 2 Keys - <b>f:</b> Move to, <b>r:</b> Reset, <b>Shift+r:</b> Reset Camera,\n<b>Shift+Mouse1:</b> Translate, <b>Shift+Ctrl+Mouse1:</b> Zoom</p>");
    if(defaultViewTypeBox->currentIndex() == SINGLE)
    {
        newImage->setDefaultView(defaultViewBox->currentIndex());
        newImage->setView(defaultViewBox->currentIndex());
    }
    newImage->setDefaultOrientation(defaultOrientationTypeBox->currentIndex());
    newImage->show();

    foreach(QAction *currAct, imageExtsActions) ///Add extension actions
    {
        newImage->addExtensionAction(currAct);
    }

    ///Add common actions/menus
    newImage->addToContextMenu(menuWindowList);

    ///Connect inter-model connections
    connect(newImage, SIGNAL( imageToSurface(vtkSmartPointer<vtkImageData>, const float) ), this, SLOT( imageToSurface(vtkSmartPointer<vtkImageData>, const float) ));
    connect(newImage, SIGNAL( imageToPolyData(vtkSmartPointer<vtkImageData>) ), this, SLOT( imageToPolyData(vtkSmartPointer<vtkImageData>) ));
    connect(newImage, SIGNAL( imageToPseudoImage(vectorImageType::Pointer) ), this, SLOT( imageToPseudoImage(vectorImageType::Pointer) ));
    connect(newImage, SIGNAL( imageToVectorField(vectorImageType::Pointer, floatImageType::Pointer, int, float) ), this, SLOT( imageToVectorField(vectorImageType::Pointer, floatImageType::Pointer, int, float) ));
    connect(newImage, SIGNAL( imageToTensorField(vectorImageType::Pointer, floatImageType::Pointer, int, float) ), this, SLOT( imageToTensorField(vectorImageType::Pointer, floatImageType::Pointer, int, float) ));
    connect(newImage, SIGNAL( imageToStreamLines(vectorImageType::Pointer, floatImageType::Pointer) ), this, SLOT( imageToStreamLines(vectorImageType::Pointer, floatImageType::Pointer) ));
    connect(newImage, SIGNAL( imageToVolume(vtkSmartPointer<vtkImageData>, bool) ), this, SLOT( imageToVolume(vtkSmartPointer<vtkImageData>, bool) ));
    connect(newImage, SIGNAL( imageToPlot(vtkSmartPointer<vtkImageData>, int) ), this, SLOT( imageToPlot(vtkSmartPointer<vtkImageData>, int) ));
    connect(newImage, SIGNAL( tableToPlot(vtkSmartPointer<vtkTable>, QString) ), this, SLOT( tableToPlot(vtkSmartPointer<vtkTable>, QString) ));
    connect(newImage, SIGNAL(imageAvailable(vtkImageData*, QString )), this, SLOT(display(vtkImageData*, QString )));
    connect(newImage, SIGNAL(modelAvailable(vtkPolyData*, QString )), this, SLOT(display(vtkPolyData*, QString )));
    connect(newImage, SIGNAL(closing(QWidget *)), this, SLOT(cleanUpOnClose(QWidget *)));
    connect(this, SIGNAL(updatedImportFromMenu(QMenu*)), newImage, SLOT(createCustomConnections(QMenu*)));
    connect(this, SIGNAL(updatedImportFromMenu(QActionGroup*)), newImage, SLOT(setCustomActionGroup(QActionGroup*)));
    Connector->Connect(newImage->GetRenderWindow(),
                       vtkCommand::ModifiedEvent,
                       this,
                       SLOT( transferViewToWindows(vtkObject*, unsigned long, void*, void*, vtkCommand*) ),
                       NULL, 1.0); //High Priority

    newImage->addToContextMenu(actionCompare);
    newImage->addToContextMenu(menuWindowList);
    newImage->addToContextMenu(importFromMenu);
    newImage->addToContextMenu(actionLinkWindows);
    newImage->appendToContextMenu(actionSave);

    ///Options
    if(!humanGlyph)
        newImage->disableOrient();
    if(!interpolationImages)
        newImage->disableInterpolateDisplay();
    if(!orientationImages)
        newImage->disableApplyOrientDisplay();

    imageWindows.append(newImage);
    update();
    emit displayed(newImage);
}

void milxQtMain::display(vtkImageData* newImage, QString nameOfImage)
{
    printDebug("Displaying Raw ImageData");
    QPointer<milxQtImage> image = new milxQtImage; //List deletion
    image->setName(nameOfImage);
    if(defaultViewTypeBox->currentIndex() == SINGLE)
    {
        image->setDefaultView(defaultViewBox->currentIndex());
        image->setView(defaultViewBox->currentIndex());
    }
    image->SetInput(newImage);
    image->setDefaultOrientation(defaultOrientationTypeBox->currentIndex());
    image->generateImage();
    display(image);
}

void milxQtMain::display(milxQtModel* newModel)
{
    if(!newModel->isLoaded())
    {
        printError("Model/Surface not loaded properly. Aborting Display.");
        return; //Check if model is valid
    }

    ///Set the Model into the viewer
    addModel(newModel);
    newModel->setToolTip("<p style='white-space:pre'>VTK Model Viewer 2 Keys - <b>f:</b> Move to, <b>r:</b> Reset, <b>Shift+r:</b> Reset Camera,\n<b>Shift+Mouse1:</b> Translate, <b>Shift+Ctrl+Mouse1:</b> Zoom</p>");
    newModel->setDefaultView(defaultViewBox->currentIndex());
    newModel->setDefaultOrientation(defaultOrientationTypeBox->currentIndex());
    newModel->setView(defaultViewBox->currentIndex());
    newModel->refresh();
    if(newModel->GetScalars())
        newModel->colourMapToJet();
    newModel->show();

    foreach(QAction *currAct, modelExtsActions) ///Add extension actions
    {
        newModel->addExtensionAction(currAct);
    }

    ///Connect custom menus and actions
    //newModel->createCustomMenu( windowActionList(importFromMenu, true, false) ); //Virtual call to any model subclasses
    //connect(this, SIGNAL(updatedImportFromMenu(QMenu*)), newModel, SLOT(copyToContextMenu(QMenu*)));
    connect(this, SIGNAL(updatedImportFromMenu(QMenu*)), newModel, SLOT(createCustomConnections(QMenu*)));
    connect(this, SIGNAL(updatedImportFromMenu(QActionGroup*)), newModel, SLOT(setCustomActionGroup(QActionGroup*)));
    connect(newModel, SIGNAL(imageAvailable(vtkImageData*, QString )), this, SLOT(display(vtkImageData*, QString )));
    connect(newModel, SIGNAL(modelAvailable(vtkPolyData*, QString )), this, SLOT(display(vtkPolyData*, QString )));
    connect(newModel, SIGNAL(closing(QWidget *)), this, SLOT(cleanUpOnClose(QWidget *)));
    Connector->Connect(newModel->GetRenderWindow(),
                       vtkCommand::ModifiedEvent,
                       this,
                       SLOT( transferViewToWindows(vtkObject*, unsigned long, void*, void*, vtkCommand*) ),
                       NULL, 1.0); //High Priority

    ///Connect inter-model connections
    connect(newModel, SIGNAL( surfaceToImage(vtkSmartPointer<vtkPolyData>) ), this, SLOT( voxeliseSurface(vtkSmartPointer<vtkPolyData>) ));

    ///Add common actions
    newModel->addToContextMenu(actionCompare);
    newModel->addToContextMenu(menuWindowList);
    newModel->addToContextMenu(importFromMenu);
    newModel->addToContextMenu(actionLinkWindows);
    newModel->appendToContextMenu(actionSave);

    ///Options
    newModel->background(whiteBackground);
    if(!humanGlyph)
        newModel->disableOrient();
    if(interpolationModels)
        newModel->enableInterpolateDisplay();
    if(scalarBarModels)
        newModel->scaleDisplay(true);

    modelWindows.append(newModel);
    update();
    emit displayed(newModel);
}

void milxQtMain::display(vtkPolyData* newModel, QString nameOfModel)
{
    printDebug("Displaying Raw PolyData");
    QPointer<milxQtModel> mdl = new milxQtModel; //List deletion
    mdl->setName(nameOfModel);
    mdl->SetInput(newModel);
    mdl->generateModel();
    display(mdl);
}

void milxQtMain::display(milxQtUnifiedWindow* newUni)
{
    addUnifiedWindow(newUni);
    newUni->setToolTip("<p style='white-space:pre'>VTK Model Viewer 2 Keys - <b>f:</b> Move to, <b>r:</b> Reset, <b>Shift+r:</b> Reset Camera,\n<b>Shift+Mouse1:</b> Translate, <b>Shift+Ctrl+Mouse1:</b> Zoom</p>");
    newUni->setDefaultView(defaultViewBox->currentIndex());
    newUni->setView(defaultViewBox->currentIndex());
    newUni->setDefaultOrientation(defaultOrientationTypeBox->currentIndex());
    newUni->show();

    connect(newUni, SIGNAL(imageAvailable(milxQtImage *)), this, SLOT(display(milxQtImage *)), Qt::UniqueConnection);
    connect(newUni, SIGNAL(modelAvailable(milxQtModel *)), this, SLOT(display(milxQtModel *)), Qt::UniqueConnection);
//    connect(newUni, SIGNAL(modelAvailable(QWidget *)), workspaces->currentWidget(), SLOT(setActiveWindow(QWidget *))); ///\todo Temp Fix: model result is not strictly active

    ///Options
    newUni->background(whiteBackground);
    if(!humanGlyph)
        newUni->disableOrient();

    update();
    emit displayed(newUni);
}

void milxQtMain::display(vtkPolyDataCollection *modelCollection, QStringList &filenames)
{
    bool success = false;
    size_t n = modelCollection->GetNumberOfItems();

    ///Check plugins for collection support
    foreach (QPointer<milxQtPluginInterface> loadedPlugin, plugins)
    {
        if(loadedPlugin->hasCollectionSupport())
        {
            printInfo(loadedPlugin->name() + " plugin has collection support.");

            QMessageBox msgBox;
            msgBox.setText(loadedPlugin->name() + " plugin has collection support.");
            msgBox.setInformativeText("Do you want to load the collection with this plugin?");
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::NoToAll);
            msgBox.setDefaultButton(QMessageBox::Yes);
            int ret = msgBox.exec();

            if(ret == QMessageBox::No) ///User says no, move to next plugin etc.
                continue;
            else if(ret == QMessageBox::NoToAll)
                break;

            loadedPlugin->SetInputCollection(modelCollection, filenames);

            QPointer<milxQtRenderWindow> renWin = loadedPlugin->genericResult();
            if(renWin)
            {
                renWin->setName("Collection");
                renWin->generateRender();
                display(renWin);
                success = true;
            }

            QPointer<milxQtModel> model = loadedPlugin->modelResult();
            if(model)
            {
                model->setName("Collection Model");
                model->generateModel();
                display(model);
                success = true;
            }

            QPointer<milxQtImage> image = loadedPlugin->imageResult();
            if(image)
            {
                image->setName("Collection Image");
                image->generateImage();
                display(image);
                success = true;
            }
        }

        if(success)
        {
            loadedPlugin->update();
            break;
        }
    }

    if(success) //!< if Plugin loaded collection then return,
        return;

    ///Init Colours
    vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
    lookupTable->SetTableRange(0.0, n+1);
    lookupTable->Build();

    QPointer<milxQtRenderWindow> allModels = new milxQtRenderWindow; //hierarchical deletion
    QList< QPointer<milxQtModel> > alignedModels;
    modelCollection->InitTraversal(); //Must be called before traversing, or crashes
    for(size_t j = 0; j < n; j ++)
    {
        double colour[3];

        ///Get Colour
        lookupTable->GetColor(j, colour); //!< Pull colour for data

        ///Build model
        QPointer<milxQtModel> alignedModel = new milxQtModel; //smart deletion
        alignedModel->setName( QString::number(j) );
        alignedModel->setConsole(console);
        alignedModel->setDefaultView(defaultViewBox->currentIndex());
        alignedModel->setView(defaultViewBox->currentIndex());
        alignedModel->SetInput(modelCollection->GetNextItem());
        alignedModel->generatePoints(colour[0], colour[2], colour[1]);
        alignedModel->SetOpacity(0.2);

        alignedModels.append(alignedModel);

        allModels->AddActor(alignedModel->GetActor()); //!< Add data to general display
    }

    ///Ask for mean meshes
    QMessageBox msgBox;
    msgBox.setText("Open the mean data file?");
    msgBox.setInformativeText("Would you like to open the mean file?\nThis will be shown as a solid model.");
    msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Ok);
    int ret = msgBox.exec();

    if(ret == QMessageBox::Ok)
    {
        if(open())
            allModels->AddActor(activeModel()->GetActor()); //!< Add data to general display
    }

    allModels->setName("Collection");
    allModels->setConsole(console);
    allModels->generateRender();
    display(allModels);
}

//Protected
//From Qt Examples for Recent Files
void milxQtMain::setCurrentFile(const QString &fileName)
{
    QFileInfo fi(fileName);

    if(!fileName.isEmpty())
        setWindowTitle(tr("%1 - %2").arg(strippedName(fileName)).arg(tr("milxQt")));
    else
        setWindowTitle(tr("sMILX"));

    QSettings settings("Shekhar Chandra", "milxQt");
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > MaxRecentFiles)
        files.removeLast();

    settings.setValue("recentFileList", files);
    settings.setValue("recentPath",fi.absolutePath());

    foreach (QWidget *widget, QApplication::topLevelWidgets())
    {
        milxQtMain *mainWin = qobject_cast<milxQtMain *>(widget);
        if (mainWin)
            mainWin->updateRecentFileActions();
    }
}

//From Qt Examples for Recent Files
void milxQtMain::updateRecentFileActions()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    QStringList files = settings.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), (int)MaxRecentFiles);

    for (int i = 0; i < numRecentFiles; ++i)
    {
        QString text = tr("&%1 %2").arg(i + 1).arg(strippedName(files[i]));
        actionsRecentFile[i]->setText(text);
        actionsRecentFile[i]->setData(files[i]);
        actionsRecentFile[i]->setToolTip(files[i]);
        actionsRecentFile[i]->setStatusTip(files[i]);
        actionsRecentFile[i]->setShortcut("Ctrl+" + text.setNum( (i+1)%numRecentFiles ));
        actionsRecentFile[i]->setVisible(true);
    }
    for (int j = numRecentFiles; j < MaxRecentFiles; ++j)
        actionsRecentFile[j]->setVisible(false);

    actionRecentFileSeparator->setVisible(numRecentFiles > 0);
}

QActionGroup* milxQtMain::updateWindowMenu()
{
    menuWindows->clear();
    menuWindows->addAction(actionCascade);
    menuWindows->addAction(actionTile);
    menuWindows->addAction(actionTileVertically);
    menuWindows->addAction(actionTileHorizontally);
    menuWindows->addSeparator()->setText(tr("Docked Windows"));
    foreach(QAction *dockAct, dockActions)
    {
        menuWindows->addAction(dockAct);
    }
    menuWindows->addSeparator()->setText(tr("Windows"));

    return windowActionList(menuWindows, true, true);
}

QActionGroup* milxQtMain::updateWindowListMenu(bool applyMapper)
{
    menuWindowList->clear();

    printDebug("Update Window List Window");
    QActionGroup* grp = windowActionList(menuWindowList, true, applyMapper);
    emit updatedWindowListMenu(grp);
    emit updatedWindowListMenu(menuWindowList);

    return grp;
}

QActionGroup* milxQtMain::updateImportFromMenu(bool applyMapper)
{
    importFromMenu->clear();

    printDebug("Update Import From Window");
    QActionGroup* grp = windowActionList(importFromMenu, true, applyMapper);
    emit updatedImportFromMenu(grp);
    emit updatedImportFromMenu(importFromMenu);

    return grp;
}

void milxQtMain::updateWindowsWithValue(int value)
{
    initialiseWindowTraversal();

    while(currentWindow())
    {
        milxQtWindow *win = currentWindow();
        if(isImage(win))
        {
          milxQtImage *img = qobject_cast<milxQtImage *>(win);
          img->setLevel(100-value);
        }
        nextWindow();
    }

    imageLevelSlider->setStatusTip("Image Contrast at "+QString::number(value)+"%");
    imageLevelSlider->setToolTip(QString::number(value)+"%");
}

void milxQtMain::updateWindowsWithAutoLevel()
{
  initialiseWindowTraversal();

  while (currentWindow())
  {
    milxQtWindow *win = currentWindow();
    if (isImage(win))
    {
      milxQtImage *img = qobject_cast<milxQtImage *>(win);
      img->autoLevel();
    }
    nextWindow();
  }
}

void milxQtMain::updateWindowsWithRefresh()
{
    initialiseWindowTraversal();

    while (currentWindow())
    {
        milxQtRenderWindow *win = nextRenderWindow();
        if (isImage(win))
        {
            milxQtImage *img = qobject_cast<milxQtImage *>(win);
            img->refresh();
        }
        else
            win->refresh();
    }
}

void milxQtMain::updateWindowsWithCursors()
{
    initialiseWindowTraversal();

    while (currentWindow())
    {
        milxQtWindow *win = currentWindow();
        if (isImage(win))
        {
            milxQtImage *img = qobject_cast<milxQtImage *>(win);
            img->enableCrosshair();
        }
        nextWindow();
    }
}

void milxQtMain::updateWindowsWithView(int value)
{
    initialiseWindowTraversal();

    while(currentWindow())
    {
        milxQtWindow *win = currentWindow();
        milxQtRenderWindow *rndWin = qobject_cast<milxQtRenderWindow *>(win);
        rndWin->setView(value);
        nextWindow();
    }
}

void milxQtMain::updateWindowsWithViewType(int value)
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();

    if(windows.isEmpty())
        return;

    QMessageBox msgBox;
    if(value == 0)
        msgBox.setText("Window View type changed to Single");
    else
        msgBox.setText("Window View type changed to Multiple");
    msgBox.setInformativeText("This will take effect once you reload your data.");
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
}

void milxQtMain::updateWindowsWithViewOrientation(int value)
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();

    if(windows.isEmpty())
        return;

    QMessageBox msgBox;
    if(value == 0)
        msgBox.setText("View Orientation Convention changed to Radiological");
    else
        msgBox.setText("View Orientation Convention changed to Neurological");
    msgBox.setInformativeText("This will take effect once you reload your data.");
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
}

QActionGroup* milxQtMain::windowActionList(QMenu *menuForList, bool groupTogether, bool applyMapper)
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();
    milxQtWindow *win = NULL;
    QActionGroup *winGp = new QActionGroup(this);
    QString text;

    foreach(QWidget *currentWindow, windows)
    {
        if(!currentWindow)
            continue;

        win = qobject_cast<milxQtWindow *>(currentWindow);

        if(win == 0)
            continue;

        text = win->strippedNamePrefix();
//        text = tr("%1 %2").arg(j + 1).arg(childVPlt->strippedNamePrefix());

//        QAction *action  = menuForList->addAction(text);
        QAction *action  = new QAction(currentWindow); //parent to be used later
        action->setText(text);
        menuForList->addAction(action);
        if(groupTogether)
        {
            action->setCheckable(true);
            action->setChecked(currentWindow == qobject_cast<QWorkspace *>(workspaces->currentWidget())->activeWindow());
            winGp->addAction(action);
        }

        if(applyMapper)
        {
            connect(action, SIGNAL(triggered()), windowMapper, SLOT(map()));
            windowMapper->setMapping(action, currentWindow);
        }
    }

    return winGp;
}

int milxQtMain::getNumberOfImageWindows()
{
    int n = getNumberOfWindows();

    int m = 0;
    initialiseWindowTraversal();
    for(int j = 0; j < n; j ++)
    {
        milxQtWindow *win = nextWindow();
        if(isImage(win))
            m ++;
    }
    printDebug("Number of Image Windows: " + QString::number(m));

    return m;
}

int milxQtMain::getNumberOfModelWindows()
{
    int n = getNumberOfWindows();

    int m = 0;
    initialiseWindowTraversal();
    for(int j = 0; j < n; j ++)
    {
        milxQtWindow *win = nextWindow();
        if(isModel(win))
            m ++;
    }
    printDebug("Number of Model Windows: " + QString::number(m));

    return m;
}

milxQtWindow* milxQtMain::currentWindow()
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();
    milxQtWindow *win = NULL;

    if(windowIterator < windows.size())
        win = qobject_cast<milxQtWindow *>(windows[windowIterator]);

    return win;
}

milxQtWindow* milxQtMain::nextWindow()
{
  QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();
  milxQtWindow *win = NULL;

  if(windowIterator < windows.size())
    {
      win = qobject_cast<milxQtWindow *>(windows[windowIterator]);
      windowIterator ++;
    }

  return win;
}

milxQtRenderWindow* milxQtMain::nextRenderWindow()
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();
    milxQtRenderWindow *win = NULL;

    while(windowIterator < windows.size())
    {
        if(isRender(windows[windowIterator]))
        {
            win = qobject_cast<milxQtRenderWindow *>(windows[windowIterator]);
            windowIterator ++;
            break;
        }

        windowIterator ++;
    }

    return win;
}

milxQtModel* milxQtMain::nextModel()
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();
    milxQtModel *win = NULL;

    while(windowIterator < windows.size())
    {
        if(isModel(windows[windowIterator]))
        {
            win = qobject_cast<milxQtModel *>(windows[windowIterator]);
            windowIterator ++;
            break;
        }

        windowIterator ++;
    }

    return win;
}

milxQtImage* milxQtMain::nextImage()
{
    QWidgetList windows = qobject_cast<QWorkspace *>(workspaces->currentWidget())->windowList();
    milxQtImage *win = NULL;

    printDebug("Number of Windows: " + QString::number(windows.size()));
    while(windowIterator < windows.size())
    {
        if(isImage(windows[windowIterator]))
        {
            win = qobject_cast<milxQtImage *>(windows[windowIterator]);
            windowIterator ++;
            break;
        }

        windowIterator ++;
    }
    printDebug("Window Iterator Value: " + QString::number(windowIterator));

    return win;
}

void milxQtMain::transferViewToWindows(vtkObject *obj, unsigned long, void *client_data, void *, vtkCommand *command)
{
    if(!actionLinkWindows->isChecked()) //if all windows not linked, exit
        return;

    printDebug("Updating Views in Other Windows");
    actionLinkWindows->setChecked(false); //prevent cycle calls
    // get render window
    vtkRenderWindow* rndWindow = vtkRenderWindow::SafeDownCast(obj);
    vtkRenderer* rnd = rndWindow->GetRenderers()->GetFirstRenderer();

    ///Copy camera of similar windows
    bool srcIsImage = isActiveImage();
    bool srcIsModel = isActiveModel();

    if(rnd)
    {
        vtkCamera* srcCamera = rnd->GetActiveCamera(); //!< Get callers camera

        initialiseWindowTraversal();

        for(int j = 0; j < getNumberOfWindows(); j ++) //!< For all windows, copy camera
        {
            QPointer<milxQtRenderWindow> window = nextRenderWindow();
            vtkRenderer* distRnd = window->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
            vtkCamera* distCamera = distRnd->GetActiveCamera();

            if(distCamera != srcCamera) //!< Copy camera over
            {
                if(srcIsImage && isImage(window))
                {
                    //Image-to-Image specific
                    QPointer<milxQtImage> img = qobject_cast<milxQtImage *>(window);
                    img->setSlice(activeImage()->getSlice());
                    img->setView(activeImage()->getView());
                    if(img->isCrosshair() && activeImage()->isCrosshair())
                        img->setCrosshairPosition(activeImage()->getCrosshairPosition());
                }
                /*else if(srcIsImage && isModel(window))
                {
                    //Image-to-Model specific
                    QPointer<milxQtModel> mdl = qobject_cast<milxQtModel *>(window);
                    mdl->setView(activeImage()->getView());
                }*/
                /*else if(srcIsModel && isImage(window))
                {
                    //Model-to-Image specific
                    QPointer<milxQtImage> img = qobject_cast<milxQtImage *>(window);
                    img->setView(activeModel()->getView());
                }*/ //Copy Camera is better

                if( !(srcIsImage && isModel(window)) && !(srcIsModel && isImage(window))) //no image-model
                    distCamera->DeepCopy(srcCamera);
            }

            window->Render(); //!< refresh the destination display, quick
        }

//        rndWindow->Render(); //!< refresh the source display, quick
    }

    actionLinkWindows->setChecked(true); //restore
}

void milxQtMain::dataMenu()
{
    milxQtRenderWindow *renWin = activeRender();

    if(renWin)
        renWin->createMenu(menuData);
}

void milxQtMain::imageToSurface(vtkSmartPointer<vtkImageData> img, const float value)
{
    //Check if value provided
    int contourNumber = 1;
    if(value == numeric_limits<float>::max())
        contourNumber = -1;
    else
        printInfo("Using iso value provided");

    QPointer<milxQtModel> model = new milxQtModel; //List deletion
        model->setConsole(console);
        model->generateIsoSurface(img, contourNumber, value);
        model->generateModel();
    display(model);
}

void milxQtMain::imageToPolyData(vtkSmartPointer<vtkImageData> img)
{
    QPointer<milxQtModel> model = new milxQtModel; //List deletion
        model->setConsole(console);
        model->generatePolyDataFromImage(img);
        model->generateModel();
    display(model);
}

void milxQtMain::imageToPseudoImage(vectorImageType::Pointer img)
{
//    typedef itk::VectorImage<floatPixelType, milx::imgDimension> vectorImageType;

    ///Setup vector field
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
//    vectorImageType::Pointer rescaledImg = milx::Image<vectorImageType>::RescaleIntensities(img, 0, 255);
    vtkSmartPointer<vtkImageData> fieldRaw = milx::Image<vectorImageType>::ConvertITKVectorImageToVTKImage(img);
    vtkSmartPointer<vtkImageData> field = milx::Image<vectorImageType>::ApplyOrientationToVTKImage(fieldRaw, img, matrix, false, true);

//    field->GetPointData()->SetActiveScalars("ImageScalars");
    field->GetPointData()->SetActiveVectors("ImageScalars");

    QPointer<milxQtImage> pseudoimage = new milxQtImage; //List deletion
        pseudoimage->setConsole(console);
        pseudoimage->SetInput(field);
        pseudoimage->generateImage();
    display(pseudoimage);
}

void milxQtMain::imageToVectorField(vectorImageType::Pointer img, floatImageType::Pointer magImg, int subsampleFactor, float scaling)
{
    bool ok1 = false;
    if(subsampleFactor == 0)
    {
        subsampleFactor = QInputDialog::getInt(this, tr("Please Provide the sub-sample factor of the data"),
                                              tr("Sub-sample Factor:"), 8, 1, 1000, 1, &ok1);
        if(!ok1)
            return;
    }

    ///Subsample image to reduce computation costs
    floatImageType::SizeType subsampleSizes;
        subsampleSizes.Fill(subsampleFactor);

    vectorImageType::Pointer imgSubSampled = milx::Image<vectorImageType>::SubsampleImage(img, subsampleSizes);
    floatImageType::Pointer magImgSubSampled = milx::Image<floatImageType>::SubsampleImage(magImg, subsampleSizes);
    std::cout << "Subsampled size: " << imgSubSampled->GetLargestPossibleRegion().GetSize() << std::endl;

    ///Need to flip y-axis because of VTK's CG coordinate system
//    vectorImageType::Pointer imgRealSpace = milx::Image<vectorImageType>::FlipImage(imgSubSampled, false, true, false);
//    floatImageType::Pointer magImgRealSpace = milx::Image<floatImageType>::FlipImage(magImgSubSampled, false, true, false);

    ///Create Polydata to hold vector field
    vtkSmartPointer<vtkPoints> fieldPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> magnitudes = vtkSmartPointer<vtkFloatArray>::New();
        magnitudes->SetNumberOfComponents(1);
    vtkSmartPointer<vtkFloatArray> vectors = vtkSmartPointer<vtkFloatArray>::New();
        vectors->SetNumberOfComponents(3);

    typedef itk::Point<double, 3> InputImagePointType;
    itk::ImageRegionConstIteratorWithIndex<vectorImageType> imageIterator(imgSubSampled, imgSubSampled->GetLargestPossibleRegion());
    itk::ImageRegionConstIteratorWithIndex<floatImageType> magImageIterator(magImgSubSampled, magImgSubSampled->GetLargestPossibleRegion());

    InputImagePointType point;
    while(!imageIterator.IsAtEnd())
    {
      coordinate position, direction;
      float value;

      //std::cout << "Index: " << imageIterator.GetIndex() << " value: " << imageIterator.Get() << std::endl;
      magImgSubSampled->TransformIndexToPhysicalPoint(magImageIterator.GetIndex(), point);

      position[0] = point[0];
      position[1] = point[1];
      position[2] = point[2];

      direction.fill(0.0);
      if(magImageIterator.Get() != 0) //if not masked
      {
        direction[0] = imageIterator.Get()[0];
        direction[1] = imageIterator.Get()[1];
        if(img->GetNumberOfComponentsPerPixel() > 2)
            direction[2] = imageIterator.Get()[2];
      }

      fieldPoints->InsertNextPoint(position.data_block());

      value = magImageIterator.Get();
      magnitudes->InsertNextTuple1(value);
      vectors->InsertNextTuple3(direction[0], direction[1], direction[2]);

      ++imageIterator;
      ++magImageIterator;
    }

    vtkSmartPointer<vtkPolyData> field = vtkSmartPointer<vtkPolyData>::New();
        field->SetPoints(fieldPoints);
        field->GetPointData()->SetScalars(magnitudes);
        field->GetPointData()->SetVectors(vectors);
        field->Modified();

    QPointer<milxQtModel> model = new milxQtModel; //List deletion
        model->setConsole(console);
        model->SetInput(field);
        model->generateModel();
        model->generateVectorField(scaling);

    display(model);
}

void milxQtMain::imageToTensorField(vectorImageType::Pointer img, floatImageType::Pointer magImg, int subsampleFactor, float scaling)
{
    bool ok1 = false;
    if(subsampleFactor == 0)
    {
        subsampleFactor = QInputDialog::getInt(this, tr("Please Provide the sub-sample factor of the data"),
                                              tr("Sub-sample Factor:"), 8, 1, 1000, 1, &ok1);
        if(!ok1)
            return;
    }

    ///Subsample image to reduce computation costs
    const size_t components = img->GetNumberOfComponentsPerPixel();
    floatImageType::SizeType subsampleSizes;
        subsampleSizes.Fill(subsampleFactor);

    vectorImageType::Pointer imgSubSampled = milx::Image<vectorImageType>::SubsampleImage(img, subsampleSizes);
    floatImageType::Pointer magImgSubSampled = milx::Image<floatImageType>::SubsampleImage(magImg, subsampleSizes);
    std::cout << "Subsampled size: " << imgSubSampled->GetLargestPossibleRegion().GetSize() << std::endl;

    ///Need to flip y-axis because of VTK's CG coordinate system
//    vectorImageType::Pointer imgRealSpace = milx::Image<vectorImageType>::FlipImage(imgSubSampled, false, true, false);
//    floatImageType::Pointer magImgRealSpace = milx::Image<floatImageType>::FlipImage(magImgSubSampled, false, true, false);

    ///Create Polydata to hold vector field
    vtkSmartPointer<vtkPoints> fieldPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> magnitudes = vtkSmartPointer<vtkFloatArray>::New();
        magnitudes->SetNumberOfComponents(1);
    vtkSmartPointer<vtkFloatArray> tensors = vtkSmartPointer<vtkFloatArray>::New();
        tensors->SetNumberOfComponents(9);

    typedef itk::Point<double, 3> InputImagePointType;
    itk::ImageRegionConstIteratorWithIndex<vectorImageType> imageIterator(imgSubSampled, imgSubSampled->GetLargestPossibleRegion());
    itk::ImageRegionConstIteratorWithIndex<floatImageType> magImageIterator(magImgSubSampled, magImgSubSampled->GetLargestPossibleRegion());

    InputImagePointType point;
    while(!imageIterator.IsAtEnd())
    {
      double position[3];
      float value;

      //std::cout << "Index: " << imageIterator.GetIndex() << " value: " << imageIterator.Get() << std::endl;
      magImgSubSampled->TransformIndexToPhysicalPoint(magImageIterator.GetIndex(), point);

      position[0] = point[0];
      position[1] = point[1];
      position[2] = point[2];

      vtkSmartPointer<vtkTensor> tens = vtkSmartPointer<vtkTensor>::New();
      if(magImageIterator.Get() != 0) //if not masked
      {
          if(components == 6)
          {
              tens->SetComponent(0,0, imageIterator.Get()[0]);
              tens->SetComponent(1,1, imageIterator.Get()[1]);
              tens->SetComponent(2,2, imageIterator.Get()[2]);
              tens->SetComponent(0,1, imageIterator.Get()[3]); //symmetric
              tens->SetComponent(1,0, imageIterator.Get()[3]); //symmetric
              tens->SetComponent(0,2, imageIterator.Get()[4]); //symmetric
              tens->SetComponent(2,0, imageIterator.Get()[4]); //symmetric
              tens->SetComponent(1,2, imageIterator.Get()[5]); //symmetric
              tens->SetComponent(2,1, imageIterator.Get()[5]); //symmetric
          }
          else
          {
              tens->SetComponent(0,0, imageIterator.Get()[0]);
              tens->SetComponent(1,1, imageIterator.Get()[1]);
              tens->SetComponent(2,2, imageIterator.Get()[2]);
              tens->SetComponent(0,1, imageIterator.Get()[3]);
              tens->SetComponent(1,0, imageIterator.Get()[4]);
              tens->SetComponent(0,2, imageIterator.Get()[5]);
              tens->SetComponent(2,0, imageIterator.Get()[6]);
              tens->SetComponent(1,2, imageIterator.Get()[7]);
              tens->SetComponent(2,1, imageIterator.Get()[8]);
          }
      }
      else //zero vector
      {
          tens->SetComponent(0,0, 0.0);
          tens->SetComponent(1,1, 0.0);
          tens->SetComponent(2,2, 0.0);
          tens->SetComponent(0,1, 0.0); //symmetric
          tens->SetComponent(1,0, 0.0); //symmetric
          tens->SetComponent(0,2, 0.0); //symmetric
          tens->SetComponent(2,0, 0.0); //symmetric
          tens->SetComponent(1,2, 0.0); //symmetric
          tens->SetComponent(2,1, 0.0); //symmetric
      }

      fieldPoints->InsertNextPoint(position);

      value = magImageIterator.Get();
      magnitudes->InsertNextTuple1(value);
      tensors->InsertNextTuple(tens->T);

      ++imageIterator;
      ++magImageIterator;
    }

    vtkSmartPointer<vtkPolyData> field = vtkSmartPointer<vtkPolyData>::New();
        field->SetPoints(fieldPoints);
        field->GetPointData()->SetScalars(magnitudes);
        field->GetPointData()->SetTensors(tensors);

    QPointer<milxQtModel> model = new milxQtModel; //List deletion
        model->setConsole(console);
        model->SetInput(field);
        model->generateModel();
        model->generateTensorField(scaling);
    display(model);
}

void milxQtMain::imageToStreamLines(vectorImageType::Pointer img, floatImageType::Pointer magImg, size_t subsampleFactor)
{
    bool ok1 = false;
    if(subsampleFactor == 0)
    {
        subsampleFactor = QInputDialog::getInt(this, tr("Please Provide the sub-sample factor of the data"),
                                              tr("Sub-sample Factor:"), 8, 1, 1000, 1, &ok1);
        if(!ok1)
            return;
    }

    //Which dimension is sliced?
    floatImageType::SizeType imageSize = magImg->GetLargestPossibleRegion().GetSize();
    size_t sliceDimension = 0;
    for(size_t j = 0; j < floatImageType::ImageDimension; j ++)
    {
      if(imageSize[j] == 1)
        sliceDimension = j;
    }

    ///Check mag image if its a 2D slice
    const size_t sliceSubsampleFactor = subsampleFactor;
    floatImageType::SizeType subsampleSizes;
        subsampleSizes.Fill(subsampleFactor);
    floatImageType::SizeType sliceSubsampleSizes;
        sliceSubsampleSizes.Fill(sliceSubsampleFactor);
        sliceSubsampleSizes[sliceDimension] = 1;

    ///Subsample image to reduce computation costs
//    const size_t components = img->GetNumberOfComponentsPerPixel();
    vectorImageType::Pointer imgSubSampled = milx::Image<vectorImageType>::SubsampleImage(img, subsampleSizes);
    floatImageType::Pointer magImgSubSampled = milx::Image<floatImageType>::SubsampleImage(magImg, sliceSubsampleSizes);
    cout << "Slice Information: " << endl;
    milx::Image<floatImageType>::Information(magImg);

    ///Need to flip y-axis because of VTK's CG coordinate system
//    vectorImageType::Pointer imgRealSpace = milx::Image<vectorImageType>::FlipImage(imgSubSampled, false, true, false);
//    floatImageType::Pointer magImgRealSpace = milx::Image<floatImageType>::FlipImage(magImgSubSampled, false, true, false);

    ///Create Polydata to hold vector field
    printDebug("Setting up seed points for streamlines using the current slice");
    vtkSmartPointer<vtkPoints> slicePoints = vtkSmartPointer<vtkPoints>::New();

    typedef itk::Point<double, 3> InputImagePointType;
    typedef itk::Point<double, 2> InputSlicePointType;
    itk::ImageRegionConstIteratorWithIndex<vectorImageType> imageIterator(imgSubSampled, imgSubSampled->GetLargestPossibleRegion());
    itk::ImageRegionConstIteratorWithIndex<floatImageType> magImageIterator(magImgSubSampled, magImgSubSampled->GetLargestPossibleRegion());

    ///Setup seed points
    InputImagePointType point;
    while(!magImageIterator.IsAtEnd())
    {
        double position[3];

        magImgSubSampled->TransformIndexToPhysicalPoint(magImageIterator.GetIndex(), point);

        position[0] = point[0];
        position[1] = point[1];
        position[2] = point[2];

        slicePoints->InsertNextPoint(position);

        ++magImageIterator;
    }

    vtkSmartPointer<vtkPolyData> slice = vtkSmartPointer<vtkPolyData>::New();
        slice->SetPoints(slicePoints);

    ///Setup vector field
    printDebug("Converting ITK Vector Image to VTK Vector Image");
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkImageData> fieldRaw = milx::Image<vectorImageType>::ConvertITKVectorImageToVTKImage(imgSubSampled);
    printDebug("Applying orientation to VTK Vector Image");
    if(!fieldRaw || fieldRaw->GetNumberOfScalarComponents() < 3)
        printWarning("Possible problem with the VTK Vector image");
    vtkSmartPointer<vtkImageData> field = milx::Image<vectorImageType>::ApplyOrientationToVTKImage(fieldRaw, imgSubSampled, matrix, false, false);

    field->GetPointData()->SetActiveVectors("ImageScalars");

    printDebug("Computing the streamlines");
    QPointer<milxQtModel> model = new milxQtModel; //List deletion
        model->setConsole(console);
        model->SetInput(slice);
        model->generateStreamLines(field);
    display(model);
}

void milxQtMain::imageToVolume(vtkSmartPointer<vtkImageData> img, bool eightbit)
{
    printDebug("Attempting Volume Rendering");
    QPointer<milxQtPlot> plot = new milxQtPlot; //List deletion
        plot->setName("Volume Rendering");
        plot->setConsole(console);
        plot->setPlotTypeVolume();
        plot->SetSource(img, eightbit, false);
        plot->generatePlot();//emits result, not connected
    printDebug("Successful. Displaying...");
    display( qobject_cast<milxQtRenderWindow *>(plot) ); //volume not really model so need rw
}

void milxQtMain::imageToPlot(vtkSmartPointer<vtkImageData> img, int displaceAxis)
{
    printDebug("Attempting Surface plot of image data");
    QPointer<milxQtPlot> plot = new milxQtPlot; //List deletion
        plot->setName("Image Surface");
        plot->setConsole(console);
        plot->setPlotTypeSurface(displaceAxis);
        plot->SetSource(img, false, true);
        plot->generatePlot();//emits result, not connected
    printDebug("Successful. Displaying...");
    display(plot);
}

void milxQtMain::tableToPlot(vtkSmartPointer<vtkTable> tbl, QString title)
{
    printDebug("Attempting plot of table data");
    QPointer<milxQtPlot> plot = new milxQtPlot; //List deletion
        plot->setName(title);
        plot->setConsole(console);
        plot->SetSource(tbl);
        plot->generatePlot();//emits result, not connected
    printDebug("Successful. Displaying...");

    if(tbl->GetNumberOfColumns() < 10)
        display( qobject_cast<milxQtRenderWindow *>(plot) );
    else //surface plot
        display(plot);
}

void milxQtMain::voxeliseSurface(vtkSmartPointer<vtkPolyData> surface)
{
    QPointer<milxQtImage> image = new milxQtImage; //List deletion
        image->generateVoxelisedSurface(surface);
        image->setConsole(console);
        image->generateImage();
        display(image);
}

void milxQtMain::imagesMix()
{
    int n = getNumberOfImageWindows();
    if(getNumberOfWindows() < 2 || n < 2)
    {
        printError("Need more than 1 window open to blend.");
        return;
    }

    initialiseWindowTraversal();
    QPointer<milxQtImage> firstImg = nextImage();

    //Create slicers and the dialog
    QVector< QSlider* > sliders;
    QVector< QCheckBox* > checkBoxes;
    QSlider *opacitySldr = new QSlider(this);
    opacitySldr->setMinimum(0);
    opacitySldr->setMaximum(10);
    opacitySldr->setValue(10);
    QCheckBox *initialCheckbox = new QCheckBox(this);
    initialCheckbox->setChecked(true);
//    initialCheckbox->setDisabled(true);
    initialCheckbox->setToolTip(firstImg->strippedBaseName());
    QVBoxLayout *initialLayout = new QVBoxLayout(this);
    initialLayout->addWidget(initialCheckbox);
    initialLayout->addWidget(opacitySldr);
    QHBoxLayout *sliderLayout = new QHBoxLayout(this);
    sliderLayout->addLayout(initialLayout);
    sliders.push_back(opacitySldr);
    checkBoxes.push_back(initialCheckbox);

    if(firstImg->GetLookupTable() == NULL)
    {
        QMessageBox msgBox;
        msgBox.setText("Colormap not set?");
        msgBox.setInformativeText("Image " + firstImg->strippedBaseName() + " has no colormap set. Setting to a default map.");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();

		//Use default cmap
		if (firstImg->is8BitImage())
			firstImg->colourMapToHSV();
		else
			firstImg->colourMapToGray();
    }

    for(int j = 1; j < n; j ++) //!< For all windows, do operation
    {
        QPointer<milxQtImage> secondImg = nextImage();

        if(secondImg->GetLookupTable() == NULL)
        {
            QMessageBox msgBox;
            msgBox.setText("Colormap not set?");
            msgBox.setInformativeText("Image " + secondImg->strippedBaseName() + " has no colormap set. Setting to a default map.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.setDefaultButton(QMessageBox::Ok);
            msgBox.exec();

			//Use default cmap
			if (secondImg->is8BitImage() || secondImg->is32BitImage()) //follow up images expected to be a label
				secondImg->colourMapToHSV();
			else
				secondImg->colourMapToGray(); 
        }

        QVBoxLayout *layout = new QVBoxLayout(this);
        QSlider *opacitySldr2 = new QSlider(this);
        opacitySldr2->setMinimum(0);
        opacitySldr2->setMaximum(10);
        opacitySldr2->setValue(5);
        QCheckBox *checkbox = new QCheckBox(this);
        checkbox->setChecked(true);
        checkbox->setToolTip(secondImg->strippedBaseName());
        layout->addWidget(checkbox);
        layout->addWidget(opacitySldr2);
        sliderLayout->addLayout(layout);
        sliders.push_back(opacitySldr2);
        checkBoxes.push_back(checkbox);
    }
//    QVBoxLayout *btnLayout = new QVBoxLayout(slidersDlg);
    QHBoxLayout *dlgLayout = new QHBoxLayout(this);
    QVBoxLayout *btnlayout = new QVBoxLayout(this);
    QPushButton *btnOK = new QPushButton("OK", this);
    QPushButton *btnCancel = new QPushButton("Cancel", this);
    btnlayout->addWidget(btnOK);
    btnlayout->addWidget(btnCancel);
    dlgLayout->addLayout(sliderLayout);
    dlgLayout->addLayout(btnlayout);

    QDialog *slidersDlg = new QDialog(this);
    slidersDlg->setLayout(dlgLayout);
    connect(btnOK, SIGNAL(clicked()), slidersDlg, SLOT(accept()));
    connect(btnCancel, SIGNAL(clicked()), slidersDlg, SLOT(reject()));
    int ret = slidersDlg->exec();

    if(ret == QDialog::Rejected)
        return;

    QVector<float> opacities;
    for(size_t j = 0; j < sliders.size(); j ++)
    {
        float value = sliders[j]->value()*0.1;
        if(!checkBoxes[j]->isChecked())
            value = 0.0;
        opacities.push_back(value);
    }

    imagesBlend(opacities);
}

void milxQtMain::imagesBlend(QVector<float> opacities)
{
    int n = getNumberOfImageWindows();
    if(getNumberOfWindows() < 2 || n < 2)
    {
        printError("Need more than 1 window open to blend.");
        return;
    }

    if(opacities.empty() || opacities.size() < n)
    {
        printError("Size of opacities list does not match number of open images.");
        return;
    }

    initialiseWindowTraversal();
    QPointer<milxQtImage> firstImg = nextImage();
    //vtkSmartPointer<vtkImageData> ucharData1 = firstImg->GetWindowLevel()->GetOutput();
    vtkSmartPointer<vtkImageMapToColors> filterColorsImage = vtkSmartPointer<vtkImageMapToColors>::New();
    filterColorsImage->SetLookupTable(firstImg->GetLookupTable());
#if VTK_MAJOR_VERSION <= 5
    filterColorsImage->SetInput(firstImg->GetOutput());
#else
    filterColorsImage->SetInputData(firstImg->GetOutput());
#endif
    filterColorsImage->PassAlphaToOutputOn();
    filterColorsImage->Update();
    vtkSmartPointer<vtkImageData> ucharData1 = filterColorsImage->GetOutput();

    int initialExtent[6];
    firstImg->GetOutput()->GetExtent(initialExtent);

    // Combine the images (blend takes multiple connections on the 0th input port)
    emit working(-1);
    vtkSmartPointer<vtkImageBlend> blend = vtkSmartPointer<vtkImageBlend>::New();
        linkProgressEventOf(blend);
        blend->SetOpacity(0,opacities[0]);
    #if VTK_MAJOR_VERSION <= 5
        blend->AddInput(ucharData1);
    #else
        blend->AddInputData(ucharData1);
    #endif
//        blend->SetBlendModeToCompound();
        blend->SetBlendModeToNormal();

    for(int j = 1; j < n; j ++) //!< For all windows, do operation
    {
        QPointer<milxQtImage> secondImg = nextImage();

        if(opacities[j] == 0.0)
            continue;

        printInfo("Blending with Image: " + secondImg->strippedName() + " with Opacity: " + QString::number(opacities[j]));
//        vtkSmartPointer<vtkImageData> ucharData2 = secondImg->GetWindowLevel()->GetOutput();
        vtkSmartPointer<vtkImageMapToColors> filterColorsOverlay = vtkSmartPointer<vtkImageMapToColors>::New();
        filterColorsOverlay->SetLookupTable(secondImg->GetLookupTable());
     #if VTK_MAJOR_VERSION <= 5
        filterColorsOverlay->SetInput(secondImg->GetOutput());
     #else
        filterColorsOverlay->SetInputData(secondImg->GetOutput());
     #endif
        filterColorsOverlay->PassAlphaToOutputOn();
        filterColorsOverlay->Update();
        vtkSmartPointer<vtkImageData> ucharData2 = filterColorsOverlay->GetOutput();

        if(!secondImg->GetLookupTable())
            printWarning("Colourmap is not set. Please set a colour map to ensure proper blending.");

        int actualExtent[6];
        secondImg->GetOutput()->GetExtent(actualExtent);

        if(initialExtent[1] == actualExtent[1] && initialExtent[3] == actualExtent[3] && initialExtent[5] == actualExtent[5])
        {
        #if VTK_MAJOR_VERSION <= 5
            blend->AddInput(ucharData2);
        #else
            blend->AddInputData(ucharData2);
        #endif
            blend->SetOpacity(j,opacities[j]);
        }
        else
            printError("Images are not the same size. Skipping.");
    }

    printInfo("Blending");
    blend->Update();
    printDebug("Number of components: " + QString::number(blend->GetOutput()->GetNumberOfScalarComponents()));
    emit done(-1);

    QPointer<milxQtImage> blendResult = new milxQtImage;
        blendResult->setName("Blended Images");
        blendResult->setConsole(console);
        blendResult->setData(blend->GetOutput());
        blendResult->generateImage();

    display(blendResult);
}

void milxQtMain::imagesAdd()
{
    if(getNumberOfWindows() == 0 || imageWindows.size() < 1)
        return;

    emit working(-1);
    initialiseWindowTraversal();
    QPointer<milxQtImage> firstImg = nextImage();

    QPointer<milxQtImage> resultImg = new milxQtImage;
    printInfo("Assigning " + firstImg->getName());
    resultImg->setName("Sum Image");
    resultImg->setConsole(console);
    resultImg->setData(firstImg, true);
    for(int j = 1; j < imageWindows.size(); j ++) //!< For all windows, do operation
    {
        QPointer<milxQtImage> img = nextImage();

        printInfo("Adding " + img->getName());
        resultImg->add(img);
    }
    emit done(-1);

    display(resultImg);
}

void milxQtMain::imagesSubtract()
{
    if(getNumberOfWindows() == 0 || imageWindows.size() < 1)
        return;

    emit working(-1);
    initialiseWindowTraversal();
    QPointer<milxQtImage> firstImg = nextImage();

    QPointer<milxQtImage> resultImg = new milxQtImage;
    printInfo("Assigning " + firstImg->getName());
    resultImg->setName("Difference Image");
    resultImg->setConsole(console);
    resultImg->setData(firstImg, true);
    for(int j = 1; j < imageWindows.size(); j ++) //!< For all windows, do operation
    {
        QPointer<milxQtImage> img = nextImage();

        printInfo("Subtracting " + img->getName());
        resultImg->subtract(img);
    }
    emit done(-1);

    display(resultImg);
}

void milxQtMain::imagesMultiply()
{
  if (getNumberOfWindows() == 0 || imageWindows.size() < 1)
    return;

  emit working(-1);
  initialiseWindowTraversal();
  QPointer<milxQtImage> firstImg = nextImage();

  QPointer<milxQtImage> resultImg = new milxQtImage;
  printInfo("Assigning " + firstImg->getName());
  resultImg->setName("Product Image");
  resultImg->setConsole(console);
  resultImg->setData(firstImg, true);
  for (int j = 1; j < imageWindows.size(); j++) //!< For all windows, do operation
  {
    QPointer<milxQtImage> img = nextImage();

    printInfo("Multiplying " + img->getName());
    resultImg->multiply(img);
  }
  emit done(-1);

  display(resultImg);
}

void milxQtMain::imagesConvolve()
{
    if(getNumberOfWindows() == 0 || imageWindows.size() < 1)
        return;
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
    emit working(-1);
    initialiseWindowTraversal();
    QPointer<milxQtImage> firstImg = nextImage();

    QPointer<milxQtImage> resultImg = new milxQtImage;
    printInfo("Assigning " + firstImg->getName());
    resultImg->setName("Convolution");
    resultImg->setConsole(console);
    resultImg->setData(firstImg, true);
    for(int j = 1; j < imageWindows.size(); j ++) //!< For all windows, do operation
    {
        QPointer<milxQtImage> img = nextImage();

        printInfo("Convolving with " + img->getName());
        resultImg->convolve(img);
    }
    emit done(-1);

    display(resultImg);
#endif
}

void milxQtMain::imagesMergeLabels()
{
    if(getNumberOfWindows() == 0 || imageWindows.size() < 1)
        return;
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
    emit working(-1);
    initialiseWindowTraversal();
    std::vector< charImageType::Pointer > images;
    std::vector<unsigned char> values;
    for(int j = 0; j < imageWindows.size(); j ++) //!< For all windows, do operation
    {
        QPointer<milxQtImage> img = nextImage();

        if(!img->is8BitImage() && !img->isFloatingPointImage())
        {
            printWarning("Ignoring non-labelled or float image " + img->getName());
            continue;
        }
        else if(img->isFloatingPointImage())
        {
            printWarning("Found float image " + img->getName());
            printWarning("Casting to 8-bit image for merge");
            images.push_back( milx::Image<floatImageType>::CastImage<charImageType>(img->GetFloatImage()) );
        }
        else
        {
            images.push_back(img->GetCharImage());
            printInfo("Merging with " + img->getName());
        }

        values.push_back(j);
    }

    if(images.empty())
    {
        printError("No labelled image found. Ignoring operation.");
        return;
    }

    charImageType::Pointer mergedImg = milx::Image<charImageType>::MergeLabelledImages(images);

    QPointer<milxQtImage> resultImg = new milxQtImage;
        resultImg->setName("Merged Labels");
        resultImg->setData(mergedImg, true);
        resultImg->setConsole(console);
        resultImg->generateImage();
    emit done(-1);

    display(resultImg);
#endif
}

void milxQtMain::imagesAverage()
{
    if(getNumberOfWindows() == 0 || imageWindows.size() < 1)
        return;

    imagesAdd();
    activeImage()->setName("Average Image");
    activeImage()->scale( 1.0/(getNumberOfWindows()-1) );
}

void milxQtMain::unify()
{
    ///Up cast and unify appropriately
    printInfo("Updating Multi-display");
    if(isActiveImage())
    {
        printInfo("Adding Image ...");
        addToUnifiedWindow(activeImage());
    }
    else if(isActiveModel())
    {
        printInfo("Adding Model ...");
        addToUnifiedWindow(activeModel());
    }

    printDebug("Uni Displayed");
    display(currentUnifiedWindow);
}

void milxQtMain::update()
{
    if(!imageWindows.isEmpty())
        menuImages->setDisabled(false);
    else
        menuImages->setDisabled(true);

    if(!imageWindows.isEmpty() || !modelWindows.isEmpty())
        menuData->setDisabled(false);
    else
        menuData->setDisabled(true);
}

void milxQtMain::linkProgressEventOf(vtkObject * obj)
{
    Connector->Connect(obj,
                       vtkCommand::ProgressEvent,
                       this,
                       SLOT( updateQtEvents() ),
                       NULL, 1.0); //High Priority
}

void milxQtMain::commonChildProperties(QWidget *widget)
{
    milxQtRenderWindow *customWindow = qobject_cast<milxQtRenderWindow *>(widget);

    if(customWindow == 0) //Not milxQt display
        widget->setAttribute(Qt::WA_DeleteOnClose);
    else
    {
        if(customWindow->isDeletableOnClose()) ///Check if deletable on close
            customWindow->setAttribute(Qt::WA_DeleteOnClose);
        customWindow->setConsole(console);
    }

    widget->setFocusPolicy(Qt::StrongFocus);

    customWindow->GetRenderWindow()->SetSize(subWindowSize, subWindowSize);
    customWindow->resize(subWindowSize, subWindowSize);
}

void milxQtMain::createMenu()
{
    menuBar = new QMenuBar(this);
    //File
    menuFile = new QMenu(menuBar);
    actionOpen = new QAction(this);
    actionOpenSeries = new QAction(this);
    actionOpenCollect = new QAction(this);
    actionSave = new QAction(this);
    actionSaveScreen = new QAction(this);
    actionCloseActive = new QAction(this);
    actionCloseAll = new QAction(this);
    for (int i = 0; i < MaxRecentFiles; ++i)
    {
        actionsRecentFile[i] = new QAction(this);
        actionsRecentFile[i]->setVisible(false);
        connect(actionsRecentFile[i], SIGNAL(triggered()), this, SLOT(openRecentFile()));
    }
    actionExit = new QAction(this);
    //New
    actionNewTab = new QAction(this);
    //Data
    menuData = new QMenu(menuBar);
    //Images
    menuImages = new QMenu(menuBar);
    actionBlendImages = new QAction(this);
    actionAddImages = new QAction(this);
    actionAverageImages = new QAction(this);
    actionSubtractImages = new QAction(this);
    actionMultiplyImages = new QAction(this);
    actionConvolveImages = new QAction(this);
    actionMergeLabels = new QAction(this);
    //Window
    menuWindows = new QMenu(menuBar);
    actionLinkWindows = new QAction(this);
    actionCascade = new QAction(this);
    actionTile = new QAction(this);
    actionTileVertically = new QAction(this);
    actionTileHorizontally = new QAction(this);
    menuWindowList = new QMenu(menuBar);
    importFromMenu = new QMenu(this);
    //Help
    menuHelp = new QMenu(menuBar);
    actionPreferences = new QAction(this);
    actionContents = new QAction(this);
    actionControls = new QAction(this);
    actionAbout = new QAction(this);

    ///Setup Exit Action and File Menus
    ///File
    menuBar->addAction(menuFile->menuAction());
    menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
    actionNewTab->setIcon(QIcon(":/resources/toolbar/new_tab.png"));
    actionNewTab->setText(QApplication::translate("MainWindow", "New Tab", 0, QApplication::UnicodeUTF8));
    actionNewTab->setShortcut(tr("Ctrl+t"));
    menuFile->addAction(actionNewTab);
    actionOpen->setIcon(QIcon(":/resources/toolbar/open.png"));
    actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
    actionOpen->setShortcut(tr("Ctrl+o"));
    menuFile->addAction(actionOpen);
    actionOpenSeries->setIcon(QIcon(":/resources/toolbar/open_series.png"));
    actionOpenSeries->setText(QApplication::translate("MainWindow", "Open DICOM Series", 0, QApplication::UnicodeUTF8));
    actionOpenSeries->setShortcut(tr("Ctrl+Alt+o"));
    menuFile->addAction(actionOpenSeries);
    actionOpenCollect->setIcon(QIcon(":/resources/toolbar/open_collection.png"));
    actionOpenCollect->setText(QApplication::translate("MainWindow", "Open Model Collection", 0, QApplication::UnicodeUTF8));
    actionOpenCollect->setShortcut(tr("Ctrl+Shift+o"));
    menuFile->addAction(actionOpenCollect);
    actionSave->setIcon(QIcon(":/resources/toolbar/save.png"));
    actionSave->setText(QApplication::translate("MainWindow", "Save", 0, QApplication::UnicodeUTF8));
    actionSave->setShortcut(tr("Ctrl+s"));
    menuFile->addAction(actionSave);
    actionSaveScreen->setIcon(QIcon(":/resources/toolbar/screenshot.png"));
    actionSaveScreen->setText(QApplication::translate("MainWindow", "Save Screenshot", 0, QApplication::UnicodeUTF8));
    actionSaveScreen->setShortcut(tr("Shift+Ctrl+s"));
    menuFile->addAction(actionSaveScreen);
    menuFile->addSeparator();
    //Closing
    menuFile->addAction(actionCloseActive);
    actionCloseActive->setIcon(QIcon(":/resources/toolbar/close.png"));
    actionCloseActive->setText(QApplication::translate("MainWindow", "Close", 0, QApplication::UnicodeUTF8));
    actionCloseActive->setShortcut(tr("Ctrl+w"));
    menuFile->addAction(actionCloseActive);
    actionCloseAll->setIcon(QIcon(":/resources/toolbar/close_all.png"));
    actionCloseAll->setText(QApplication::translate("MainWindow", "Close All", 0, QApplication::UnicodeUTF8));
    actionCloseAll->setShortcut(tr("Ctrl+Shift+w"));
    menuFile->addAction(actionCloseAll);
    menuFile->addSeparator();
    //Recent file list
    for (int i = 0; i < MaxRecentFiles; ++i)
        menuFile->addAction(actionsRecentFile[i]);
    actionRecentFileSeparator = menuFile->addSeparator();
    updateRecentFileActions(); //Do not move up
    //Exit
    menuFile->addAction(actionPreferences);
    actionPreferences->setIcon(QIcon(":/resources/toolbar/settings.png"));
    actionPreferences->setText(QApplication::translate("MainWindow", "Preferences ...", 0, QApplication::UnicodeUTF8));
    actionPreferences->setShortcut(tr("Ctrl+p"));
    actionExit->setIcon(QIcon(":/resources/toolbar/exit.png"));
    actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
    actionExit->setShortcut(tr("Ctrl+x"));
    menuFile->addAction(actionExit);
    ///Data
    menuBar->addAction(menuData->menuAction());
    menuData->setTitle(QApplication::translate("MainWindow", "Data", 0, QApplication::UnicodeUTF8));
    ///Images
    menuImages->setTitle(QApplication::translate("MainWindow", "Images", 0, QApplication::UnicodeUTF8));
    menuBar->addAction(menuImages->menuAction());
    actionBlendImages->setText(QApplication::translate("MainWindow", "Blend Images", 0, QApplication::UnicodeUTF8));
    actionBlendImages->setShortcut(tr("Ctrl+b"));
    menuImages->addAction(actionBlendImages);
    actionAddImages->setText(QApplication::translate("MainWindow", "Add Images", 0, QApplication::UnicodeUTF8));
    actionAddImages->setShortcut(tr("Ctrl+a"));
    menuImages->addAction(actionAddImages);
    actionAverageImages->setText(QApplication::translate("MainWindow", "Average Images", 0, QApplication::UnicodeUTF8));
    actionAverageImages->setShortcut(tr("Ctrl+Shift+a"));
    menuImages->addAction(actionAverageImages);
    actionSubtractImages->setText(QApplication::translate("MainWindow", "Difference Images", 0, QApplication::UnicodeUTF8));
    actionSubtractImages->setShortcut(tr("Ctrl+d"));
    menuImages->addAction(actionSubtractImages);
    actionMultiplyImages->setText(QApplication::translate("MainWindow", "Multiply Images", 0, QApplication::UnicodeUTF8));
    actionMultiplyImages->setShortcut(tr("Ctrl+d"));
    menuImages->addAction(actionMultiplyImages);
    actionConvolveImages->setText(QApplication::translate("MainWindow", "Convolve Images", 0, QApplication::UnicodeUTF8));
    actionConvolveImages->setShortcut(tr("Ctrl+c"));
    menuImages->addAction(actionConvolveImages);
    actionMergeLabels->setText(QApplication::translate("MainWindow", "Merge Labels", 0, QApplication::UnicodeUTF8));
    actionMergeLabels->setShortcut(tr("Ctrl+l"));
    menuImages->addAction(actionMergeLabels);
    ///Windows
    menuBar->addAction(menuWindows->menuAction());
    menuWindows->setTitle(QApplication::translate("MainWindow", "Windows", 0, QApplication::UnicodeUTF8));
    actionLinkWindows->setText(QApplication::translate("MainWindow", "Link All Windows", 0, QApplication::UnicodeUTF8));
    actionLinkWindows->setIcon(QIcon(":/resources/toolbar/link.png"));
    actionLinkWindows->setCheckable(true);
    menuWindows->addAction(actionLinkWindows);
    actionCascade->setIcon(QIcon(":/resources/toolbar/cascade.png"));
    actionCascade->setText(QApplication::translate("MainWindow", "Cascade", 0, QApplication::UnicodeUTF8));
    menuWindows->addAction(actionCascade);
    actionTile->setIcon(QIcon(":/resources/toolbar/tile.png"));
    actionTile->setText(QApplication::translate("MainWindow", "Tile", 0, QApplication::UnicodeUTF8));
    menuWindows->addAction(actionTile);
    actionTileVertically->setIcon(QIcon(":/resources/toolbar/tilev.png"));
    actionTileVertically->setText(QApplication::translate("MainWindow", "Tile Vertically", 0, QApplication::UnicodeUTF8));
    menuWindows->addAction(actionTileVertically);
    actionTileHorizontally->setIcon(QIcon(":/resources/toolbar/tileh.png"));
    actionTileHorizontally->setText(QApplication::translate("MainWindow", "Tile Horizontally", 0, QApplication::UnicodeUTF8));
    menuWindows->addAction(actionTileHorizontally);
    updateWindowMenu();
    connect(menuWindows, SIGNAL(aboutToShow()), this, SLOT(updateWindowMenu()));
    ///Help
    menuHelp->setTitle(QApplication::translate("MainWindow", "Help", 0, QApplication::UnicodeUTF8));
    menuBar->addAction(menuHelp->menuAction());
    actionContents->setText(QApplication::translate("MainWindow", "sMILX Help Contents", 0, QApplication::UnicodeUTF8));
    actionContents->setIcon(QIcon(":/resources/toolbar/help.png"));
    actionContents->setShortcut(tr("F1"));
    menuHelp->addAction(actionContents);
    actionControls->setText(QApplication::translate("MainWindow", "Controls", 0, QApplication::UnicodeUTF8));
    actionControls->setShortcut(tr("F2"));
    actionControls->setIcon(QIcon(":/resources/toolbar/controls.png"));
    menuHelp->addAction(actionControls);
    actionAbout->setText(QApplication::translate("MainWindow", "About", 0, QApplication::UnicodeUTF8));
    actionAbout->setShortcut(tr("Ctrl+h"));
    menuHelp->addAction(actionAbout);

    ///Common Actions/Menus
    actionCompare = new QAction(this);
    actionCompare->setText(QApplication::translate("MainWindow", "Compare", 0, QApplication::UnicodeUTF8));
    actionCompare->setShortcut(tr("Ctrl+u"));
    actionCompare->setDisabled(true); ///\todo enable compare when fixed
    menuWindowList->setTitle(QApplication::translate("MainWindow", "Switch Window To", 0, QApplication::UnicodeUTF8));
    connect(menuWindowList, SIGNAL(aboutToShow()), this, SLOT(updateWindowListMenu()));
    importFromMenu->setTitle(QApplication::translate("MainWindow", "Import View From", 0, QApplication::UnicodeUTF8));
    connect(importFromMenu, SIGNAL(aboutToShow()), this, SLOT(updateImportFromMenu()));

#if !(ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
    actionConvolveImages->setDisabled(true);
    actionMergeLabels->setDisabled(true);
#endif

    ///Image Toolbar actions
    actionImageText = new QAction(this);
    actionImageText->setText(QApplication::translate("MainWindow", "Text", 0, QApplication::UnicodeUTF8));
    actionImageText->setShortcut(tr("Ctrl+Alt+t"));

    statusBar()->showMessage(tr("Ready"));
    setMenuBar(menuBar);
}

void milxQtMain::createComboBoxes()
{
    defaultViewBox = new QComboBox(this);
    defaultViewBox->insertItem(AXIAL, "View: Axial");
    defaultViewBox->insertItem(CORONAL, "View: Coronal");
    defaultViewBox->insertItem(SAGITTAL, "View: Sagittal"); ///Assigned is deliberately reversed, \todo Why need to reverse?
    defaultViewBox->setEditable(false);
    QObject::connect(defaultViewBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateWindowsWithView(int)));

    defaultViewTypeBox = new QComboBox(this);
    defaultViewTypeBox->insertItem(SINGLE, "View Type: Single");
    defaultViewTypeBox->insertItem(SCANNER, "View Type: Multiple");
    defaultViewTypeBox->setEditable(false);
    QObject::connect(defaultViewTypeBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateWindowsWithViewType(int)));

    defaultOrientationTypeBox = new QComboBox(this);
    defaultOrientationTypeBox->insertItem(RADIOLOGICAL, "Orient. Type: Radiological");
    defaultOrientationTypeBox->insertItem(NEUROLOGICAL, "Orient. Type: Neurological");
    defaultOrientationTypeBox->setEditable(false);
    QObject::connect(defaultOrientationTypeBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateWindowsWithViewOrientation(int)));
}

void milxQtMain::createToolBars()
{
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(actionNewTab);
    fileToolBar->addAction(actionOpen);
    fileToolBar->addAction(actionOpenSeries);
    fileToolBar->addAction(actionOpenCollect);
    fileToolBar->addAction(actionSave);
    fileToolBar->addAction(actionSaveScreen);
    fileToolBar->addAction(actionCloseActive);
    fileToolBar->addAction(actionCloseAll);
    fileToolBar->setObjectName("FileToolBar");

//    editToolBar = addToolBar(tr("Edit"));
//    editToolBar->addAction(cutAct);
//    editToolBar->addAction(copyAct);
//    editToolBar->addAction(pasteAct);

    windowToolBar = addToolBar(tr("Window"));
    windowToolBar->addAction(actionTile);
    windowToolBar->addAction(actionTileVertically);
    windowToolBar->addAction(actionTileHorizontally);
    windowToolBar->addAction(actionCascade);
    windowToolBar->addAction(actionLinkWindows);
    windowToolBar->addAction(actionConsole);
    windowToolBar->addAction(actionContents);
    windowToolBar->addAction(actionControls);
    windowToolBar->addAction(actionPreferences);
    windowToolBar->setObjectName("WindowToolBar");

    defaultToolBar = addToolBar(tr("Default"));
    QLabel *defaultBarLabel = new QLabel(this);
//    defaultBarLabel->setText("View: ");
    defaultToolBar->addWidget(defaultBarLabel);
    defaultToolBar->addWidget(defaultViewBox);
    defaultToolBar->addWidget(defaultViewTypeBox);
    defaultToolBar->addWidget(defaultOrientationTypeBox);
    defaultToolBar->setObjectName("DefaultToolBar");

    //Sliders etc.
    imageToolBar = new QToolBar(tr("Image"), this);
    imageToolBar->setObjectName("Image");
    imageLevelSlider = new QSlider(Qt::Vertical, this);
    imageLevelSlider->setMinimum(0);
    imageLevelSlider->setMaximum(100);
    imageLevelSlider->setSingleStep(1);
    imageLevelSlider->setValue(50);
    imageLevelSlider->setTickPosition(QSlider::TicksRight);
    imageLevelSlider->setTickInterval(1);
    QObject::connect(imageLevelSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWindowsWithValue(int)));
    imageLevelButton = new QPushButton(tr(""), this);
    imageLevelButton->setIcon(QIcon(":/resources/toolbar/intensity.png"));
    QObject::connect(imageLevelButton, SIGNAL(clicked()), this, SLOT(updateWindowsWithAutoLevel()));
    refreshButton = new QPushButton(tr(""), this);
    refreshButton->setIcon(QIcon(":/resources/toolbar/refresh.png"));
    QObject::connect(refreshButton, SIGNAL(clicked()), this, SLOT(updateWindowsWithRefresh()));
    cursorButton = new QPushButton(tr(""), this);
    cursorButton->setIcon(QIcon(":/resources/toolbar/crosshairs_2D.png"));
    QObject::connect(cursorButton, SIGNAL(clicked()), this, SLOT(updateWindowsWithCursors()));
    /*imageLevelDial = new QDial(this);
//    imageLevelDial->setFloatable(true);
    imageLevelDial->setMinimum(0);
    imageLevelDial->setMaximum(100);
    imageLevelDial->setSingleStep(1);
    imageLevelDial->setValue(50);*/
    addToolBar(Qt::LeftToolBarArea, imageToolBar);
//    imageToolBar->addAction(actionImageText);
    imageToolBar->addWidget(refreshButton);
    imageToolBar->addWidget(imageLevelButton);
    imageToolBar->addWidget(cursorButton);
    imageToolBar->addWidget(imageLevelSlider);
//    imageToolBar->addWidget(imageLevelDial);
}

void milxQtMain::createConnections()
{
    qRegisterMetaType< milxQtRenderWindow* >("milxQtRenderWindow*"); ///Need for ensuring cross thread signalling
    qRegisterMetaType< milxQtModel* >("milxQtModel*");
    qRegisterMetaType< milxQtImage* >("milxQtImage*");

    QObject::connect(workspaces, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTab(int)));
    QObject::connect(windowMapper, SIGNAL(mapped(QWidget *)), this, SLOT(setActiveWindow(QWidget *)));
    ///Actions
    ///File
    QObject::connect(actionNewTab, SIGNAL(activated()), this, SLOT(newTab()));
    QObject::connect(actionOpen, SIGNAL(activated()), this, SLOT(open()));
    QObject::connect(actionOpenSeries, SIGNAL(activated()), this, SLOT(openSeries()));
    QObject::connect(actionOpenCollect, SIGNAL(activated()), this, SLOT(openCollection()));
    QObject::connect(actionSave, SIGNAL(activated()), this, SLOT(save()));
    QObject::connect(actionSaveScreen, SIGNAL(activated()), this, SLOT(saveScreen()));
    QObject::connect(actionCloseActive, SIGNAL(activated()), this, SLOT(closeTabActiveWindow()));
    QObject::connect(actionCloseAll, SIGNAL(activated()), this, SLOT(closeTabAllWindows()));
    QObject::connect(actionExit, SIGNAL(activated()), this, SLOT(close()));
    ///Data
    QObject::connect(menuData, SIGNAL(aboutToShow()), this, SLOT(dataMenu()));
    ///Images
    QObject::connect(actionBlendImages, SIGNAL(activated()), this, SLOT(imagesMix()));
    QObject::connect(actionAddImages, SIGNAL(activated()), this, SLOT(imagesAdd()));
    QObject::connect(actionAverageImages, SIGNAL(activated()), this, SLOT(imagesAverage()));
    QObject::connect(actionSubtractImages, SIGNAL(activated()), this, SLOT(imagesSubtract()));
    QObject::connect(actionMultiplyImages, SIGNAL(activated()), this, SLOT(imagesMultiply()));
  #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
    QObject::connect(actionConvolveImages, SIGNAL(activated()), this, SLOT(imagesConvolve()));
  #endif
    QObject::connect(actionMergeLabels, SIGNAL(activated()), this, SLOT(imagesMergeLabels()));
    ///Windows
    //actionLinkWindows is not connected because its used as a Boolean in the transferViewToWindows() member
    QObject::connect(actionCascade, SIGNAL(activated()), this, SLOT(cascadeTab()));
    QObject::connect(actionTile, SIGNAL(activated()), this, SLOT(tileTab()));
    QObject::connect(actionTileVertically, SIGNAL(activated()), this, SLOT(tileTabVertically()));
    QObject::connect(actionTileHorizontally, SIGNAL(activated()), this, SLOT(tileTabHorizontally()));
    ///Help
    QObject::connect(actionContents, SIGNAL(activated()), this, SLOT(helpContents()));
    QObject::connect(actionPreferences, SIGNAL(activated()), this, SLOT(preferences()));
    QObject::connect(actionControls, SIGNAL(activated()), this, SLOT(controls()));
    QObject::connect(actionAbout, SIGNAL(activated()), this, SLOT(about()));

    ///Common actions
    QObject::connect(actionCompare, SIGNAL(activated()), this, SLOT(unify()));
}

void milxQtMain::createProgressBar()
{
    progressCallCount = 1; //done decrements
    progressBar = new QProgressBar(this); //hierarchical deletion
    progressBar->setMaximumSize(128, progressBar->maximumSize().height());
    statusBar()->addPermanentWidget(progressBar);
    done(-1);
}

void milxQtMain::setupTooltips()
{
    actionNewTab->setToolTip("Create a new tab for display");
    actionNewTab->setStatusTip("Create a new tab for display");
    actionOpen->setToolTip("Open a file (all supported formats, including those by plugins) for display");
    actionOpen->setStatusTip("Open a file (all supported formats, including those by plugins) for display");
    actionOpenSeries->setToolTip("Open a DICOM series by choosing a directory");
    actionOpenSeries->setStatusTip("Open a DICOM series by choosing a directory");
    actionOpenCollect->setToolTip("Open a model collection by choosing a set of files");
    actionOpenCollect->setStatusTip("Open a model collection by choosing a set of files");
    actionSave->setToolTip("Saves the data in the active window (as a supported file format, including those in plugins)");
    actionSave->setStatusTip("Saves the data in the active window (as a supported file format, including those in plugins)");
    actionSaveScreen->setToolTip("Save a screenshot of the active window exactly as displayed");
    actionSaveScreen->setStatusTip("Save a screenshot of the active window exactly as displayed");
    actionCloseActive->setToolTip("Close active window");
    actionCloseActive->setStatusTip("Close active window");
    actionCloseAll->setToolTip("Close all windows in the tab");
    actionCloseAll->setStatusTip("Close all windows in the tab");
    actionCascade->setToolTip("Cascade windows in current tab");
    actionCascade->setStatusTip("Cascade windows in current tab");
    actionTile->setToolTip("Tile windows in current tab");
    actionTile->setStatusTip("Tile windows in current tab");
    actionTileVertically->setToolTip("Tile windows vertically in current tab");
    actionTileVertically->setStatusTip("Tile windows vertically in current tab");
    actionTileHorizontally->setToolTip("Tile windows horizontally in current tab");
    actionTileHorizontally->setStatusTip("Tile windows horizontally in current tab");
    actionConsole->setToolTip("Hide/Show console docking window");
    actionConsole->setStatusTip("Hide/Show console docking window");
    actionLinkWindows->setToolTip("Link the view in all windows");
    actionLinkWindows->setStatusTip("Link the view in all windows");
    actionContents->setToolTip("Browse the help");
    actionContents->setStatusTip("Browse the help");
    actionPreferences->setToolTip("Customise the program settings");
    actionPreferences->setStatusTip("Customise the program settings");

    actionCompare->setToolTip("Compare the data in window with others by placing it into a multi-display window");
    actionCompare->setStatusTip("Compare the data in window with others by placing it into a multi-display window");

    refreshButton->setToolTip("Refresh the window views, levels and camera to default.");
    refreshButton->setStatusTip("Refresh the window views, levels and camera to default.");
    imageLevelButton->setToolTip("Auto level the image intensities based on inter-quartile ranges.");
    imageLevelButton->setStatusTip("Auto level the image intensities based on inter-quartile ranges.");
    cursorButton->setToolTip("Enable crosshairs for all images.");
    cursorButton->setStatusTip("Enable crosshairs for all images.");
}

void milxQtMain::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    QMenu* contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection

    contextMenu->addAction(actionNewTab);
    contextMenu->addAction(actionOpen);
    contextMenu->addAction(actionOpenSeries);
    contextMenu->addAction(actionOpenCollect);
    contextMenu->addAction(actionSave);
    contextMenu->addAction(actionSaveScreen);
    contextMenu->addSeparator();
    contextMenu->addAction(actionPreferences);
    contextMenu->addAction(actionCloseActive);
    contextMenu->addAction(actionCloseAll);
    contextMenu->addSeparator();
    contextMenu->addAction(actionCascade);
    contextMenu->addAction(actionTile);
    contextMenu->addAction(actionTileVertically);
    contextMenu->addAction(actionTileHorizontally);
    contextMenu->addSeparator();
    contextMenu->addAction(actionExit);

    contextMenu->exec(currentEvent->globalPos());
}

void milxQtMain::dragEnterEvent(QDragEnterEvent *currentEvent)
{
    if(currentEvent->mimeData()->hasFormat("text/uri-list") || currentEvent->mimeData()->hasFormat("text/plain"))
        currentEvent->acceptProposedAction();
}

void milxQtMain::dropEvent(QDropEvent *currentEvent)
{
    QList<QUrl> urlsList = currentEvent->mimeData()->urls();
    QString tmp;

    for(int j = 0; j < urlsList.size(); j ++)
    {
        if(urlsList[j].isValid())
        {
#ifdef Q_WS_WIN
            tmp = urlsList[j].path().remove(0,1); //!< Remove leading forward slash
#else
            tmp = urlsList[j].path();
#endif
            printInfo("Dropped Path: " + tmp);
            loadFile(tmp);
        }
    }

    currentEvent->acceptProposedAction();
}

void milxQtMain::closeEvent(QCloseEvent *event)
{
	writeSettings();
    event->accept();
}

bool milxQtMain::loadPlugins()
{
    QDir pluginsDir(qApp->applicationDirPath());
    QString sharedObjectSuffix = "";
#if defined(Q_OS_WIN)
    if (pluginsDir.dirName().toLower() == "debug" || pluginsDir.dirName().toLower() == "release")
        pluginsDir.cdUp();
    sharedObjectSuffix = "dll"; //!<\todo use better cross-platform method here
#elif defined(Q_OS_MAC)
    if (pluginsDir.dirName() == "MacOS")
    {
        pluginsDir.cdUp();
        pluginsDir.cdUp();
        pluginsDir.cdUp();
    }
    sharedObjectSuffix = "dylib"; //!<\todo use better cross-platform method here
#else
    sharedObjectSuffix = "so"; //!<\todo use better cross-platform method here
#endif
    pluginsDir.cdUp();
    pluginsDir.cd("plugins");
    printInfo("Plugin Path to be searched: " + pluginsDir.path());
    
    QStringList pluginFileNames = pluginsDir.entryList(QDir::Files);
    foreach (QString fileName, pluginFileNames)
    {
        QFileInfo fi(fileName);

    #if defined(Q_OS_WIN)
        if(!fi.completeSuffix().contains(sharedObjectSuffix, Qt::CaseInsensitive))
            continue;
    #else
//        printDebug("File: " + fileName + ", ext: " + fi.completeSuffix());
        if(!fi.completeSuffix().contains(sharedObjectSuffix, Qt::CaseSensitive))
            continue;
    #endif

        printDebug("Attempting to Load " + fileName);
        QPluginLoader pluginLoader(pluginsDir.absoluteFilePath(fileName));
        QObject *plugin = pluginLoader.instance();
        if (plugin)
        {
            printDebug("Instanced ... ");
            milxQtPluginFactory *pluginFactory = qobject_cast<milxQtPluginFactory *>(plugin);
            if (pluginFactory)
            {
                QPointer<milxQtPluginInterface> loadedPlugin = pluginFactory->newPlugin(this); //hierarchical deletion
                plugins.append( loadedPlugin );

                connect(loadedPlugin, SIGNAL(resultAvailable(milxQtRenderWindow*)), this, SLOT(display(milxQtRenderWindow*)));
                connect(loadedPlugin, SIGNAL(resultAvailable(milxQtModel*)), this, SLOT(display(milxQtModel*)));
                connect(loadedPlugin, SIGNAL(resultAvailable(milxQtImage*)), this, SLOT(display(milxQtImage*)));
                connect(loadedPlugin, SIGNAL( resultAvailable(vtkPolyDataCollection*, QStringList&) ), this, SLOT( display(vtkPolyDataCollection*, QStringList&) ));

                if(loadedPlugin->isThreaded())
                {
                    ///Connect cross-thread signals and slots
                    connect(loadedPlugin, SIGNAL(finished()), loadedPlugin, SLOT(postStartTasks()));
                }

                if(loadedPlugin->isDockable()) //!< If a dock plugin, load the dock widget
                {
                    if(!loadedPlugin->isConsole())
                    {
                        addDockWidget(loadedPlugin->dockDefaultArea(), loadedPlugin->dockWidget());
                        dockActions.append(loadedPlugin->dockWidget()->toggleViewAction());
    //                    loadedPlugin->start(); //!< Start the plugins event loop
                    }
                    else
                    {
                        console->setTab(loadedPlugin->dockWidget()->widget());
                    }
                }

                if(loadedPlugin->isExtension()) //!< If an extension, load it
                {
                    QString tmpName = "[" + loadedPlugin->name() + " Extension]";

                    QAction *extAction = new QAction(this);
                    extAction->setText(QApplication::translate("Plugin", tmpName.toStdString().c_str(), 0, QApplication::UnicodeUTF8));
                    extAction->setShortcut(tr("Alt+n"));

                    if(loadedPlugin->genericResult())
                        renderExtsActions.append(extAction);

                    if(loadedPlugin->modelResult())
                        modelExtsActions.append(extAction);

                    if(loadedPlugin->imageResult())
                        imageExtsActions.append(extAction);

                    connect(extAction, SIGNAL(triggered()), loadedPlugin, SLOT(loadExtension()));
                }

                foreach(QMenu *menuToAdd, loadedPlugin->addToMenuBar())
                {
                    menuBar->insertAction(menuWindows->menuAction(), menuToAdd->menuAction());
                }

                loadedPlugin->setConsole(console);

                //Update support formats strings
                foreach (QString ext, loadedPlugin->saveExtensions())
                {
                    saveSupport += " *" + ext;
                }
                foreach (QString ext, loadedPlugin->openExtensions())
                {
                    openSupport += " *" + ext;
                }

                printInfo(plugins.last()->name() + " plugin was successfully loaded.");
            }
            //printDebug("Unloading Plugin Loader.");
            //pluginLoader.unload();
        }
        else
            printError("Failed to load plugin: " + pluginLoader.errorString());
    }
    printDebug("Support Save Formats: " + saveSupport);
    printDebug("Support Open Formats: " + openSupport);

    return false;
}

void milxQtMain::writeSettings()
{
    QSettings settings("Shekhar Chandra", "milxQt");

	if(resettingInterface)
		return;

    settings.beginGroup("milxQtMain");
    settings.setValue("size", size());
    settings.setValue("pos", pos());
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
    settings.setValue("defaultView", defaultViewBox->currentIndex());
    settings.setValue("defaultViewType", defaultViewTypeBox->currentIndex());
    settings.setValue("defaultOrientationType", defaultOrientationTypeBox->currentIndex());
    //Options
    settings.setValue("whiteBackground", whiteBackground);
    settings.setValue("humanGlyph", humanGlyph);
    settings.setValue("subWindowSize", subWindowSize);
    settings.setValue("maxProcessors", maxProcessors);
    settings.setValue("magnifyFactor", magnifyFactor);
    settings.setValue("timestamping", timestamping);
    settings.setValue("interpolationImages", interpolationImages);
    settings.setValue("orientationImages", orientationImages);
    settings.setValue("interpolationModels", interpolationModels);
    settings.setValue("scalarBarModels", scalarBarModels);
    settings.endGroup();

	//resettingInterface = false;
}

void milxQtMain::resetSettings()
{
	QSettings settings("Shekhar Chandra", "milxQt");

	//New defaults
	QSize desktopSize = qApp->desktop()->availableGeometry().size();
	int newWidth = 2.0*desktopSize.width() / 3.0 + 0.5;
	int newHeight = 4.0*desktopSize.height() / 5.0 + 0.5;
	int xOffset = (desktopSize.width() - newWidth) / 2.0;
	int yOffset = (desktopSize.height() - newHeight) / 2.0;
	int defaultViewMode = 2; //axial
	int defaultViewTypeMode = 0; //1-multi-view
	int defaultOrientationTypeMode = 0; //radiological

	settings.beginGroup("milxQtMain");
	settings.setValue("size", QSize(newWidth, newHeight));
	settings.setValue("pos", QPoint(xOffset, yOffset));
	settings.remove("geometry");
	settings.remove("windowState");
	settings.setValue("defaultView", defaultViewMode);
	settings.setValue("defaultViewType", defaultViewTypeMode);
	settings.setValue("defaultOrientationType", defaultOrientationTypeMode);
	settings.endGroup();

	printInfo("Toolbars, window positions etc. have been reset.");
	resettingInterface = true;

	QMessageBox msgBox;
	msgBox.setText("Need to restart to take effect");
	msgBox.setInformativeText("Other changes to the preferences have been ignored.");
	msgBox.setStandardButtons(QMessageBox::Ok);
	int ret = msgBox.exec();
}

void milxQtMain::readSettings()
{
    QSettings settings("Shekhar Chandra", "milxQt");

    //New defaults
    QSize desktopSize = qApp->desktop()->availableGeometry().size();
    int newWidth = 2.0*desktopSize.width()/3.0 + 0.5;
    int newHeight = 4.0*desktopSize.height()/5.0 + 0.5;
    int xOffset = (desktopSize.width()-newWidth)/2.0;
    int yOffset = (desktopSize.height()-newHeight)/2.0;
    int defaultViewMode = 2; //axial
    int defaultViewTypeMode = 0; //1-multi-view
    int defaultOrientationTypeMode = 0; //radiological

    settings.beginGroup("milxQtMain");
    resize( settings.value("size", QSize(newWidth, newHeight)).toSize() );
    move( settings.value("pos", QPoint(xOffset, yOffset)).toPoint() );
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
    defaultViewBox->setCurrentIndex( settings.value("defaultView", defaultViewMode).toInt() );
    defaultViewTypeBox->setCurrentIndex( settings.value("defaultViewType", defaultViewTypeMode).toInt() );
    defaultOrientationTypeBox->setCurrentIndex( settings.value("defaultOrientationType", defaultOrientationTypeMode).toInt() );
    //Options
    whiteBackground = settings.value("whiteBackground", whiteBackground).toBool();
    humanGlyph = settings.value("humanGlyph", humanGlyph).toBool();
    subWindowSize = settings.value("subWindowSize", subWindowSize).toInt();
    maxProcessors = settings.value("maxProcessors", maxProcessors).toInt();
    magnifyFactor = settings.value("magnifyFactor", magnifyFactor).toInt();
    timestamping = settings.value("timestamping", timestamping).toBool();
    interpolationImages = settings.value("interpolationImages", interpolationImages).toBool();
    orientationImages = settings.value("orientationImages", orientationImages).toBool();
    interpolationModels = settings.value("interpolationModels", interpolationModels).toBool();
    scalarBarModels = settings.value("scalarBarModels", scalarBarModels).toBool();

    ///Handle saving dock positions/areas etc.
    restoreDockWidget(console->dockWidget());
    console->setTimestamps(timestamping);
	resettingInterface = false;

    settings.endGroup();
}

