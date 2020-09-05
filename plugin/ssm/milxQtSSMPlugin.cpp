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
#include "milxQtSSMPlugin.h"

#include <qplugin.h>

//Qt
#include <QMenu>
#include <QFileDialog>
#include <QGroupBox>
#include <QFormLayout>
#include <QWizardPage>
#include <QLineEdit>
#include <QListWidget>
#include <QWizard>

//SMILI
#include <milxQtFile.h>

milxQtSSMPlugin::milxQtSSMPlugin(QObject *theParent) : milxQtPluginInterface(theParent)
{
    ///Up cast parent to milxQtMain
    MainWindow = qobject_cast<milxQtMain *>(theParent);

    threaded = false;
    dockable = true;
    consoleWindow = false;
    extension = false;
    pluginName = "Statistical Shape Model (SSM)";
    dataName = "";
    modelManagerCreated = false;
    caseManagerCreated = false;

    modelTabIndex = 0;
    caseTabIndex = 0;

    manager = new milxQtManager(MainWindow);
        manager->hide();
    dock = new QDockWidget(tr("Manager"), MainWindow);
        dock->setFeatures(QDockWidget::AllDockWidgetFeatures);
        dock->setWidget(manager);
        dock->setObjectName("Manager");

    //focus variables
    txtAtlasName = new QLineEdit;
    btnAtlasName = new QPushButton;
    comboSurfaceNames = new QListWidget;
    btnSurfaceNames = new QPushButton;
    btnClearSurfaceNames = new QPushButton;

    createActions();
    createMenu();
    createWizard();
    createConnections();
}

milxQtSSMPlugin::~milxQtSSMPlugin()
{
    shapes.clear();
    robustShapes.clear();
    cerr << "SSM Destroyed" << endl;
}

QString milxQtSSMPlugin::name()
{
    return pluginName;
}

QString milxQtSSMPlugin::openFileSupport()
{
    QString openSSMExt = "Shape Model Files (*.ssm *.rssm)";

    return openSSMExt;
}

QStringList milxQtSSMPlugin::openExtensions()
{
    QStringList exts;

    exts.append(".ssm");
    exts.append(".rssm");

    return exts;
}

QStringList milxQtSSMPlugin::saveExtensions()
{
    QStringList exts;

    exts.append(".ssm");
    exts.append(".rssm");

    return exts;
}

QString milxQtSSMPlugin::saveFileSupport()
{
    QString saveSSMExt = "Shape Model Files (*.ssm *.rssm)";

    return saveSSMExt;
}

void milxQtSSMPlugin::SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames)
{
    cout << "Loading Collection via the SSM Plugin" << endl;
    addShapeModel( new milxQtShapeModel );
        currentModel = shapes.last();
        shapes.last()->setConsole(console);
        shapes.last()->SetInputCollection(collection, filenames);
        connect(shapes.last(), SIGNAL(resultAvailable(milxQtRenderWindow*)), MainWindow, SLOT(display(milxQtRenderWindow*)));
        connect(shapes.last(), SIGNAL(resultAvailable(milxQtModel*)), MainWindow, SLOT(display(milxQtModel*)));
        connect(shapes.last(), SIGNAL(collectionAvailable(vtkPolyDataCollection*, QStringList&)), this, SLOT(passOnCollection(vtkPolyDataCollection*, QStringList&)));
        cout << "Loaded Collection as a Normal SSM." << endl;

    ///Update case list in manager
    QStringList headingList;

    headingList << "Index" << "Case ID"; //!< Case Browser
    if(!caseManagerCreated)
    {
        caseTabIndex = manager->newTab("Case Browser", headingList);
        caseManagerCreated = true;
    }
    else
        manager->clearTab(caseTabIndex);

    cout << "Loading cases into browser." << endl;
    QList< int > cases = shapes.last()->getCaseIDs();
    for(int j = 0; j < cases.size(); j ++)
    {
        QStringList caseList;

        caseList << QString::number(j) << "Case " + QString::number(cases[j]);

        manager->addItem(caseTabIndex, caseList, Qt::NoItemFlags);
    }

    manager->show();
    cout << "Done." << endl;
}

void milxQtSSMPlugin::SetInputCollection(vtkPolyDataCollection* collection, vtkPolyData *atlasSurface, QStringList &filenames)
{
  cout << "Loading Collection via the SSM Plugin" << endl;
  addShapeModel( new milxQtRobustShapeModel );
      currentModel = robustShapes.last();
      robustShapes.last()->setConsole(console);
      robustShapes.last()->SetInputCollection(collection, atlasSurface, filenames);
      connect(robustShapes.last(), SIGNAL(resultAvailable(milxQtRenderWindow*)), MainWindow, SLOT(display(milxQtRenderWindow*)));
      connect(robustShapes.last(), SIGNAL(resultAvailable(milxQtModel*)), MainWindow, SLOT(display(milxQtModel*)));
      connect(robustShapes.last(), SIGNAL(collectionAvailable(vtkPolyDataCollection*, QStringList&)), this, SLOT(passOnCollection(vtkPolyDataCollection*, QStringList&)));
      cout << "Loaded Collection as a Robust SSM." << endl;

  if(!robustShapes.last()->isLoaded())
  {
      cout << "Model failed to be loaded. See log/terminal output." << endl;
      return;
  }

  ///Update case list in manager
  QStringList headingList;

  headingList << "Index" << "Case ID"; //!< Case Browser
  if(!caseManagerCreated)
    {
      caseTabIndex = manager->newTab("Case Browser", headingList);
      caseManagerCreated = true;
    }
  else
    manager->clearTab(caseTabIndex);

  cout << "Loading cases into browser." << endl;
  QList< int > cases = robustShapes.last()->getCaseIDs();
  for(int j = 0; j < cases.size(); j ++)
    {
      QStringList caseList;

      caseList << QString::number(j) << "Case " + QString::number(cases[j]);

      manager->addItem(caseTabIndex, caseList, Qt::NoItemFlags);
    }

  manager->show();
  cout << "Done." << endl;
}

void milxQtSSMPlugin::open(QString filename)
{
    if(filename.contains(".rssm", Qt::CaseInsensitive))
    {
        cout << "Loading as Robust Shape Model" << endl;
        addShapeModel( new milxQtRobustShapeModel );
            currentModel = robustShapes.last();
            robustShapes.last()->setConsole(console);
            robustShapes.last()->openModel(filename);
            connect(robustShapes.last(), SIGNAL(resultAvailable(milxQtRenderWindow*)), MainWindow, SLOT(display(milxQtRenderWindow*)));
            connect(robustShapes.last(), SIGNAL(resultAvailable(milxQtModel*)), MainWindow, SLOT(display(milxQtModel*)));
            connect(robustShapes.last(), SIGNAL(collectionAvailable(vtkPolyDataCollection*, QStringList&)), this, SLOT(passOnCollection(vtkPolyDataCollection*, QStringList&)));
    }
    else if(filename.contains(".ssm", Qt::CaseInsensitive))
    {
        cout << "Loading as Normal Shape Model" << endl;
        addShapeModel( new milxQtShapeModel );
            currentModel = shapes.last();
            shapes.last()->setConsole(console);
            shapes.last()->openModel(filename);
            connect(shapes.last(), SIGNAL(resultAvailable(milxQtRenderWindow*)), MainWindow, SLOT(display(milxQtRenderWindow*)));
            connect(shapes.last(), SIGNAL(resultAvailable(milxQtModel*)), MainWindow, SLOT(display(milxQtModel*)));
            connect(shapes.last(), SIGNAL(collectionAvailable(vtkPolyDataCollection*, QStringList&)), this, SLOT(passOnCollection(vtkPolyDataCollection*, QStringList&)));
    }
    else
    {
        cout << "File extension not supported." << endl;
        return;
    }

    dataName = filename;
    cout << "Loaded SSM File." << endl;
}

void milxQtSSMPlugin::save(QString filename)
{
    cout << "Saving as an SSM file" << endl;
    if(filename.contains(".rssm", Qt::CaseInsensitive))
        robustShapes.last()->saveModel(filename);
    else if(filename.contains(".ssm", Qt::CaseInsensitive))
        shapes.last()->saveModel(filename);
    else
    {
        cout << "File Extension not support!" << endl;
        return;
    }

    dataName = filename;
    cout << "Done." << endl;
}

milxQtRenderWindow* milxQtSSMPlugin::genericResult()
{
    if(!isPluginRobustWindow(currentModel))
    {
        cout << "Creating Normal SSM Model for display." << endl;
            shapes.last()->generateMeanModel();
            shapes.last()->setName(dataName);
            shapes.last()->generateRender();

        return shapes.last();
    }

    cout << "Creating Robust SSM Model for display." << endl;
        robustShapes.last()->generateMeanModel();
        robustShapes.last()->setName(dataName);
        robustShapes.last()->generateRender();

    return robustShapes.last();
}

milxQtModel* milxQtSSMPlugin::modelResult()
{
    if(!isPluginRobustWindow(currentModel))
    {
        QPointer<milxQtModel> model = shapes.last()->getMeanModel();
            model->setName(dataName + " - Mean Model");
            model->setDeletableOnClose(false);

        return model; //!< return mean model result
    }

    QPointer<milxQtModel> model = robustShapes.last()->getMeanModel();
        model->setName(dataName + " - Mean Model");
        model->setDeletableOnClose(false);

    return model; //!< return mean model result
}

bool milxQtSSMPlugin::isPluginRobustWindow(QWidget *window) //not in plugin-interface
{
    if(pluginRobustWindow(window) == 0)
        return false;
    else
        return true;
}

bool milxQtSSMPlugin::isPluginWindow(QWidget *window)
{
  if(pluginWindow(window) == 0)
    return false;
  else
    return true;
}

milxQtRobustShapeModel* milxQtSSMPlugin::pluginRobustWindow(QWidget *window) //not in plugin-interface
{
    if(window)
    {
        milxQtRobustShapeModel *ssmWindow = qobject_cast<milxQtRobustShapeModel *>(window);
        if(ssmWindow)
            return ssmWindow;
    }

    return 0;
}

milxQtShapeModel* milxQtSSMPlugin::pluginWindow(QWidget *window)
{
  if(window)
    {
      milxQtShapeModel *ssmWindow = qobject_cast<milxQtShapeModel *>(window);
      if(ssmWindow)
        return ssmWindow;
    }

  return 0;
}

void milxQtSSMPlugin::preStartTasks()
{
    cout << "Pre Start Tasks" << endl;
    if(dataName.contains(".rssm", Qt::CaseInsensitive))
    {
        cout << "Loading as Robust Shape Model" << endl;
        addShapeModel( new milxQtRobustShapeModel );
        currentModel = robustShapes.last();
    }
    else if(dataName.contains(".ssm", Qt::CaseInsensitive))
    {
        cout << "Loading as Normal Shape Model" << endl;
        addShapeModel( new milxQtShapeModel );
        currentModel = shapes.last();
    }
}

void milxQtSSMPlugin::postStartTasks()
{
    cout << "Post Start Tasks" << endl;
    if(dataName.contains(".rssm", Qt::CaseInsensitive))
    {
        robustShapes.last()->generateMeanModel();
        robustShapes.last()->generateRender();
    }
    else if(dataName.contains(".ssm", Qt::CaseInsensitive))
    {
        shapes.last()->generateMeanModel();
        shapes.last()->generateRender();
    }

    cout << "Creating Shape Model Display" << endl;
    milxQtRenderWindow *rnd = genericResult();
    cout << "Creating Mean Display" << endl;
    milxQtModel *mdl = modelResult();

    cout << "Updating Plugin" << endl;
    update();

    emit resultAvailable(rnd);
    emit resultAvailable(mdl);
}

void milxQtSSMPlugin::run()
{
    QMutexLocker locker(&mutex); //Lock memory

    ///Execute own thread work here
    int duration = 0;
    QTime *timer = new QTime;

    timer->start();
    cout << "Executing load." << endl;
		if(dataName.contains(".rssm", Qt::CaseInsensitive))
        {
            robustShapes.last()->setName(dataName);
            robustShapes.last()->openModel(dataName);
            robustShapes.last()->generateSSM();
        }
        else if(dataName.contains(".ssm", Qt::CaseInsensitive))
        {
            shapes.last()->setName(dataName);
            shapes.last()->openModel(dataName);
            shapes.last()->generateSSM();
        }

//		cout << "Creating Shape Model Display" << endl;
//        milxQtRenderWindow *rnd = genericResult();
//		cout << "Creating Mean Display" << endl;
//        milxQtModel *mdl = modelResult();
//
//		cout << "Updating Plugin" << endl;
//		update();
//
//        emit resultAvailable(rnd);
//        emit resultAvailable(mdl);
    cout << "Completed load." << endl;

    duration = timer->elapsed();
    cout << "Computation took: " << QString::number(duration/1000.0).toStdString() << " secs" << endl;
    delete timer;

    //exec();
}

void milxQtSSMPlugin::update()
{
    if(shapes.isEmpty() && modelManagerCreated)
        manager->clearTab(modelTabIndex);
    else if(!isPluginRobustWindow(currentModel))
        updateManager(shapes.last());
    else
        updateManager(robustShapes.last());

    if(!shapes.isEmpty() || !robustShapes.isEmpty())
        actionMultiModel->setDisabled(false);
    else
        actionMultiModel->setDisabled(true);
}

void milxQtSSMPlugin::updateManager(QWidget *newWindow)
{
    milxQtShapeModel *newModel = pluginWindow(newWindow);
    ///If fails, then neither new or old SSM so dont bother doing anything

    if(!newModel)
        return;

    if(!modelManagerCreated) //!< If model manager not created, create it
    {
        QStringList headings;

        headings << "Model Name" << "# Shapes" << "Weight" << "";

        modelTabIndex = manager->newTab("Model Browser", headings);

        modelManagerCreated = true;
    }

    ///If valid SSM, dont know if its new so check all windows and updated
    ///model manager tab
    manager->clearTab(modelTabIndex);
//    MainWindow->initialiseWindowTraversal();
//    for(int j = 0; j < MainWindow->getNumberOfWindows(); j ++)
//    {
//        milxQtShapeModel *newWindow = pluginWindow(MainWindow->nextWindow());
//
//        if(!newWindow)
//            continue;
//
//        QStringList modelEntry;
//        size_t n = newWindow->GetNumberOfShapes();
//
//        modelEntry << newWindow->strippedBaseName() << QString::number(n) << "1.0";
//        manager->addItem(modelTabIndex, modelEntry);
//    }

    ///Assume this plugin is the only loader of SSMs
    ///If loaded SSM, update the manager display
    foreach(QPointer<milxQtShapeModel> currModel, shapes)
    {
        QStringList modelEntry;
//        QPushButton *editButton = new QPushButton;
        size_t n = currModel->GetNumberOfShapes();

//        editButton->setText("Edit");

        modelEntry << currModel->strippedName() << QString::number(n) << "1.0"; //Last entry is the widget
//        manager->addItem(modelTabIndex, modelEntry, editButton, modelEntry.size()-1); //Put button at end
        manager->addItem(modelTabIndex, modelEntry, Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable); //Put button at end
    }
}

void milxQtSSMPlugin::multiModel()
{
    if(shapes.isEmpty())
        return;

    hybridShapeModel = new milxQtShapeModel(MainWindow);

    const int n = shapes.last()->GetNumberOfShapes();

    cout << "Creating a Multi-Model by Fusing " << n << " Shapes together." << endl;
    vtkPolyDataCollection* collection = vtkPolyDataCollection::New();
    for(int j = 0; j < n; j ++)
    {
        milxQtModel *hybridModel = new milxQtModel; //!< Deleted by SSM class

        foreach(QPointer<milxQtShapeModel> shape, shapes)
        {
            hybridModel->AddInput( shape->GetShape(j) );
//            hybridModel->AddInput( shape->GetAlignedShape(j) );

            qApp->processEvents(); ///Keep UI responsive
        }

        hybridModel->generateModel();
//        MainWindow->display(hybridModel);

//        hybridShapeModel->AddShape( hybridModel->GetOutput() );
        collection->AddItem( hybridModel->GetOutput() );

        qApp->processEvents(); ///Keep UI responsive
    }
    hybridShapeModel->SetInputCollection(collection);

    cout << "Display Hybrid Shape Model with " << hybridShapeModel->GetNumberOfShapes() << " shapes" << endl;
    addShapeModel(hybridShapeModel);
    currentModel = shapes.last();
    MainWindow->display( genericResult() );
    cout << "Displayed Hybrid Shape." << endl;
    MainWindow->display( modelResult() );

    update();
}

void milxQtSSMPlugin::focusedModel()
{
    int ret = wizard.exec();

    if(ret == QDialog::Rejected)
    {
        wizard.restart();
        return;
    }

    atlasFilename = txtAtlasName->text();

    //open atlas
    cout << "Opening atlas model" << endl;
    vtkSmartPointer<vtkPolyData> atlasMesh;
    milx::File::OpenModel(atlasFilename.toStdString(), atlasMesh);
    QPointer<milxQtModel> atlasModel = new milxQtModel;
    atlasModel->setName(atlasFilename);
    atlasModel->SetInput(atlasMesh);
    atlasModel->generateModel();
    emit resultAvailable(atlasModel);

    //open collection
    vtkPolyDataCollection* modelCollection = vtkPolyDataCollection::New(); //RSSM takes ownership, but not fully?
    /*std::vector<std::string> filenames;
    for (int j = 0; j < surfaceFilenames.size(); j ++)
        filenames.push_back(surfaceFilenames.at(j).toStdString());
    milx::File::OpenModelCollection(filenames, collection);*/
    QPointer<milxQtFile> file = new milxQtFile;
    if(!file->openModelCollection(modelCollection, surfaceFilenames))
    {
        cout << "Unable to read selected files." << endl;
        return;
    }

    //create model
    cout << "Computing Model..." << endl;
    SetInputCollection(modelCollection, atlasMesh, surfaceFilenames);

    //emit results for display
    emit resultAvailable(genericResult()); //shape model
    emit resultAvailable(modelResult()); //mean shape

    wizard.restart();
    update();
}

void milxQtSSMPlugin::showAtlasFileDialog()
{
    QSettings settings("Shekhar Chandra", "milxQt");
    //Add supported file entry
    QString exts = openModelExts.c_str();
    QString path = settings.value("recentPath").toString();

    QFileDialog *fileOpener = new QFileDialog;
    atlasFilename = fileOpener->getOpenFileName(&wizard,
                            tr("Select Atlas Surface to Open"),
                            path,
                            tr(exts.toStdString().c_str()) ); //!< \todo Check and validate extensions support at Open in Main class
    txtAtlasName->setText(atlasFilename);
}

void milxQtSSMPlugin::showSurfacesFileDialog()
{
  QSettings settings("Shekhar Chandra", "milxQt");
  //Add supported file entry
  QString exts = openModelExts.c_str();
  QString path = settings.value("recentPath").toString();

  QFileDialog *fileOpener = new QFileDialog;
  surfaceFilenames = fileOpener->getOpenFileNames(&wizard,
                                              tr("Select Training Surfaces to Open"),
                                              path,
                                              tr(exts.toStdString().c_str()) ); //!< \todo Check and validate extensions support at Open in Main class
  comboSurfaceNames->insertItems(0, surfaceFilenames);
}

void milxQtSSMPlugin::closedSSM(QWidget *win)
{
    milxQtShapeModel *tmpModel = pluginWindow(win);
    shapes.removeAll(tmpModel);
    update();
}

void milxQtSSMPlugin::passOnCollection(vtkPolyDataCollection *modelCollection, QStringList &filenames)
{
    cout << "Passing On Collection Signal" << endl;
    emit resultAvailable(modelCollection, filenames);
}

void milxQtSSMPlugin::createActions()
{
    actionMultiModel = new QAction(MainWindow);
    actionMultiModel->setText("Create &Multi-Model");
    actionMultiModel->setShortcut(tr("Ctrl+m"));
    actionFocusModel = new QAction(MainWindow);
    actionFocusModel->setText("Create &Focused Model");
    actionFocusModel->setShortcut(tr("Ctrl+f"));
}

void milxQtSSMPlugin::createMenu()
{
    menuSSM = new QMenu(MainWindow);
    menuSSM->setTitle("Shape Models");
    menuSSM->addAction(actionMultiModel);
    menuSSM->addAction(actionFocusModel);
    actionMultiModel->setDisabled(true);

    menuToAdd.append(menuSSM);
}

void milxQtSSMPlugin::createWizard()
{
    //use wizard to ask
    //intro
    QWizardPage *introPage = new QWizardPage;
    introPage->setTitle("Introduction");

    QLabel *label1 = new QLabel("This wizard will help you create a focused model from your training surfaces.");
    label1->setWordWrap(true);
    QLabel *label2 = new QLabel("You will be asked to provide an atlas/template mesh with scalars representing the weights for the model.");
    label2->setWordWrap(true);
    QLabel *label3 = new QLabel("You will be asked to provide training surfaces (with correspondences pre-computed) for the model.");
    label3->setWordWrap(true);
    QVBoxLayout *introLayout = new QVBoxLayout;
    introLayout->addWidget(label1);
    introLayout->addWidget(label2);
    introLayout->addWidget(label3);
    introPage->setLayout(introLayout);

    //ask for atlas mesh
    QWizardPage *atlasPage = new QWizardPage;
    atlasPage->setTitle("Atlas Surface with Weights");

    QLabel *label4 = new QLabel("Please provide an atlas/template mesh with scalars representing the weights for the model.");
    QLabel *label4b = new QLabel("Atlas/template mesh with scalars must have correspondence with training surfaces.");
    QLabel *label4c = new QLabel("Atlas/template mesh with scalars must not have zero values. Replace any with very small values.");
    label4->setWordWrap(true);
    //txtAtlasName class var
    //btnAtlasName class var
    btnAtlasName->setText("Browse...");
    QHBoxLayout *atlasNameLayout = new QHBoxLayout;
    atlasNameLayout->addWidget(txtAtlasName);
    atlasNameLayout->addWidget(btnAtlasName);
    QGroupBox *atlasGroupBox = new QGroupBox("Atlas/Template Surface Filename");
    atlasGroupBox->setLayout(atlasNameLayout);
    QVBoxLayout *atlasLayout = new QVBoxLayout;
    atlasLayout->addWidget(label4);
    atlasLayout->addWidget(label4b);
    atlasLayout->addWidget(label4c);
    atlasLayout->addWidget(atlasGroupBox);
    atlasPage->setLayout(atlasLayout);

    //ask for training surfaces
    QWizardPage *surfacesPage = new QWizardPage;
    surfacesPage->setTitle("Training Surfaces");

    QLabel *label5 = new QLabel("Please provide training surfaces (with correspondences pre-computed) for the model.");
    label5->setWordWrap(true);
    //comboSurfaceNames class var
    //btnSurfaceNames class var
    //btnClearSurfaceNames class var
    btnSurfaceNames->setText("Add...");
    btnClearSurfaceNames->setText("Clear");
    QFormLayout *surfaceNamesLayout = new QFormLayout;
    surfaceNamesLayout->addRow(btnClearSurfaceNames, btnSurfaceNames);
    surfaceNamesLayout->addRow(comboSurfaceNames);
    QGroupBox *surfacesGroupBox = new QGroupBox("Atlas/Template Surface Filename");
    surfacesGroupBox->setLayout(surfaceNamesLayout);
    QVBoxLayout *surfacesLayout = new QVBoxLayout;
    surfacesLayout->addWidget(label5);
    surfacesLayout->addWidget(surfacesGroupBox);
    surfacesPage->setLayout(surfacesLayout);

    //wizard class var
    wizard.addPage(introPage);
    wizard.addPage(atlasPage);
    wizard.addPage(surfacesPage);

    wizard.setWindowTitle("Focused Model Wizard");
}

void milxQtSSMPlugin::createConnections()
{
//    connect(MainWindow, SIGNAL(windowActivated(QWidget*)), this, SLOT(updateManager(QWidget*)));
    connect(actionMultiModel, SIGNAL(activated()), this, SLOT(multiModel()));
    connect(actionFocusModel, SIGNAL(activated()), this, SLOT(focusedModel()));

    connect(btnAtlasName, SIGNAL(clicked()), this, SLOT(showAtlasFileDialog()));
    connect(btnSurfaceNames, SIGNAL(clicked()), this, SLOT(showSurfacesFileDialog()));
    connect(btnClearSurfaceNames, SIGNAL(clicked()), comboSurfaceNames, SLOT(clear()));
}

//Q_EXPORT_PLUGIN2(SSMPlugin, milxQtSSMPluginFactory);
