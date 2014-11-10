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
#ifndef MILXQTSSMPLUGIN_H
#define MILXQTSSMPLUGIN_H

#include <QList>

#include <milxQtPluginInterface.h>
#include <milxQtMain.h>

#include "milxQtShapeModel.h"
#include "milxQtRobustShapeModel.h"
#include "milxQtManager.h"

/** \brief The Statistical Shape Model (SSM) plugin for milxQt.

    Allows the creation, display and modifications of SSMs from vtkPolyData collections or via the *.ssm files.

    The plugin also has two managers, a model and a case manager. The former shows a list of all SSMs loaded and their weights (for usage with
    multi-model SSMs) and the latter shows the case list for the current model if created via collections.
 */
class MILXQT_PLUGIN_EXPORT milxQtSSMPlugin : public milxQtPluginInterface
{
    Q_OBJECT

public:
    milxQtSSMPlugin(QObject *theParent = 0);
    ~milxQtSSMPlugin();

    QString name();

    inline bool hasOpenSupport()
    {   return true;    }
    QString openFileSupport();
    QStringList openExtensions();
    inline bool hasSaveSupport()
    {   return true;    }
    QString saveFileSupport();
    QStringList saveExtensions();

    inline bool hasCollectionSupport()
    {   return true;    }
    void SetInputCollection(vtkPolyDataCollection* collection, QStringList &filenames);
    void SetInputCollection(vtkPolyDataCollection* collection, vtkPolyData *atlasSurface, QStringList &filenames);

    void open(QString filename);
    void save(QString filename);

    milxQtRenderWindow* genericResult();
    milxQtModel* modelResult();
    inline milxQtImage* imageResult()
    {   return NULL;    } //No image result
    inline QDockWidget* dockWidget()
    {   return dock;    }
    inline Qt::DockWidgetArea dockDefaultArea()
    {   return Qt::RightDockWidgetArea;    }

    bool isPluginRobustWindow(QWidget *window);
    bool isPluginWindow(QWidget *window);

    /**
    Casts window to a ShapeModel class after performing relevant checks
    */
    milxQtRobustShapeModel* pluginRobustWindow(QWidget *window);
    milxQtShapeModel* pluginWindow(QWidget *window);

public slots:
    virtual void loadExtension() {}
    virtual void update();
    void updateManager(QWidget *newWindow);
    void multiModel();
    void focusedModel();
    void showAtlasFileDialog();
    void showSurfacesFileDialog();
    void closedSSM(QWidget *win);
    void passOnCollection(vtkPolyDataCollection *modelCollection, QStringList &filenames);

    void preStartTasks();
    void postStartTasks();

protected:
    bool modelManagerCreated; //!< Model Manager tab has been created?
    bool caseManagerCreated; //!< Model Manager tab has been created?
    int modelTabIndex;
    int caseTabIndex;

    QList< QPointer<milxQtShapeModel> > shapes;
    QList< QPointer<milxQtRobustShapeModel> > robustShapes;
    QWidget *currentModel; //!< Current model being viewed/processed
    milxQtShapeModel *hybridShapeModel; //!< Hybrid model
    QPointer<milxQtManager> manager; //!< Manager widget
    QPointer<QDockWidget> dock; //!< Dock widget

    QPointer<milxQtMain> MainWindow;

    QMenu* menuSSM; //!< SSM menu
    //----SSM---- (hierarchical deletion)
    QAction* actionMultiModel; //!< hybrid model action
    QAction* actionFocusModel; //!< focus model action

    //Focus model variables
    QWizard wizard;
    QLineEdit *txtAtlasName;
    QPushButton *btnAtlasName;
    QString atlasFilename;
    QListWidget *comboSurfaceNames;
    QPushButton *btnSurfaceNames;
    QPushButton *btnClearSurfaceNames;
    QStringList surfaceFilenames;

    void run();
    inline void addShapeModel(milxQtShapeModel *newShapeModel)
    {
        cout << "Adding Normal Shape Model to System." << endl;
        shapes.append(newShapeModel);
        currentModel = qobject_cast<QWidget *>(newShapeModel);
        if(isPluginWindow(currentModel))
            cout << "Successfully added normal model to system" << endl;
        connect(shapes.last(), SIGNAL(closing(QWidget *)), this, SLOT(closedSSM(QWidget *)));
    }
    inline void addShapeModel(milxQtRobustShapeModel *newShapeModel)
    {
        cout << "Adding Robust Shape Model to System." << endl;
        robustShapes.append(newShapeModel);
        currentModel = qobject_cast<QWidget *>(newShapeModel);
        connect(robustShapes.last(), SIGNAL(closing(QWidget *)), this, SLOT(closedSSM(QWidget *)));
    }
    void createActions();
    void createMenu();
    void createWizard();
    void createConnections();

private:

};

class MILXQT_PLUGIN_EXPORT milxQtSSMPluginFactory: public QObject, public milxQtPluginFactory
{
    Q_OBJECT
    Q_INTERFACES(milxQtPluginFactory)

public:
    milxQtPluginInterface* newPlugin(QObject *theParent = 0)
    {   return new milxQtSSMPlugin(theParent);  }
};

#endif // MILXQTSSMPLUGIN_H
