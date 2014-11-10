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
#ifndef MILXQTSHAPEMODEL_H
#define MILXQTSHAPEMODEL_H

//milxSSM3D
#include <itkStatisticalShapeTimepointModel.h>

#include <milxQtAliases.h>
#include <milxQtRenderWindow.h>
#include <milxQtModel.h>

typedef itk::StatisticalShapeModelBase<double> ShapeModelBaseType;
typedef itk::StatisticalShapeTimepointModel<double> ShapeModelType;

class MILXQT_PLUGIN_EXPORT milxQtShapeModel : public milxQtRenderWindow
{
    Q_OBJECT

public:
    /** Default constructor */
    milxQtShapeModel(QWidget *theParent = 0, bool contextSystem = true);
    /** Default destructor */
    virtual ~milxQtShapeModel();

    /*!
        \fn milxQtShapeModel::isLoaded()
        \brief Is the model loaded successfully?
    */
    inline virtual bool isLoaded()
    {   return m_loaded;    }

public slots:
    /*!
        \fn milxQtShapeModel::LoadModel(const QString filename)
        \brief Loads a model as an SSM (*.ssm) file. MILX-MSK like call.

        See itkStatisticalShapeModel.h for SSM file format details.
    */
    void LoadModel(const QString filename);
    /*!
        \fn milxQtShapeModel::openModel(const QString filename)
        \brief Loads a model as an SSM (*.ssm) file. Alias same as LoadModel().

        See itkStatisticalShapeModel.h for SSM file format details.
    */
    inline virtual void openModel(const QString filename)
    {
        LoadModel(filename);
    }
    /*!
        \fn milxQtShapeModel::SaveModel(const QString filename)
        \brief Saves a model as an SSM (*.ssm) file. MILX-MSK like call.

        See itkStatisticalShapeModel.h for SSM file format details.
    */
    inline virtual void SaveModel(const QString filename)
    {
        m_StandardSSM->Update();
        m_StandardSSM->SaveModel(filename.toStdString().c_str());
    }
    /*!
        \fn milxQtShapeModel::saveModel(const QString filename)
        \brief Saves a model as an SSM (*.ssm) file.

        See itkStatisticalShapeModel.h for SSM file format details.
    */
    inline virtual void saveModel(const QString filename)
    {
        SaveModel(filename);
    }

    /*!
        \fn milxQtShapeModel::SetInputCollection(vtkPolyDataCollection* meshes)
        \brief Uses the collection of polydata to create shape model
    */
    virtual void SetInputCollection(vtkPolyDataCollection* meshes);
    /*!
        \fn milxQtShapeModel::SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames)
        \brief Uses the collection of polydata to create shape model

        This version extracts case IDs from the filenames provided and displays them in the case browser.
    */
    virtual void SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames);
    /*!
        \fn milxQtShapeModel::SetPrecision(float modePrecision)
        \brief Assigns the reduction factor of the modes, i.e the number of modes required to quantify modePrecision percent of the model.
    */
    inline void SetPrecision(float modePrecision)
    {
        m_SSM->SetPrecision(modePrecision);
    }
    /*!
        \fn milxQtShapeModel::AddShape(vtkPolyData *shape)
        \brief Adds the shape to the SSM, once done adding shapes, call generateSSM() to get the SSM.
    */
    inline void AddShape(vtkPolyData *shape)
    {   m_SSM->AddShape(shape);  }
    /*!
        \fn milxQtShapeModel::GetShape(unsigned int index)
        \brief Gets the shape at index in the SSM.
    */
    inline vtkPolyData* GetShape(unsigned int index)
    {   return m_SSM->GetShape(index);  }
    /*!
        \fn milxQtShapeModel::GetAlignedShape(unsigned int index)
        \brief Gets the shape at index in the SSM.
    */
    inline vtkPolyData* GetAlignedShape(unsigned int index)
    {   return m_SSM->GetProcrustesAlignedSurface(index);  }
    /*!
        \fn milxQtShapeModel::GetProcrustesAlignedSurface(unsigned int index)
        \brief Gets the shape at index in the SSM. Same as GetAlignedShape().
    */
    inline vtkPolyData* GetProcrustesAlignedSurface(unsigned int index)
    {   return m_SSM->GetProcrustesAlignedSurface(index);  }

    /*!
        \fn milxQtShapeModel::GetMeanShape()
        \brief Returns the mean shape of the Statistical Shape Model. Assumes SSM has been generated (via generateSSM()).
    */
    inline vtkPolyData* GetMeanShape()
    {
        if(m_modelled) return m_SSM->GetMeanShape();
        return NULL;
    }
    /*!
        \fn milxQtShapeModel::getMeanModel()
        \brief Returns the mean model of the Statistical Shape Model. Assumes SSM has been generated (via generateSSM()).
    */
    QPointer<milxQtModel> getMeanModel();
    /*!
        \fn milxQtShapeModel::GetNumberOfModes()
        \brief Returns the number of modes in the Statistical Shape Model according to the precision set. Assumes SSM has been generated (via generateSSM()).
    */
    inline int GetNumberOfModes()
    {
        if(m_modelled) return m_SSM->GetNumberOfModes();
        return -1;
    }
    /*!
        \fn milxQtShapeModel::GetNumberOfShapes()
        \brief Returns the number of shapes in the Statistical Shape Model. Assumes SSM has been generated (via generateSSM()).
    */
    inline int GetNumberOfShapes()
    {   return m_SSM->GetNumberOfShapes();   }
    /*!
        \fn milxQtShapeModel::getCaseIDs()
        \brief Returns the case IDs for the SSM. Assumes that an SSM has been loaded.
    */
    inline QList< int > getCaseIDs()
    {   return m_caseIDs;   }
    /**
      \brief Get the data, depending on whats been generated, return that dataset
    */
    virtual vtkDataSet* GetDataSet();

    /*!
        \fn milxQtShapeModel::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);

    /**
        \fn milxQtShapeModel::generateSSM()
        Generate the Statistical Shape Model. Assumes SetInputCollection() has already been called.
    */
    virtual void generateSSM();
    /**
        \fn milxQtShapeModel::generateMeanModel(vtkSmartPointer<vtkPolyData> shape = NULL)
        Generates the display of the mean Model. Assumes generateSSM() has already been called.
    */
    virtual void generateMeanModel(vtkSmartPointer<vtkPolyData> shape = NULL);
    /**
        \fn milxQtShapeModel::generateModels(const bool display = true)
        Generate the display of the original Models. Assumes generateSSM() has already been called.
    */
    void generateModels(const bool display = true);
    /**
        \fn milxQtShapeModel::generateAlignedModels(const bool display = true)
        Generate the display of the aligned Models within the SSM framework. Assumes generateSSM() has already been called.

        If display is true, all the models are rendered into window. Could be slow for a large set of models.
    */
    void generateAlignedModels(const bool display = true);
    /**
        \fn milxQtShapeModel::generateModes()
        Generate the display of the PCA modes within the SSM framework. Assumes atleast generateSSM() has already been called.

        Code for this member was contributed by Kaikai Shen, 2010.
    */
    virtual void generateModes();
    virtual void generateCollectionBasedOnMode();
    /**
        \fn milxQtShapeModel::generateVectorField()
        Generate the display of the PCA modes within the SSM framework as a vector field. Assumes atleast generateModes() has already been called.

        Code for this member was contributed by Kaikai Shen, 2010.
    */
    void generateVectorField();
    /**
        \fn milxQtShapeModel::generateTensorField()
        Generate the display of the PCA modes within the SSM framework as a tensor field. Assumes atleast generateModes() has already been called.

        Code for this member was contributed by Kaikai Shen, 2010.
    */
    void generateTensorField();
    /**
        \fn milxQtShapeModel::generateCorrespondences()
        Generate the display of the correspondences within the SSM framework as a hedgehog plot. Assumes atleast generateSSM() has already been called.
    */
    void generateCorrespondences();

    /*!
        \fn milxQtShapeModel::mean()
        \brief Show the mean model in the window.
    */
    void mean();
    /*!
        \fn milxQtShapeModel::aligned()
        \brief Show the aligned models in the window as a set of coloured point clouds.
    */
    void aligned();
    /*!
        \fn milxQtShapeModel::original()
        \brief Show the original models in the window as a set of coloured point clouds.
    */
    void original();
    /*!
        \fn milxQtShapeModel::modesAsVectors()
        \brief Show the modes variation as a vector field.
    */
    void modesAsVectors();
    /*!
        \fn milxQtShapeModel::modesAsTensors()
        \brief Show the modes variation as a ellipsoidal tensor field.
    */
    void modesAsTensors();
    /*!
        \fn milxQtShapeModel::modesAsCollection()
        \brief Show the modes variation as a collection of surfaces.
    */
    void modesAsCollection();
    /*!
        \fn milxQtShapeModel::procrustes()
        \brief Change the procrustes mode of the SSM (causes recomputation of the model).
    */
    void procrustes();
    /*!
        \fn milxQtShapeModel::alignment()
        \brief Output each aligned shape as a PNG image in a directory and prefix provided by the user (internally via file dialog).
    */
    void alignment();
    /*!
        \fn milxQtShapeModel::alignedMeshes()
        \brief Output each aligned shape as a poly data file in a directory and prefix provided by the user (internally via file dialog).
    */
    void alignedMeshes();
    /*!
        \fn milxQtShapeModel::originalMeshes()
        \brief Output each original shape as a poly data file in a directory and prefix provided by the user (internally via file dialog).
    */
    void originalMeshes();
    /*!
        \fn milxQtShapeModel::pointIds()
        \brief Output each aligned shape coloured by the point IDs as a PNG image in a directory and prefix provided by the user (internally via file dialog).
    */
    void pointIds();
    /*!
        \fn milxQtShapeModel::coordinates()
        \brief Output each aligned shape coloured by the coordinates as a PNG image in a directory and prefix provided by the user (internally via file dialog).

        Same regions (i.e. regions with correspondences) will have the same colour so that mis-matched points show up as coloured points in the wrong areas of the model.
    */
    void coordinates();
    /*!
        \fn milxQtShapeModel::replaceOriginal()
        \brief Replace the original (unaligned) meshes with the aligned ones. The new meshes are unnormalised as compared to the aligned meshes which are normalised.

        Useful when aligned meshes are accessed by GetShape() of the SSM class or unnormalised results are required.
        The rescaling of each aligned mesh is done based on the centroid size of the corresponding original mesh.
    */
    void replaceOriginal();
    /*!
        \fn milxQtShapeModel::correspondences()
        \brief Compare points from each shape to measure the correspondence as a hedgehog plot.
    */
    void correspondences();

    /*!
        \fn milxQtShapeModel::compactness()
        \brief Plot the compactness of the SSM.

        The compactness is basically how many modes are needed to represent what proportion of the model.
        High compactness is having a small number of modes representing a large amount of the model.
    */
    void compactness();
    /*!
        \fn milxQtShapeModel::specificity()
        \brief Plot the specificity of the SSM.

        The specificity is given a randoming generated b-vector (shape parameters), how close would it be to the training population?
    */
    void specificity();
    /*!
        \fn milxQtShapeModel::generalisability()
        \brief Plot the generalisability of the SSM.

        The generalisability is basically how good the model is at reconstruction using leave-one-out approach.
    */
    void generalisability();
    /*!
        \fn milxQtShapeModel::eigenvalues()
        \brief Plot the eigenvalues of the model.
    */
    void eigenvalues();
    /*!
        \fn milxQtShapeModel::eigenmodes()
        \brief Plot the primary eigenmodes of the model for each training shape.
    */
    void eigenmodes();
    /*!
        \fn milxQtShapeModel::parameters()
        \brief Plot the training shape parameters of the model for each training shape (leave 1 out).
    */
    void parameters();

signals:
    /*!
        \fn milxQtShapeModel::surfaceToImage(vtkPolyDataCollection*, QStringList&)
        \brief Emit signal to compute surface to image process.
    */
    void collectionAvailable(vtkPolyDataCollection*, QStringList&);
    /*!
        \fn milxQtShapeModel::resultAvailable(milxQtRenderWindow*)
        \brief Send signal that Resultant render window is available for showing.
    */
    void resultAvailable(milxQtRenderWindow*);
    /*!
        \fn milxQtShapeModel::resultAvailable(milxQtModel*)
        \brief Send signal that Resultant model is available for showing.
    */
    void resultAvailable(milxQtModel*);

protected:
    bool m_loaded; //!< Shapes loaded?
    bool m_modelled; //!< SSM Modelled?
    bool m_meaned; //!< Mean generated?
    bool m_shapesModelled; //!< Shapes were generated?
    bool m_alignedShapesModelled; //!< Shapes were generated?
    bool m_modes; //!< Modes been generated?
    bool m_vector; //!< Vectors generated?
    bool m_tensor; //!< Tensors generated?
    bool m_correspond; //!< Hedgehog generated?
    int m_mode; //!< Mode being viewed

    //----Context----
    QAction* actionMean; //!< mean action
    QAction* actionAligned; //!< aligned points action
    QAction* actionOriginal; //!< original points action
    QAction* actionModesAsVectors; //!< modes action
    QAction* actionModesAsTensors; //!< modes action
    QAction* actionModesAsCollection; //!< modes collection action
    //-------------
    QAction* actionAlignment; //!< alignment check action
    QAction* actionOriginalMeshes; //!< original meshes output action
    QAction* actionAlignedMeshes; //!< alignmed meshes output action
    QAction* actionPointIDs; //!< point ids check action
    QAction* actionCoordinates; //!< coordinates check action
    QAction* actionCorrespond; //!< correspondence check action
    //-------------
    QAction* actionReplaceOriginal; //!< correspondence check action
    //----Plot----
    QMenu* plotMenu; //!< Plot menu for SSM analysis
    QAction* actionCompact; //!< compactness plot action
    QAction* actionSpecificity; //!< specificity plot action
    QAction* actionGeneralise; //!< generalisibility error plot action
    QAction* actionValues; //!< primary modes plot action
    QAction* actionModes; //!< primary modes plot action
    QAction* actionParameters; //!< training shape parameters plot action
    //----Procrustes----
    QMenu* alignMenu; //!< Plot menu for SSM alignment
    QActionGroup* alignGroup; //!< Procrustes align gp action
    QAction* actionRigid; //!< rigid Procrustes action
    QAction* actionSimilarity; //!< similarity Procrustes action
    QAction* actionAffine; //!< affine Procrustes action

    ShapeModelType::Pointer m_StandardSSM; //!< The Statistical Shape Model
    ShapeModelBaseType::Pointer m_SSM; //!< The Statistical Shape Model Base Class

    QList< int > m_caseIDs; //!< A list of case IDs corresponding to the models in the SSM
    QList< QPointer<milxQtModel> > m_models; //!< Original Models in the SSM
    QList< QPointer<milxQtModel> > m_alignedModels; //!< Aligned Models in the SSM

    QPointer<milxQtModel> m_meanModel; //!< Model of the mean
    QPointer<milxQtModel> m_modesVectorModel; //!< Vector Model of the modes
    QPointer<milxQtModel> m_modesTensorModel; //!< Tensor Model of the modes
    QPointer<milxQtModel> m_correspondences; //!< Hedgehog display of the correspondences

    /*!
        \fn milxQtShapeModel::updateLookupTable()
        \brief Update the colour maps depending on whats been generated for the model.
    */
    virtual void updateLookupTable();
    virtual void reset();
    void outputSnapshots(const QString filename);
    /*!
        \fn milxQtShapeModel::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
        \fn milxQtShapeModel::createConnections()
        \brief Create the connections for context menu etc.
    */
    void createConnections();
    /*!
        \fn milxQtShapeModel::setupTooltips()
        \brief Assign tooltips to each action for users benefit. The tooltips explains the function of each action
    */
    void setupTooltips();
    /*!
    	\fn milxQtShapeModel::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);

private:

};

#endif // MILXQTSHAPEMODEL_H
