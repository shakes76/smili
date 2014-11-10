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
#ifndef milxQtRobustShapeModel_H
#define milxQtRobustShapeModel_H

//milxSSM3D
#include <itkRobustStatisticalShapeModel.h>

#include <milxQtAliases.h>
#include <milxQtShapeModel.h>

typedef itk::RobustStatisticalShapeModel<double> RobustShapeModelType;

class MILXQT_PLUGIN_EXPORT milxQtRobustShapeModel : public milxQtShapeModel
{
    Q_OBJECT

public:
    /** Default constructor */
    milxQtRobustShapeModel(QWidget *theParent = 0, bool contextSystem = true);
    /** Default destructor */
    virtual ~milxQtRobustShapeModel();

public slots:
    /*!
        \fn milxQtRobustShapeModel::LoadModel(const QString filename)
        \brief Loads a robust model as an SSM (*.rssm) file. MILX-MSK like call.

        See itkRobustStatisticalShapeModel.h for RSSM file format details.
    */
    void LoadModel(const QString filename);
    /*!
        \fn milxQtRobustShapeModel::openModel(const QString filename)
        \brief Loads a robust model as an SSM (*.rssm) file. Alias same as LoadModel().

        See itkRobustStatisticalShapeModel.h for RSSM file format details.
    */
    inline virtual void openModel(const QString filename)
    {
        LoadModel(filename);
    }
    /*!
        \fn milxQtRobustShapeModel::SaveModel(const QString filename)
        \brief Saves a robust model as an SSM (*.rssm) file. MILX-MSK like call.

        See itkRobustStatisticalShapeModel.h for RSSM file format details.
    */
    inline virtual void SaveModel(const QString filename)
    {
        m_RobustSSM->Update();
        m_RobustSSM->SaveCompactModel(filename.toStdString().c_str());
    }
    /*!
        \fn milxQtRobustShapeModel::saveModel(const QString filename)
        \brief Saves a robust model as an SSM (*.rssm) file.

        See itkRobustStatisticalShapeModel.h for RSSM file format details.
    */
    inline virtual void saveModel(const QString filename)
    {
        SaveModel(filename);
    }

    /*!
        \fn milxQtRobustShapeModel::SetInputCollection(vtkPolyDataCollection* meshes)
        \brief Uses the collection of polydata to create shape model
    */
    virtual void SetInputCollection(vtkPolyDataCollection* meshes);
    /*!
        \fn milxQtRobustShapeModel::SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames)
        \brief Uses the collection of polydata to create shape model

        This version extracts case IDs from the filenames provided and displays them in the case browser.
    */
    virtual void SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames);
    /*!
        \fn milxQtRobustShapeModel::SetInputCollection(vtkPolyDataCollection* meshes, vtkPolyData *atlasSurface, QStringList &filenames)
        \brief Uses the collection of polydata to create focused shape model with weights from atlasSurface

        This version extracts case IDs from the filenames provided and displays them in the case browser.
    */
    virtual void SetInputCollection(vtkPolyDataCollection* meshes, vtkPolyData *atlasSurface, QStringList &filenames);

    /*!
        \fn milxQtRobustShapeModel::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);

    /**
        \fn milxQtRobustShapeModel::generateSSM()
        Generate the Statistical Shape Model. Assumes SetInputCollection() has already been called.
    */
    virtual void generateSSM();
    /**
        \fn milxQtRobustShapeModel::generateMeanModel(vtkSmartPointer<vtkPolyData> shape = NULL)
        Generates the display of the mean Model. Assumes generateSSM() has already been called.
    */
    virtual void generateMeanModel(vtkSmartPointer<vtkPolyData> shape = NULL);
    /**
        \fn milxQtRobustShapeModel::generateModes()
        Generate the display of the PCA modes within the SSM framework. Assumes atleast generateSSM() has already been called.

        Code for this member was contributed by Kaikai Shen, 2010.
    */
    virtual void generateModes();
    virtual void generateCollectionBasedOnMode();
    /**
        \fn milxQtRobustShapeModel::generateCorrespondences()
        Generate the display of the correspondences within the SSM framework as a hedgehog plot. Assumes atleast generateSSM() has already been called.
    */
    virtual void generateCorrespondences();

    /*!
        \fn milxQtRobustShapeModel::compactness()
        \brief Plot the compactness of the SSM.

        The compactness is basically how many modes are needed to represent what proportion of the model.
        High compactness is having a small number of modes representing a large amount of the model.
    */
    void compactness();
    /*!
        \fn milxQtRobustShapeModel::specificity()
        \brief Plot the specificity of the SSM.

        The specificity is given a randoming generated b-vector (shape parameters), how close would it be to the training population?
    */
    void specificity();
    /*!
        \fn milxQtRobustShapeModel::generalisability()
        \brief Plot the generalisability of the SSM.

        The generalisability is basically how good the model is at reconstruction using leave-one-out approach.
    */
    void generalisability();
    /*!
        \fn milxQtRobustShapeModel::eigenvalues()
        \brief Plot the eigenvalues of the model.
    */
    void eigenvalues();
    /*!
        \fn milxQtRobustShapeModel::eigenmodes()
        \brief Plot the primary eigenmodes of the model for each training shape.
    */
    void eigenmodes();
    /*!
        \fn milxQtRobustShapeModel::parameters()
        \brief Plot the training shape parameters of the model for each training shape (leave 1 out).
    */
    void parameters();

signals:
    /*!
        \fn milxQtRobustShapeModel::surfaceToImage(vtkPolyDataCollection*, QStringList&)
        \brief Emit signal to pass on collection.
    */
    void collectionAvailable(vtkPolyDataCollection*, QStringList&);
    /*!
        \fn milxQtRobustShapeModel::resultAvailable(milxQtRenderWindow*)
        \brief Send signal that Resultant render window is available for showing.
    */
    void resultAvailable(milxQtRenderWindow*);
    /*!
        \fn milxQtRobustShapeModel::resultAvailable(milxQtModel*)
        \brief Send signal that Resultant model is available for showing.
    */
    void resultAvailable(milxQtModel*);

protected:
    RobustShapeModelType::Pointer m_RobustSSM; //!< The Statistical Shape Model
    QList< int > m_caseIDs; //!< A list of case IDs corresponding to the models in the SSM

    virtual void reset();
    /*!
        \fn milxQtRobustShapeModel::createActions()
        \brief Create the actions for context menu etc.
    */
    void createActions();
    /*!
    	\fn milxQtRobustShapeModel::contextMenuEvent(QContextMenuEvent *event)
    	\brief The context menu setup member
    */
    void contextMenuEvent(QContextMenuEvent *event);

private:

};

#endif // milxQtRobustShapeModel_H
