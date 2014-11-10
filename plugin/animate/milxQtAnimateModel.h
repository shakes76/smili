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
#ifndef MILXQTANIMATEMODEL_H
#define MILXQTANIMATEMODEL_H

#include <QList>

#include <vtkPolyDataCollection.h>

#include <milxQtAliases.h>
#include <milxQtModel.h>

class MILXQT_PLUGIN_EXPORT milxQtAnimateModel : public milxQtModel
{
    Q_OBJECT

public:
    milxQtAnimateModel(QWidget *theParent = 0, bool contextSystem = true);
    virtual ~milxQtAnimateModel();

public slots:
    /*!
        \fn milxQtAnimateModel::SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames)
        \brief Animate the collection provided.

        This version extracts case IDs from the filenames provided and displays them in the case browser.
        If idIndex is given, the case IDs will be read and expected at this index (where index is the number of integer substrings found in the filename).
    */
    void SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames, const int idIndex = -1);
    /*!
        \fn milxQtAnimateModel::createMenu(QMenu *menu)
        \brief Create the menu for the data in this object. Used for context menu and file menus.
    */
    virtual void createMenu(QMenu *menu);
    void reset();
    void updateAnimation();
    void startAnimation();
    void pauseAnimation();
    void interval(int newInterval = 0);
    void intervalRotation(int newInterval = 0);
    void movie(QString filename = "", int frames = 0);
    inline void setInterval(const int delay)
    {   m_interval = delay; }
    inline int getInterval()
    {   return m_interval;  }
    inline size_t getCurrentFrame()
    {   return m_currentID; }
    inline bool isLoaded()
    {   return m_loaded;        }

protected:
    bool m_loaded; //!< Shapes loaded?
    bool m_pause; //!< Pause animation?
    int m_currentID; //!< Current frame or ID being displayed
    int m_interval; //!< delay between frames in ms
    int m_rotationInterval; //!< delay between frames in ms

    //Context Menu
    //menu defined in milxQtWindow
    //------------------
    QAction* startAct; //!< Action for starting animation
    QAction* pauseAct; //!< Action for stopping/pausing animation
    QAction* rotationAct; //!< Action for rotating view as well
    QAction* intervalAct; //!< Action for setting interval
    QAction* rotationIntervalAct; //!< Action for rotating view as well
    QAction* movieAct; //!< Action for recording movie

    vtkSmartPointer<vtkPolyDataCollection> m_meshes; //!< Collection of meshes to animate
    QList< int > m_caseIDs; //!< A list of case IDs corresponding to the models in the SSM

    QTimer timer; //!< Timer for the animation loop

private:
    void createActions();
    void createConnections();
    void contextMenuEvent(QContextMenuEvent *currentEvent);

};

#endif // MILXQTANIMATEMODEL_H
