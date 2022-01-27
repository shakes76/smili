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
#include "milxQtAnimateModel.h"

#include <time.h>
//Qt 
#include <QMenu>
#include <QDialog>
#include <QInputDialog>
#include <QFileDialog>
#include <QFormLayout>
#include <QComboBox>
#include <QPushButton>
//VTK
#include <vtkWindowToImageFilter.h>
#include <vtkFFMPEGWriter.h>
#include <vtkCamera.h>

milxQtAnimateModel::milxQtAnimateModel(QWidget *theParent, bool contextSystem) : milxQtModel(theParent, contextSystem)
{
    m_loaded = false;
    m_pause = true;
    m_currentID = 0;
    m_interval = 250;
    m_rotationInterval = 5;

    milxQtWindow::prefix = "aniMdl: ";

    createActions();

    createConnections();
}

milxQtAnimateModel::~milxQtAnimateModel()
{
    //dtor
}

void milxQtAnimateModel::SetInputCollection(vtkPolyDataCollection* meshes, QStringList &filenames, const int idIndex)
{
    const int n = meshes->GetNumberOfItems();

    working(-1);
    m_meshes = meshes;

    printInfo("Extracting IDs");
    ///Setup temporary combo widget for selecting case ids when there are
    ///Multiple integer matches in the filename
    int index = 0;
    QDialog *casePossibilities = new QDialog(this);
    QComboBox *integerValues = new QComboBox(this);
    QPushButton *okButton = new QPushButton(this);
        okButton->setText("Ok");
    QLabel *lblMessage = new QLabel(this);
        lblMessage->setText("Choose Case ID from possibilities.");
    QFormLayout *formLayout = new QFormLayout(this);
        formLayout->addRow(lblMessage);
        formLayout->addRow(tr("&Select Case ID from first file: "), integerValues);
        formLayout->addRow(okButton);
        casePossibilities->setLayout(formLayout);

    connect(okButton, SIGNAL(clicked(bool)), casePossibilities, SLOT(accept()));

    meshes->InitTraversal();
    milxQtModel::SetInput(meshes->GetNextItem());
    reset();
    updateAnimation();

    for(int j = 0; j < n; j ++)
    {
        QFileInfo fi(filenames[j]);
        QRegExp rx("(\\d+)", Qt::CaseSensitive, QRegExp::RegExp2); ///Capture integer expression

        ///Read all integers from filename and save in list
        int pos = 0;
        QStringList list;
        while ((pos = rx.indexIn(fi.baseName(), pos)) != -1)
        {
            list << rx.cap(1);
            pos += rx.matchedLength(); //Move
        }

        ///If more than one integer in filename, ask user to select correct integer
        if(j == 0)
        {
            if(list.size() > 1 && idIndex < 0)
            {
                printInfo("Please choose Case ID from possibilities for first file.");
                for(int k = 0; k < list.size(); k ++)
                    integerValues->addItem(list[k]);

                casePossibilities->exec();

                index = integerValues->currentIndex();
            }
            else if(list.size() > 0 && idIndex < list.size() && idIndex > 0)
            {
                index = idIndex;
            }
            else if(list.size() > 0 && idIndex < list.size() && idIndex < 0)
            {
                index = 0;
            }
            else
            {
              printError("Index provided for expected case IDs is incorrect. Ignoring operation.");
              m_loaded = false;
              return;
            }
        }

        int caseID = 0;
        printDebug("Frame Index to be used: " + QString::number(index));
        if(!list.isEmpty())
            caseID = list[index].toInt(); //!< Extract case ID from filename by using integer user provided (or only one present)

        m_caseIDs.append(caseID); //!< Add case ID to list

        qApp->processEvents(); ///Keep UI responsive
    }
    done(-1);

    qDebug() << "IDs: " << m_caseIDs << endl;

    printInfo("Starting Animation Loop");
    connect(&timer, SIGNAL(timeout()), this, SLOT(updateAnimation()));
    timer.start(m_interval);

    m_loaded = true;
}

void milxQtAnimateModel::createMenu(QMenu *menu)
{
    if(!menu)
        return;

    menu->clear();
    foreach(QAction *currAct, milxQtModel::actionsToAdd)
    {
        menu->addAction(currAct);
    }
    foreach(QMenu *currMenu, menusToAdd)
    {
        menu->addMenu(currMenu);
    }
    menu->addSeparator()->setText("Animation");
    menu->addAction(startAct);
    menu->addAction(pauseAct);
    menu->addAction(rotationAct);
    menu->addAction(intervalAct);
    menu->addAction(rotationIntervalAct);
    menu->addAction(movieAct);
    menu->addSeparator()->setText(tr("Display"));
    menu->addAction(milxQtModel::pointsAct);
    menu->addAction(milxQtModel::wireframeAct);
    menu->addAction(milxQtModel::surfaceAct);
    ///Change View of Volume
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
    milxQtModel::showMenu = menu->addMenu("Show"); //!< Only exists for the duration of the context selection
    milxQtModel::showMenu->addAction(milxQtModel::normalsAct);
    milxQtModel::showMenu->addAction(milxQtModel::centroidAct);
    milxQtModel::showMenu->addAction(milxQtModel::outlineAct);
    milxQtModel::showMenu->addAction(milxQtModel::cubeAxesAct);
    milxQtModel::showMenu->addMenu(milxQtModel::arraysMenu());
    menu->addSeparator()->setText(tr("Properties"));
    menu->addAction(milxQtModel::colourAct);
    menu->addAction(milxQtModel::interpAct);
    menu->addSeparator();
    menu->addAction(milxQtModel::scaleAct);
    menu->addMenu(milxQtRenderWindow::windowPropertiesMenu);
    menu->addAction(milxQtRenderWindow::refreshAct);
    menu->addAction(milxQtRenderWindow::resetAct);
}

void milxQtAnimateModel::reset()
{
    m_currentID = 0;
    m_meshes->InitTraversal();
}

void milxQtAnimateModel::updateAnimation()
{
    if(m_pause)
        return;

    const int n = m_meshes->GetNumberOfItems();

    m_currentID ++;
    if(m_currentID >= n) //Wrap around animation frames
        reset();

    //camera rotation
    vtkCamera* srcCamera = milxQtRenderWindow::GetRenderer()->GetActiveCamera(); //!< Get callers camera
    if(rotationAct->isChecked())
        srcCamera->Azimuth(m_rotationInterval);

    printDebug("Animating frame " + QString::number(m_currentID) + " of " + QString::number(n));
    milxQtModel::setName("ID " + QString::number(m_caseIDs[m_currentID]));

    vtkSmartPointer<vtkPolyData> mesh = m_meshes->GetNextItem();
    milxQtModel::SetInput(mesh);
    milxQtModel::refresh();
    qApp->processEvents();
}

void milxQtAnimateModel::startAnimation()
{
    m_pause = false;
    reset();
    updateAnimation();
}

void milxQtAnimateModel::pauseAnimation()
{
    if(!m_pause)
        m_pause = true;
    else
        m_pause = false;

    updateAnimation();
}

void milxQtAnimateModel::interval(int newInterval)
{
    bool ok = false;

    if(newInterval == 0)
    {
        newInterval = QInputDialog::getInt(this, tr("Please Provide the delay between frames"),
                              tr("Interval:"), m_interval, 0, 1000000, 1, &ok);
    }
    else
        ok = true;

    if(ok && newInterval > 0)
    {
        m_interval = newInterval;
        timer.start(m_interval);
    }
}

void milxQtAnimateModel::intervalRotation(int newInterval)
{
    bool ok = false;

    if(newInterval == 0)
    {
        newInterval = QInputDialog::getDouble(this, tr("Please Provide the angle increment"),
                              tr("Angle:"), 5, 0, 360, 2, &ok);
    }
    else
        ok = true;

    if(ok && newInterval > 0)
        m_rotationInterval = newInterval;
}

void milxQtAnimateModel::movie(QString filename, int frames)
{
    m_pause = true; //Force pause of animation
    timer.stop(); //Pause current animation

    if(filename.isEmpty())
    {
        QSettings settings("Shekhar Chandra", "milxQt");
        QString path = settings.value("recentPath").toString();
        QFileDialog *fileSaver = new QFileDialog(this);
        filename = fileSaver->getSaveFileName(this,
                                  tr("Select File Name to Save"),
                                  path,
                                  tr("Movie Files (*.avi)"));
    }

    if(!filename.isEmpty())
    {
        const int n = m_meshes->GetNumberOfItems();
        const int frameRate = 1000.0/m_interval;
        bool ok = false;

        if(frames == 0)
        {
            frames = static_cast<int>( QInputDialog::getInt(this, tr("Please Provide the total frames to write"),
                                  tr("Frames:"), n, 0, 8192, 1, &ok) );

            if(!ok)
                return;
        }

        milxQtRenderWindow::OffScreenRenderingOn(); //Ensure no UI stuff interfers

        vtkWindowToImageFilter *windowToImage = vtkWindowToImageFilter::New(); //Manually deleted deliberately
            windowToImage->SetInput(milxQtRenderWindow::GetRenderWindow());
            linkProgressEventOf(windowToImage);

        vtkSmartPointer<vtkFFMPEGWriter> writer = vtkSmartPointer<vtkFFMPEGWriter>::New();
            writer->SetFileName(filename.toStdString().c_str());
            writer->SetInputConnection(windowToImage->GetOutputPort());
            writer->SetRate(frameRate);
            linkProgressEventOf(writer);
            writer->Start();

        //camera rotation
        vtkCamera* srcCamera = milxQtRenderWindow::GetRenderer()->GetActiveCamera(); //!< Get callers camera

        printDebug("Movie Write Begin");
        for(int j = 0; j < frames; j ++)
        {
            printDebug("Loading Frame " + QString::number(j));
            if( (j % n) == 0)
                reset(); //cycle through collection if frames > n

            vtkPolyData *mesh = m_meshes->GetNextItem();

            milxQtModel::setName("ID " + QString::number(m_caseIDs[j % n]));
            milxQtModel::SetInput(mesh);
            milxQtModel::refresh(); //update display
            qApp->processEvents();

            windowToImage->Update();

            printDebug("Writing Frame " + QString::number(j));
            writer->SetInputConnection(windowToImage->GetOutputPort()); //Force update of pointer
            writer->Write();

            printDebug("Prep next Frame ");
            windowToImage->Delete();
                windowToImage = vtkWindowToImageFilter::New();
                windowToImage->SetInput(milxQtRenderWindow::GetRenderWindow());
                linkProgressEventOf(windowToImage);
                windowToImage->Update();

            if(rotationAct->isChecked())
                srcCamera->Azimuth(m_rotationInterval);
            printDebug("Finish Loop");
        }
        printDebug("Movie Written Successfully");
        writer->End();
        printDebug("Movie File Closed Sucessfully");
        windowToImage->Delete();

        milxQtRenderWindow::OffScreenRenderingOff();
        refresh();

        printInfo("Movie Operation Completed");
    }

    timer.start(m_interval); //Restart timer
}

void milxQtAnimateModel::createActions()
{
    contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection
    contextMenu->setTitle("Animation");

    startAct = new QAction(this);
        startAct->setText("Start/Restart");
        startAct->setShortcut(tr("Alt+s"));
    pauseAct = new QAction(this);
        pauseAct->setText("Pause/Unpause");
        pauseAct->setShortcut(tr("Alt+p"));
    intervalAct = new QAction(this);
        intervalAct->setText("Change Interval");
        intervalAct->setShortcut(tr("Alt+i"));
    rotationIntervalAct = new QAction(this);
        rotationIntervalAct->setText("Change Rotation Interval");
        rotationIntervalAct->setShortcut(tr("Shift+Alt+r"));
    rotationAct = new QAction(this);
        rotationAct->setText("Rotate View");
        rotationAct->setShortcut(tr("Alt+r"));
        rotationAct->setCheckable(true);
        rotationAct->setChecked(false);
    movieAct = new QAction(this);
        movieAct->setText("Record as Movie");
        movieAct->setShortcut(tr("Alt+m"));
}

void milxQtAnimateModel::createConnections()
{
    //Operations
    connect(startAct, SIGNAL(triggered()), this, SLOT(startAnimation()));
    connect(pauseAct, SIGNAL(triggered()), this, SLOT(pauseAnimation()));
    connect(intervalAct, SIGNAL(triggered()), this, SLOT(interval()));
    connect(rotationIntervalAct, SIGNAL(triggered()), this, SLOT(intervalRotation()));
    connect(movieAct, SIGNAL(triggered()), this, SLOT(movie()));

    milxQtRenderWindow::createConnections();
}

void milxQtAnimateModel::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    createMenu(contextMenu);

    contextMenu->exec(currentEvent->globalPos());
}
