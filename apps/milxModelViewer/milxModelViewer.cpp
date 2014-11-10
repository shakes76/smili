/*=========================================================================
Program: SMILI
Module: milxModelViewer.cxx
Author: Shekhar Chandra
Language: C++
Created: 15 August 2012

Copyright: (c) 2012 CSIRO, Australia.

This software is protected by international copyright laws.
Any unauthorised copying, distribution or reverse engineering is prohibited.

Licence:
All rights in this Software are reserved to CSIRO. You are only permitted
to have this Software in your possession and to make use of it if you have
agreed to a Software License with CSIRO.

BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/
//Qt
#include <QApplication>
#include <QMainWindow>
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
#include "milxQtFile.h"
#include "milxQtModel.h"

int main(int argc, char *argv[])
{
    QApplication app(argc,argv);
    QMainWindow mainWindow;
  
    QPixmap icon(":resources/smilx_icon.png");
    app.setWindowIcon(QIcon(icon));

    if (argc < 2)
    {
        cerr << "milxModelViewer Application:" << endl;
        cerr << "For quick and fast display of model/surface/polydata files." << endl;
        cerr << "View configuration always matches sMILX settings wherever possible." << endl;
        cerr << "Usage:" << endl;
        cerr << "<Model Filename> " << endl;
        return EXIT_FAILURE;
    }

    std::string inputSurfaceFilename = argv[1];

    //Load the vertex table from CSV file
    milxQtModel *model = new milxQtModel(&mainWindow); //app takes ownership so need to be on stack, not heap
    milxQtFile *reader = new milxQtFile(&mainWindow);
    bool success = reader->openModel(inputSurfaceFilename.c_str(), model);

    if(!success)
    {
      cerr << "Error opening model file." << endl;
      return EXIT_FAILURE;
    }

    model->generateModel();
    
    //Setup size
    QSize desktopSize = qApp->desktop()->availableGeometry().size();
    int newWidth = 2.0*desktopSize.width()/3.0 + 0.5;
    int newHeight = 4.0*desktopSize.height()/5.0 + 0.5;
    int xOffset = (desktopSize.width()-newWidth)/2.0;
    int yOffset = (desktopSize.height()-newHeight)/2.0;
    mainWindow.resize( QSize(newWidth, newHeight) );
    mainWindow.move( QPoint(xOffset, yOffset) );
    
    //Setup view to match sMILX
    QSettings settings("Shekhar Chandra", "milxQt");
    int defaultViewMode = 2; //axial
    int defaultOrientationTypeMode = 0; //radiological
    model->setView(settings.value("defaultView", defaultViewMode).toInt());
    model->setDefaultOrientation(settings.value("defaultOrientationType", defaultOrientationTypeMode).toInt());

    //Setup menu
    QMenuBar *menuBar = new QMenuBar(&mainWindow);
    //File
    QMenu *menuFile = menuBar->addMenu("File");
    model->createMenu(menuFile);
    QAction *actionExit = menuFile->addAction("Exit");
    QObject::connect(actionExit, SIGNAL(activated()), &mainWindow, SLOT(close()));
    mainWindow.setMenuBar(menuBar);

    model->colourMapToJet();
    
    mainWindow.setWindowTitle("milxModelViewer");
    mainWindow.setCentralWidget(model);
    mainWindow.show();

    return app.exec();
}
