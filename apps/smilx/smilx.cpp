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
/**
    \brief This program constructs the sMILX applicaton and executes the milxQtMain
    \author Shekhar S. Chandra, 2014
*/
//Qt
#include <QSplashScreen>

#include "milxQtMain.h"

int main(int argc, char** argv)
{
    QApplication app(argc,argv);

    QPixmap icon(":resources/smilx_icon.png");
    app.setWindowIcon(QIcon(icon));

    QPixmap pixmap(":resources/smilx_splash.png");
    QSplashScreen splash(pixmap);
//        splash.setMask(pixmap.mask());
        splash.showMessage("This software is for research purposes only and is NOT approved for clinical use", Qt::AlignBottom | Qt::AlignHCenter);
        splash.show();
    app.processEvents();
	
    milxQtMain Main;
    Main.setWindowTitle("SMILX");
    Main.show();
    splash.finish(&Main);

    // Open files if provided
    QStringList files = app.arguments();
		files.erase(files.begin()); //First element is the program name
		Main.loadFiles(files);

    return app.exec();
}
