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
    \brief This program presents the documentation of SMILI via a simple web interface.
    \author Shekhar S. Chandra, 2014
*/
//Qt
#include <QApplication>
#include <QMainWindow>
#include <QWebEngineView>
#include <QWebEnginePage>
#include <QFile>
#include <QToolBar>
#include <QDebug>

int main(int argc, char* argv[])
{
    QApplication app(argc,argv);
	QMainWindow Main;
    Q_INIT_RESOURCE(smili);
	QWebEngineView *assistantView = new QWebEngineView(&Main);

    QFile file(":/resources/index.html");
    if(file.open(QIODevice::ReadOnly))
        assistantView->setHtml(file.readAll());
    else
        qDebug() << "Failed reading documentation.";
//    assistantView->setUrl(QUrl(":/resources/index.html"));
    assistantView->setWindowTitle("SMILI Assistant");
    Main.setCentralWidget(assistantView);
    Main.setUnifiedTitleAndToolBarOnMac(true);
    assistantView->show();

    //Connect

    //Quick setup toolbar
    QToolBar *toolBar = Main.addToolBar(QObject::tr("Navigation"));
    toolBar->addAction(assistantView->pageAction(QWebEnginePage::Back));
    toolBar->addAction(assistantView->pageAction(QWebEnginePage::Forward));
    toolBar->addAction(assistantView->pageAction(QWebEnginePage::Reload));
    toolBar->addAction(assistantView->pageAction(QWebEnginePage::Stop));

    Main.setWindowTitle("SMILI Assistant");
    Main.show();
    app.processEvents();

    return app.exec();
}
