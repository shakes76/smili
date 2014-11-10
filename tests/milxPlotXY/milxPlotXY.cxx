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
//Qt
#include <QApplication>
#include <QMainWindow>

#include "milxFile.h"
#include "milxQtPlot.h"

int main(int argc, char *argv[])
{
    QApplication app(argc,argv);
    QMainWindow mainWindow;

    if (argc < 2)
    {
        cerr << "milxPlotXY Test:" << endl;
        cerr << "CSV Filename " << endl;
        return EXIT_FAILURE;
    }

    std::string inputSurfaceFilename = argv[1];

    //Load the vertex table from CSV file
    vtkSmartPointer<vtkTable> table;
    bool success = milx::File::OpenDelimitedText(inputSurfaceFilename, table);

    if(!success)
    {
      cerr << "Error opening CSV file." << endl;
      return EXIT_FAILURE;
    }
    cout << "Has " << table->GetNumberOfRows() << " points" << endl;

    milxQtPlot *plot = new milxQtPlot(&mainWindow); //hierachical deletion
    plot->scatterPlot(table, 0, 1);

    mainWindow.setCentralWidget(plot);
    mainWindow.show();

    return app.exec();
}
