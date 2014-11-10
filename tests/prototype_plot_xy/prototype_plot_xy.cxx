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
///Other Headers
#include <string>
///VTK
#include <vtkSmartPointer.h>
#include <vtkDelimitedTextReader.h>
#include <vtkTable.h>
#include <vtkChartXY.h>
#include <vtkContextView.h>
#include <vtkPlotPoints.h>
#include <vtkContextScene.h>

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "prototype_plot_xy test:" << endl;
        cerr << "CSV Filename " << endl;
        return EXIT_FAILURE;
    }

    std::string inputSurfaceFilename = argv[1];
    //~ std::string inputImageFilename1 = argv[2];
    
    //Load the vertex table from CSV file
    cout << "Reading file: " << inputSurfaceFilename << endl;
    vtkSmartPointer<vtkDelimitedTextReader> csv_vert_source = vtkSmartPointer<vtkDelimitedTextReader>::New();
        csv_vert_source->SetFieldDelimiterCharacters(",");
        csv_vert_source->DetectNumericColumnsOn();
        //~ csv_vert_source->SetHaveHeaders(True);
        csv_vert_source->SetFileName(inputSurfaceFilename.c_str());
        csv_vert_source->Update();
        cout << "Has " << csv_vert_source->GetOutput()->GetNumberOfColumns() << " columns" << endl;
    
    vtkSmartPointer<vtkChartXY> templateChart = vtkSmartPointer<vtkChartXY>::New();
        templateChart->SetShowLegend(true);
    
    vtkPlot *templatePoints = templateChart->AddPlot(vtkChart::POINTS);
        templatePoints->SetInput(csv_vert_source->GetOutput(), 0, 1);
        templatePoints->SetColor(255, 0, 0, 255);
        templatePoints->SetWidth(1.0);
        vtkPlotPoints::SafeDownCast(templatePoints)->SetMarkerStyle(vtkPlotPoints::CROSS);
        //~ vtkPlotPoints::SafeDownCast(templatePoints)->SetMarkerStyle(vtkPlotPoints::PLUS);
        
    cout << "Plotting... " << endl;
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
        view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
        view->GetRenderWindow()->SetSize(400, 300);
        view->GetScene()->AddItem(templateChart);
        view->GetRenderWindow()->SetMultiSamples(0);
        view->GetInteractor()->Initialize();
        view->GetInteractor()->Start();

    return EXIT_SUCCESS;
}
