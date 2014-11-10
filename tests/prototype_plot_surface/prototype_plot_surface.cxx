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
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkWarpScalar.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "prototype_plot_surface test:" << endl;
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

    vtkSmartPointer<vtkTable> table = csv_vert_source->GetOutput();
    cout << "Has " << table->GetNumberOfColumns() << " columns" << endl;

    vtkSmartPointer<vtkImageData> samplingArray = vtkSmartPointer<vtkImageData>::New();
        samplingArray->SetNumberOfScalarComponents(1);
        samplingArray->SetScalarTypeToDouble();
        samplingArray->SetExtent(0, table->GetNumberOfRows()-1, 0, table->GetNumberOfColumns()-1, 0, 0);
        samplingArray->Update();

    //Copy data into imagedata
    for(int j = 0; j < table->GetNumberOfRows(); j ++)
        for(int k = 0; k < table->GetNumberOfColumns(); k ++)
            samplingArray->SetScalarComponentFromDouble(j, k, 0, 0, table->GetValue(j, k).ToDouble());
    samplingArray->AllocateScalars();
    samplingArray->Update();

    cout << "Warping Scalars" << endl;
    vtkSmartPointer<vtkImageDataGeometryFilter> geometry = vtkSmartPointer<vtkImageDataGeometryFilter>::New();
        geometry->SetInput(samplingArray);
        geometry->Update();

    double range[2], bounds[6];
    samplingArray->GetScalarRange(range);
    samplingArray->GetBounds(bounds);
    const double scaling = 0.25*( 1.0 / ( 2.0*(range[1]-range[0]) / (bounds[1]+bounds[3]) ) );

    vtkSmartPointer<vtkWarpScalar> warpScalar = vtkSmartPointer<vtkWarpScalar>::New();
        warpScalar->SetInputConnection(geometry->GetOutputPort());
//        warpScalar->SetInput(samplingArray);
//        warpScalar->XYPlaneOn();
//        warpScalar->UseNormalOff();
        warpScalar->SetNormal(0, 0, 1);
        warpScalar->UseNormalOn();
        warpScalar->SetScaleFactor(scaling);
        warpScalar->Update();

      // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(warpScalar->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    // Visualize
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    renderer->SetBackground(1,1,1); // Background color white

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
