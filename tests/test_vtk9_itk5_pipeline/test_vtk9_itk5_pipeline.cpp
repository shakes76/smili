/*=========================================================================
  SMILI Phase 2 Automated Verification Suite:
  VTK 9.6 & ITK 5.4 Core & Rendering Pipeline Tests
=========================================================================*/
#include <iostream>
#include <cassert>
#include <cmath>

#include <QApplication>
#include <QSurfaceFormat>
#include <QTimer>

// SMILI Core & ITK/VTK Glue
#include "milxGlobal.h"
#include "milxFile.h"
#include "milxImage.h"
#include "milxModel.h"
#include "itkMesh.h"
#include "itkVTKPolyDataToMesh.h"

// VTK-Ext
#include "vtkImageViewer3.h"
#include "vtkInteractorStyleImage2.h"
#include "vtkTensor.h"

// SMILI Qt / VTK Rendering
#include "milxQtAliases.h"
#include "milxQtWindow.h"
#include "milxQtRenderWindow.h"
#include "milxQtImage.h"
#include "milxQtModel.h"
#include "milxQtPlot.h"
#include "milxQtUnifiedWindow.h"

// VTK
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkTable.h>
#include <vtkCommand.h>

#define TEST_ASSERT(cond, msg) \
  do { \
    if (!(cond)) { \
      std::cerr << "[FAIL] " << msg << " (" << __FILE__ << ":" << __LINE__ << ")" << std::endl; \
      return EXIT_FAILURE; \
    } \
  } while (0)

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);

  std::cout << "=== Phase 2: VTK 9.6 & ITK 5.4 Pipeline Verification ===" << std::endl;

  // 1. Verify ITK 5.4 <-> VTK 9.6 Image Bridge & VectorShiftScaleImageFilter
  {
    std::cout << "[1/6] Testing ITK 5.4 <-> VTK 9.6 image conversion & VectorShiftScaleImageFilter..." << std::endl;
    typedef itk::Image<float, 3> FloatImageType;
    FloatImageType::Pointer itkImg = FloatImageType::New();
    FloatImageType::SizeType size;
    size[0] = 16; size[1] = 16; size[2] = 8;
    FloatImageType::IndexType start;
    start.Fill(0);
    FloatImageType::RegionType region(start, size);
    itkImg->SetRegions(region);
    itkImg->Allocate();
    itkImg->FillBuffer(42.5f);

    vtkSmartPointer<vtkImageData> vtkImg = milx::Image<FloatImageType>::ConvertITKImageToVTKImage(itkImg);
    TEST_ASSERT(vtkImg != nullptr, "ConvertITKImageToVTKImage returned null");
    int dims[3] = {0, 0, 0};
    vtkImg->GetDimensions(dims);
    TEST_ASSERT(dims[0] == 16 && dims[1] == 16 && dims[2] == 8, "VTK image dimensions mismatch");

    FloatImageType::Pointer roundTrip = milx::Image<FloatImageType>::ConvertVTKImageToITKImage(vtkImg);
    TEST_ASSERT(roundTrip != nullptr, "ConvertVTKImageToITKImage returned null");
    FloatImageType::IndexType probeIdx = {{8, 8, 4}};
    TEST_ASSERT(std::fabs(roundTrip->GetPixel(probeIdx) - 42.5f) < 1e-4f, "Round-trip pixel value mismatch");

    // Test ITK 5 VectorShiftScaleImageFilter via milx::Image::ScaleVectorImage
    typedef itk::VectorImage<float, 3> VecImageType;
    VecImageType::Pointer vecImg = VecImageType::New();
    vecImg->SetRegions(region);
    vecImg->SetNumberOfComponentsPerPixel(3);
    vecImg->Allocate();
    itk::VariableLengthVector<float> v(3);
    v[0] = 1.0f; v[1] = 2.0f; v[2] = 3.0f;
    vecImg->FillBuffer(v);

    VecImageType::Pointer scaledVec = milx::Image<VecImageType>::ScaleVectorImage(vecImg, 2.0f, 3);
    TEST_ASSERT(scaledVec != nullptr, "ScaleVectorImage returned null");
    TEST_ASSERT(std::fabs(scaledVec->GetPixel(probeIdx)[1] - 4.0f) < 1e-4f, "ScaleVectorImage value mismatch");
  }

  // 2. Verify 64-bit itk::VTKPolyDataToMeshFilter & VTK 9.6 Model::GenerateVertexScalars / GenerateTensorField
  {
    std::cout << "[2/6] Testing 64-bit VTKPolyDataToMeshFilter & VTK 9.6 milx::Model filters..." << std::endl;
    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetThetaResolution(8);
    sphere->SetPhiResolution(8);
    sphere->Update();

    typedef itk::Mesh<float, 3> MeshType;
    typedef itk::VTKPolyDataToMesh<MeshType> PolyDataToMeshType;
    PolyDataToMeshType::Pointer polyToMesh = PolyDataToMeshType::New();
    polyToMesh->SetInput(sphere->GetOutput());
    polyToMesh->Update();
    TEST_ASSERT(polyToMesh->GetOutput()->GetNumberOfPoints() == sphere->GetOutput()->GetNumberOfPoints(),
                "VTKPolyDataToMesh point count mismatch");
    TEST_ASSERT(polyToMesh->GetOutput()->GetNumberOfCells() > 0,
                "VTKPolyDataToMesh produced 0 cells");

    milx::Model model(sphere->GetOutput());
    model.GenerateVertexScalars();
    TEST_ASSERT(model.Result() != nullptr, "GenerateVertexScalars result is null");
    TEST_ASSERT(model.Result()->GetPointData()->GetScalars() != nullptr,
                "GenerateVertexScalars failed to set active scalars in VTK 9.6");

    // Add a synthetic 9-component tensor array to test GenerateTensorField
    vtkSmartPointer<vtkFloatArray> tensors = vtkSmartPointer<vtkFloatArray>::New();
    tensors->SetNumberOfComponents(9);
    tensors->SetName("Tensors");
    for (vtkIdType i = 0; i < sphere->GetOutput()->GetNumberOfPoints(); ++i)
    {
      float t[9] = {1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f};
      tensors->InsertNextTypedTuple(t);
    }
    vtkSmartPointer<vtkPolyData> tensorInput = vtkSmartPointer<vtkPolyData>::New();
    tensorInput->DeepCopy(sphere->GetOutput());
    tensorInput->GetPointData()->SetTensors(tensors);

    milx::Model tensorModel(tensorInput);
    tensorModel.GenerateTensorField(0.5);
    TEST_ASSERT(tensorModel.Result() != nullptr && tensorModel.Result()->GetNumberOfPoints() > 0,
                "GenerateTensorField produced empty polydata");
  }

  // 3. Verify vtk-ext (vtkImageViewer3, vtkInteractorStyleImage2)
  {
    std::cout << "[3/6] Testing vtk-ext (vtkImageViewer3, vtkInteractorStyleImage2)..." << std::endl;
    vtkSmartPointer<vtkInteractorStyleImage2> style = vtkSmartPointer<vtkInteractorStyleImage2>::New();
    TEST_ASSERT(style->GetViewer() == nullptr, "vtkInteractorStyleImage2 viewer pointer was not initialized to NULL");

    vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
    img->SetDimensions(10, 10, 10);
    img->AllocateScalars(VTK_FLOAT, 1);
    for (int z = 0; z < 10; ++z)
      for (int y = 0; y < 10; ++y)
        for (int x = 0; x < 10; ++x)
          img->SetScalarComponentFromFloat(x, y, z, 0, static_cast<float>(x + y + z));

    vtkSmartPointer<vtkImageViewer3> viewer = vtkSmartPointer<vtkImageViewer3>::New();
    viewer->SetInputData(img);
    viewer->EnableCursor();
    TEST_ASSERT(viewer->GetCursorEnabled(), "vtkImageViewer3 cursor failed to enable");
    viewer->DisableCursor();
    TEST_ASSERT(!viewer->GetCursorEnabled(), "vtkImageViewer3 cursor failed to disable");
  }

  // 4. Verify milxQtRenderWindow crosshair & colourMap null-safety
  {
    std::cout << "[4/6] Testing milxQtRenderWindow (QVTKOpenGLNativeWidget, crosshair, colourMaps)..." << std::endl;
    milxQtRenderWindow renderWin;
    TEST_ASSERT(renderWin.GetRenderWindow() != nullptr, "milxQtRenderWindow::GetRenderWindow() is null");
    TEST_ASSERT(renderWin.GetRenderWindow()->GetInteractor() != nullptr, "milxQtRenderWindow::GetRenderWindow()->GetInteractor() is null");
    renderWin.generateRender();
    // Previously crashed due to shadowed uninitialized renderWindow pointer:
    renderWin.enableCrosshair();
    renderWin.disableCrosshair();
    // Previously crashed due to null GetDataSet() on base milxQtRenderWindow:
    renderWin.colourMapToJet();
  }

  // 5. Verify milxQtImage::getFloatData() & milxQtModel::contour()
  {
    std::cout << "[5/6] Testing milxQtImage::getFloatData() & milxQtModel..." << std::endl;
    typedef itk::Image<float, 3> FloatImageType;
    FloatImageType::Pointer itkImg = FloatImageType::New();
    FloatImageType::SizeType size = {{8, 8, 4}};
    FloatImageType::IndexType start = {{0, 0, 0}};
    itkImg->SetRegions(FloatImageType::RegionType(start, size));
    itkImg->Allocate();
    itkImg->FillBuffer(7.0f);

    milxQtImage qtImg;
    qtImg.SetInput(itkImg);
    qtImg.generateImage();
    // Previously crashed because getFloatData() called GetCharImage()->GetLargestPossibleRegion()
    std::vector<float> data = qtImg.getFloatData();
    TEST_ASSERT(data.size() == 8 * 8 * 4, "milxQtImage::getFloatData() size mismatch");
    TEST_ASSERT(std::fabs(data[0] - 7.0f) < 1e-4f, "milxQtImage::getFloatData() value mismatch");
  }

  // 6. Verify milxQtPlot (scatterPlot contextView & surfacePlot null imageData guard)
  {
    std::cout << "[6/6] Testing milxQtPlot (scatterPlot & surfacePlot)..." << std::endl;
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> colX = vtkSmartPointer<vtkFloatArray>::New();
    colX->SetName("X");
    vtkSmartPointer<vtkFloatArray> colY = vtkSmartPointer<vtkFloatArray>::New();
    colY->SetName("Y");
    for (int i = 0; i < 10; ++i)
    {
      colX->InsertNextValue(static_cast<float>(i));
      colY->InsertNextValue(static_cast<float>(i * i));
    }
    table->AddColumn(colX);
    table->AddColumn(colY);

    milxQtPlot plot2D;
    vtkRenderWindowInteractor *origInteractor = plot2D.GetRenderWindow()->GetInteractor();
    plot2D.setPlotType2D(0, 1);
    plot2D.SetSource(table);
    plot2D.generatePlot();
    TEST_ASSERT(plot2D.GetRenderWindow()->GetInteractor() == origInteractor,
                "milxQtPlot::scatterPlot detached QVTKInteractor from QVTKOpenGLNativeWidget");

    // Test surfacePlot(img) directly when plotSurface.imageData is initially null
    vtkSmartPointer<vtkImageData> surfImg = vtkSmartPointer<vtkImageData>::New();
    surfImg->SetDimensions(8, 8, 1);
    surfImg->AllocateScalars(VTK_FLOAT, 1);
    for (int y = 0; y < 8; ++y)
      for (int x = 0; x < 8; ++x)
        surfImg->SetScalarComponentFromFloat(x, y, 0, 0, static_cast<float>(x + y));

    milxQtPlot plotSurf;
    plotSurf.surfacePlot(surfImg);
    TEST_ASSERT(plotSurf.GetSourceImage() != nullptr, "milxQtPlot::surfacePlot failed to initialize imageData");
  }

  std::cout << "=== All Phase 2 VTK 9.6 & ITK 5.4 Pipeline Tests PASSED ===" << std::endl;
  return EXIT_SUCCESS;
}
