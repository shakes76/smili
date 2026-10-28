# Phase 2: VTK 9.6 & ITK 5.4 Core & Rendering Pipeline Review, Fixes & Verification

## 1. Executive Summary

Phase 2 audited the VTK 9.6 rendering pipeline, `QVTKOpenGLNativeWidget` integration, custom `vtk-ext` classes (`vtkImageViewer3`, `vtkInteractorStyleImage2`, `vtkTensor`), core `milx::*` classes (`milxFile`, `milxImage`, `milxModel`, `milxDeformableModel`), and ITK 5.4 / `itkVTKGlue` / `itk-ext` bridges. Seven critical compile/runtime crash blockers and multiple high-severity rendering and data pipeline defects were identified, remediated, and verified via the new automated test suite `test_vtk9_itk5_pipeline`.

---

## 2. Issues Identified

| ID | File(s) | Issue Description | Severity |
| :--- | :--- | :--- | :--- |
| **VP-01** | [`include/itkVTKGlue/itkVTKPolyDataToMesh.txx`](../../include/itkVTKGlue/itkVTKPolyDataToMesh.txx#L161-L195) | Used `unsigned long pointIds[3]` and `(unsigned long*)cellPoints` when calling `TriangleCellType::SetPointIds(...)`. On 64-bit Windows (`LLP64`), `unsigned long` is 32-bit whereas ITK 5.4's `PointIdentifier` (`itk::IdentifierType`) is 64-bit `unsigned long long` and VTK 9.6's `vtkIdType` is 64-bit `long long`, causing a compilation failure and 64-bit pointer truncation. | **Critical** |
| **VP-02** | [`include/Qt/milxQtRenderWindow.h`](../../include/Qt/milxQtRenderWindow.h#L650-L666) & [#L935](../../include/Qt/milxQtRenderWindow.h#L935) | Protected member `vtkSmartPointer<vtkRenderWindow> renderWindow;` shadowed `milxQtWindow::renderWindow()` and was never initialized (`nullptr`). Calling `enableCrosshair()` or `disableCrosshair()` dereferenced `this->renderWindow` (`nullptr`) and crashed immediately. | **Critical** |
| **VP-03** | [`src/Qt/milxQtPlot.cpp`](../../src/Qt/milxQtPlot.cpp#L206-L253), [`include/Qt/milxQtPlot.h`](../../include/Qt/milxQtPlot.h#L305-L309) | In `milxQtPlot::scatterPlot()`, `vtkContextView` was a local variable destroyed at function exit, `SetRenderWindow` was commented out, and `view->SetInteractor(GetRenderWindow()->GetInteractor())` stole the `QVTKInteractor` from `QVTKOpenGLNativeWidget`'s `vtkGenericOpenGLRenderWindow` and rebound it to `view`'s temporary native `vtkWin32OpenGLRenderWindow`. | **Critical** |
| **VP-04** | [`src/Qt/milxQtPlot.cpp`](../../src/Qt/milxQtPlot.cpp#L301-L308) | In `milxQtPlot::surfacePlot(vtkSmartPointer<vtkImageData> img, ...)`, `if(!imageData) imageData->DeepCopy(img);` dereferenced `imageData` when it was `nullptr`. | **Critical** |
| **VP-05** | [`src/Qt/milxQtRenderWindow.cpp`](../../src/Qt/milxQtRenderWindow.cpp#L198-L202), [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp#L1032-L1044), [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp#L1878-L1888) | `contour()` called `vtkContourWidget::New()` and immediately dereferenced `contourWidget->GetRepresentation()` (which is `nullptr` until `CreateDefaultRepresentation()` is called). | **Critical** |
| **VP-06** | [`src/Qt/milxQtUnifiedWindow.cpp`](../../src/Qt/milxQtUnifiedWindow.cpp#L284-L355) | `generateCheckerBoard()` dereferenced `checkerWidget->GetRepresentation()` (`nullptr`) without calling `checkerWidget->CreateDefaultRepresentation()`, used static `SetInputData` on an un-updated `windowLevel` filter instead of `SetInputConnection`, and checked `if(unifyImages.isEmpty() && unifyImages.size() < 2)` (allowing a single image through). | **Critical** |
| **VP-07** | [`include/Qt/milxQtAliases.h`](../../include/Qt/milxQtAliases.h#L119-L139) | `static milxQtVTKInit milxInitInstance;` in a header ran across every translation unit and plugin load, forced `QSurfaceFormat::CompatibilityProfile` over VTK 9's `QVTKOpenGLNativeWidget::defaultFormat()` (`CoreProfile`), and called `vtkOutputWindow::GetInstance()` during static init instead of `vtkObject::GlobalWarningDisplayOff()`. | **High** |
| **VP-08** | [`src/milxModel.cxx`](../../src/milxModel.cxx#L996-L1026) | In VTK 9.6, `vtkGenerateIds` (replacing deprecated `vtkIdFilter` in `Model::GenerateVertexScalars`) marks its output array as `PedigreeIds` rather than active `Scalars`, leaving `GetPointData()->GetScalars()` null. | **High** |
| **VP-09** | [`src/milxModel.cxx`](../../src/milxModel.cxx#L1159-L1186) | `Model::GenerateTensorField` passed `sphere->GetOutput()` via `SetSourceData` without calling `sphere->Update()` or `SetSourceConnection(sphere->GetOutputPort())`, resulting in empty tensor ellipsoid glyphs in VTK 9.6. | **High** |
| **VP-10** | [`src/Qt/milxQtRenderWindow.cpp`](../../src/Qt/milxQtRenderWindow.cpp#L910-L1365), [`include/Qt/milxQtPlot.h`](../../include/Qt/milxQtPlot.h#L161-L168) | All 19 `colourMapTo*()` methods called `GetDataSet()->GetScalarRange(range)` without checking `if (GetDataSet())` (which is `NULL` on base `milxQtRenderWindow` and on 3D scatter plots where `milxQtPlot::GetDataSet()` returned `imageData` instead of `model.Result()`). | **High** |
| **VP-11** | [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp#L419-L437), [#L622-L639](../../src/Qt/milxQtImage.cpp#L622-L639), [#L2637-L2644](../../src/Qt/milxQtImage.cpp#L2637-L2644), [#L3530-L3534](../../src/Qt/milxQtImage.cpp#L3530-L3534) | `getFloatData()` dereferenced `GetCharImage()` (`nullptr` on float images) instead of `GetFloatImage()`; `generateVoxelisedSurface()` left `bool ok1, ok2, ok3` uninitialized when `spacing != NULL`; `checkerBoard()` 2D `integer` branch passed `imageChar`; `resize()` `floating` branch called `refImage->GetCharImage()->GetOrigin()`. | **High** |
| **VP-12** | [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp#L519), [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp#L886), [`apps/milxLabelVisualisation/milxLabelVisualisation.cpp`](../../apps/milxLabelVisualisation/milxLabelVisualisation.cpp#L212) | Passed transposed `(height, width)` to `vtkRenderWindow::SetSize(width, height)`. | **Medium** |
| **VP-13** | [`src/vtk-ext/vtkImageViewer3.cxx`](../../src/vtk-ext/vtkImageViewer3.cxx#L58-L65) & [#L351-L356](../../src/vtk-ext/vtkImageViewer3.cxx#L351-L356), [`include/vtk-ext/vtkImageViewer3.h`](../../include/vtk-ext/vtkImageViewer3.h#L82-L89), [`include/vtk-ext/vtkInteractorStyleImage2.h`](../../include/vtk-ext/vtkInteractorStyleImage2.h#L72) | `ResetWindowLevelEvent` manually mutated `UPDATE_EXTENT` instead of calling VTK 9's `UpdateWholeExtent()`; `cbk2` leaked in `InstallPipeline()`; `SetCursor()` did not check `if (cursorActor)`; `vtkInteractorStyleImage2` left `viewer` uninitialized. | **High** |
| **VP-14** | [`src/milxFile.cxx`](../../src/milxFile.cxx#L335-L350) | `File::OpenVTKTransform` deleted `buffer` when `fin == NULL` but did not `return NULL`, falling through to `::gzgets(fin, buffer, 512)` with a null file handle and freed buffer. | **High** |
| **VP-15** | [`include/Qt/milxQtRenderWindow.h`](../../include/Qt/milxQtRenderWindow.h#L178-L198), [`src/Qt/milxQtFile.cpp`](../../src/Qt/milxQtFile.cpp#L1098), [`plugin/ssm/itkStatisticalShapeModelBase.txx`](../../plugin/ssm/itkStatisticalShapeModelBase.txx#L933) | VTK 9.7.1 removes APIs marked `VTK_DEPRECATED_IN_9_5_0` / `VTK_DEPRECATED_IN_9_6_0`: `renderer->AddActor2D(actor)` (deprecated in 9.5 in favor of `AddViewProp(actor)`), `vtkOBJExporter::SetInput(renWin)` (deprecated in 9.6 in favor of `SetRenderWindow(renWin)`), and `vtkThreshold::ThresholdByUpper` (removed in VTK 9 in favor of `SetUpperThreshold` + `SetThresholdFunction`). | **Critical (VTK 9.7+)** |

---

## 3. Remediations Applied

1. **64-bit ITK 5.4 / VTK 9.6 Mesh Bridge ([`include/itkVTKGlue/itkVTKPolyDataToMesh.txx`](../../include/itkVTKGlue/itkVTKPolyDataToMesh.txx))**:
   - Replaced `unsigned long pointIds[3]` and `(unsigned long*)cellPoints` with `typename TriangleCellType::PointIdentifier pointIds[3]` in both triangle strip and polygon conversion loops.
2. **OpenGL Surface Format & `QVTKOpenGLNativeWidget` Lifecycle ([`include/Qt/milxQtAliases.h`](../../include/Qt/milxQtAliases.h), [`include/Qt/milxQtWindow.h`](../../include/Qt/milxQtWindow.h), [`include/Qt/milxQtRenderWindow.h`](../../include/Qt/milxQtRenderWindow.h))**:
   - Guarded `milxQtVTKInit` with `static bool initialized`, preserved `QVTKOpenGLNativeWidget::defaultFormat()` (OpenGL 3.2 CoreProfile), and replaced `vtkOutputWindow::GetInstance()->GlobalWarningDisplayOff()` with `vtkObject::GlobalWarningDisplayOff()`.
   - Added `milxQtWindow::GetInteractor()` helper.
   - Removed the shadowed uninitialized `vtkSmartPointer<vtkRenderWindow> renderWindow;` member in `milxQtRenderWindow.h` and updated `enableCrosshair()` / `disableCrosshair()` to use `GetRenderWindow()`.
3. **2D/3D Plotting & `vtkContextView` ([`include/Qt/milxQtPlot.h`](../../include/Qt/milxQtPlot.h), [`src/Qt/milxQtPlot.cpp`](../../src/Qt/milxQtPlot.cpp))**:
   - Added persistent `vtkSmartPointer<vtkContextView> contextView;` member to `milxQtPlot` and bound `contextView->SetRenderWindow(GetRenderWindow())` and `contextView->SetInteractor(GetRenderWindow()->GetInteractor())` before adding charts to the scene.
   - Allocated `imageData = vtkSmartPointer<vtkImageData>::New()` when null in `milxQtPlot::surfacePlot`.
   - Updated `milxQtPlot::GetDataSet()` to return `model.Result()` for both `plotTypeSurface` and `plotType3D`.
4. **VTK 9.6 Pipeline & Widget Fixes ([`src/milxModel.cxx`](../../src/milxModel.cxx), [`src/milxFile.cxx`](../../src/milxFile.cxx), [`src/Qt/milxQtRenderWindow.cpp`](../../src/Qt/milxQtRenderWindow.cpp), [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp), [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp), [`src/Qt/milxQtUnifiedWindow.cpp`](../../src/Qt/milxQtUnifiedWindow.cpp), [`src/vtk-ext/vtkImageViewer3.cxx`](../../src/vtk-ext/vtkImageViewer3.cxx))**:
   - Explicitly set `"Point IDs"` as active scalars on `vtkGenerateIds` output in `Model::GenerateVertexScalars()`.
   - Connected `sphere->GetOutputPort()` via `SetSourceConnection` in `Model::GenerateTensorField()`.
   - Added `CreateDefaultRepresentation()` calls after `vtkContourWidget::New()` and `vtkCheckerboardWidget::New()`.
   - Connected `windowLevel` and `checkerActor` via `SetInputConnection` in `milxQtUnifiedWindow::generateCheckerBoard()` and initialized `bigBounds` in `milxQtUnifiedWindow::viewToBest()`.
   - Guarded `if (GetDataSet())` across all 19 `colourMapTo*()` methods in `milxQtRenderWindow.cpp`.
   - Fixed `getFloatData()`, `generateVoxelisedSurface()`, `checkerBoard()`, `resize()`, and `SetSize(width, height)` in `milxQtImage.cpp`, `milxQtModel.cpp`, and `milxLabelVisualisation.cpp`.

---

## 4. Verification Results

Added automated CTest target `test_vtk9_itk5_pipeline` ([`tests/test_vtk9_itk5_pipeline/test_vtk9_itk5_pipeline.cpp`](../../tests/test_vtk9_itk5_pipeline/test_vtk9_itk5_pipeline.cpp)) exercising:
1. ITK 5.4 ↔ VTK 9.6 3D float image round-trip conversion and `itk::VectorShiftScaleImageFilter` (`ScaleVectorImage`).
2. 64-bit `itk::VTKPolyDataToMesh`, `milx::Model::GenerateVertexScalars()` (active scalars in VTK 9.6), and `milx::Model::GenerateTensorField()`.
3. `vtk-ext` (`vtkImageViewer3` cursor lifecycle and `vtkInteractorStyleImage2` initialization).
4. `milxQtRenderWindow` (`QVTKOpenGLNativeWidget` render window/interactor, `enableCrosshair()` / `disableCrosshair()`, and `colourMapToJet()` null-dataset safety).
5. `milxQtImage::getFloatData()` on 3D floating-point volumes.
6. `milxQtPlot` (`scatterPlot` `vtkContextView` interactor retention, non-blocking `renameTitle` when `name` is empty, and `surfacePlot` null `imageData` initialization).

### CTest Output
```text
Test project D:/Antigravity/smili-release/build-vs2022
    Start 1: cli_milxImageApp_help
1/6 Test #1: cli_milxImageApp_help .............   Passed    0.11 sec
    Start 2: cli_milxDICOMApp_help
2/6 Test #2: cli_milxDICOMApp_help .............   Passed    0.09 sec
    Start 3: cli_milxModelApp_help
3/6 Test #3: cli_milxModelApp_help .............   Passed    0.09 sec
    Start 4: cli_milxDeformableModelApp_help
4/6 Test #4: cli_milxDeformableModelApp_help ...   Passed    0.11 sec
    Start 5: cli_milxHausdorffDistance_help
5/6 Test #5: cli_milxHausdorffDistance_help ....   Passed    0.11 sec
    Start 6: test_vtk9_itk5_pipeline
6/6 Test #6: test_vtk9_itk5_pipeline ...........   Passed    0.40 sec

100% tests passed, 0 tests failed out of 6
```
