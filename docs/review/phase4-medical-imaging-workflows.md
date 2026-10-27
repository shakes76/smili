# Phase 4: Biomedical Imaging & Scientific Visualization Assessment, Fixes & Verification

**Reviewers**: `medical-imaging-researcher` (Medical Imaging Researcher & SMILI Architect), `senior-vtk-developer`, & `senior-qt-developer`  
**Target Toolchain**: ITK 5.4.3, VTK 9.6.0, Qt 6.10.2 (`msvc2022_64`)  
**Branch**: `topic-review`

---

## 1. Executive Summary

Phase 4 conducted a domain-specific biomedical imaging and scientific visualization review of SMILI and `sMILX` against the architecture described in [Chandra et al., *Computer Methods and Programs in Biomedicine*, 2018](../../papers/CMPB_Chandra2018.pdf). The review covered:
1. **Dual-Layer API Consistency** (`milx::File`, `milx::Image<T>`, `milx::Model`, `milx::DeformableModel` vs. `milxQtFile`, `milxQtImage`, `milxQtModel`, `milxQtUnifiedWindow`).
2. **Physical Coordinate Space, Anatomical Orientation & Reslicing** (Origin, Spacing, Direction Cosines, Radiological vs. Neurological conventions).
3. **2D, 3D, and 4D Multi-Component / Vector / Tensor / RGB Volumes**.
4. **Surface Meshes, Point Clouds, Vector/Tensor Fields & Deformable Models**.
5. **Multi-Window Overlays (`importFrom` / `customOperation`) & Measurement Widgets**.
6. **VNL (`vnl_vector` / `vnl_matrix`) & `std::vector` / NumPy Array Interoperability**.

All critical and high-severity domain defects were remediated and verified via a new automated CTest suite (`test_smili_imaging_workflows`).

---

## 2. Issues Identified & Remediated

### 2.1 Physical Coordinate Space, Anatomical Orientation & Resampling (`include/milxImage.h`, `include/milxDeformableModel.h`)
- **Severity**: Critical (Spatial Misregistration in Clinical Volumes & Meshes)
- **Files**:
  - [`include/milxImage.h`](../../include/milxImage.h)
  - [`include/milxDeformableModel.h`](../../include/milxDeformableModel.h)
- **Root Causes & Fixes**:
  1. **`milx::Image<T>::ResizeImage` Resetting Origin & Direction Cosines**: Both 2-argument and 3-argument overloads of `ResizeImage()` set `Size` and `OutputSpacing` on `itk::ResampleImageFilter` but omitted `SetOutputOrigin(img->GetOrigin())` and `SetOutputDirection(img->GetDirection())`, resetting any resampled NIfTI/DICOM scan to origin `(0, 0, 0)` and identity direction cosines. Also fixed `itk::IdentityTransform<double, imgDimension>` to use `TImage::ImageDimension` so 2D images are supported.
  2. **`milx::Image<T>::ApplyOrientationToVTKImage` Origin Pivot Order & 2D Bounds**: `transform->Concatenate(transformMatrix)` (Y-flip) was applied **before** `transform->Translate(-origin)` under `PostMultiply()`, shifting volumes with non-zero $o_y$ by $2 o_y$. Moved `Translate(-origin)` before `Concatenate(transformMatrix)` and guarded the $3 \times 3$ direction loop and Z-origin access with `TImage::ImageDimension` bounds checks.
  3. **`milx::DeformableModel::ApplyOrientation` Missing Centroid Restoration**: When `applyOrigin == false`, `ApplyOrientation()` translated the mesh by `-centroid` before rotating by `direction`, but never translated back by `+centroid` (and when `applyOrigin == true`, pivoted about `-centroid` instead of `-origin`). Fixed to pivot about `origin` when `applyOrigin == true` and about `centroid` (restoring `+centroid`) when `applyOrigin == false`.

### 2.2 4D Time-Series & 6/9-Component DTI Tensor Heap Buffer Overflow (`include/milxImage.h`)
- **Severity**: Critical (Heap Corruption on $>3$ Components)
- **File**: [`include/milxImage.h`](../../include/milxImage.h)
- **Root Cause**: `milx::Image<T>::ConvertITKVectorImageToVTKImage()` hardcoded `field->AllocateScalars(VTK_FLOAT, 3);` under `#if VTK_MAJOR_VERSION > 5`, and then looped `for (vtkIdType j = 0; j < img->GetNumberOfComponentsPerPixel(); j++)` calling `SetScalarComponentFromFloat(k, l, m, j, ...)`. Loading any 6- or 9-component DTI tensor volume or 4D time-series ($T > 3$) caused an immediate heap buffer overflow.
- **Fix**: Changed allocation to `field->AllocateScalars(VTK_FLOAT, img->GetNumberOfComponentsPerPixel());`. Also fixed the off-by-one max component index (`GetNumberOfComponentsPerPixel() - 1`) in `milxQtImage::component()`.

### 2.3 1-Channel vs. 3-Channel `unsigned char` Image Export & 32-Bit Integer Image Saving (`src/Qt/milxQtFile.cpp`, `include/Qt/milxQtImage.h`)
- **Severity**: Critical (3x Buffer Over-Read on 8-Bit Grayscale Segmentation Masks)
- **Files**:
  - [`src/Qt/milxQtFile.cpp`](../../src/Qt/milxQtFile.cpp)
  - [`include/Qt/milxQtImage.h`](../../include/Qt/milxQtImage.h)
- **Root Cause**:
  1. `milxQtFile::saveImage(const QString, vtkSmartPointer<vtkImageData>)` checked `typeStr == "unsigned char"` (`integerFormat = true`) and unconditionally converted the VTK image using `itk::Image<rgbPixelType, ...>` (3 bytes/pixel), reading 3x past the buffer end when saving a 1-channel 8-bit label map.
  2. `milxQtFile::saveImage(const QString, milxQtImage*)` lacked an `else if (data->is32BitImage())` branch for 32-bit `intImageType` volumes, falling through to `GetFloatImage()` (null), and `milxQtImage::isFloatingPointImage()` did not check `!integer`.
- **Fix**: Branched on `data->GetNumberOfScalarComponents() >= 3` vs. `1` (`ucharImageSliceType` / `charImageType`), added the `data->is32BitImage()` branch calling `WriteImageUsingITK<intImageType>(..., data->GetIntImage())`, and updated `isFloatingPointImage()` to `(!eightbit && !integer && !rgb && !vectorised)`.

### 2.4 VNL & `std::vector` / NumPy Data Exchange Use-After-Free & Tuple Count Bugs (`include/milxImage.h`, `src/Qt/milxQtModel.cpp`, `src/Qt/milxQtImage.cpp`)
- **Severity**: Critical (Use-After-Free & 3x Buffer Over-Read in Python/C++ Array Exchange)
- **Files**:
  - [`include/milxImage.h`](../../include/milxImage.h)
  - [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp)
  - [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp)
- **Root Causes & Fixes**:
  1. **`milx::Image<T>::ImportVectorToImage` & `ImportMatrixToImage`**:
     - Default argument `itk::SmartPointer<TImage> image = NULL` caused a C++17 `C2440` overload ambiguity in ITK 5.4 (`SmartPointer(std::nullptr_t)` vs `SmartPointer(ObjectType*)`); updated default arguments on both `ImportVectorToImage` and `ImportMatrixToImage` to `= nullptr`.
     - Used `importImageFilterWillOwnTheBuffer = false` without duplicating the output image, leaving the returned ITK image with a dangling pointer as soon as the caller's `vnl_vector` or `vnl_matrix` went out of scope. Also wrote `size[2] = 1` unconditionally even when `TImage::ImageDimension == 2`. Fixed by guarding `if (TImage::ImageDimension > 2) size[2] = 1;` and returning `DuplicateImage(importFilter->GetOutput())`.
  2. **`milxQtModel::SetInput(const std::vector<double>&)`**: Passed `numElements` (`inputPoints.size()`) instead of `numElements / 3` to `pointsArray->SetNumberOfTuples()`, creating 3x too many 3D points and reading 3x past the vector, and used `SetVoidArray(..., save=1)` which left a dangling pointer when a temporary `std::vector<double>` (such as from Shiboken6/Python) was destroyed. Fixed to `SetNumberOfTuples(numElements / numComponents)` and copied elements into `pointsArray`.
  3. **`milxQtModel::getData()`**: Downcast `GetPoints()->GetData()` to `vtkDoubleArray*` without null checking, crashing on standard `VTK_FLOAT` meshes. Rewrote using `GetPoints()->GetPoint(ptIdx, ...)` so both `float` and `double` meshes are supported.
  4. **`milxQtImage::setData(const std::vector<float>&, ...)`**: Used `SetVoidArray` (dangling pointer on temporary vectors) and left `imageFloat` null so `getFloatData()` failed. Fixed by copying values into `scalars`, populating `imageFloat = milx::Image<floatImageType>::ConvertVTKImageToITKImage(image)`, and adding a fallback in `getFloatData()`.

### 2.5 Surface Mesh & Collection Processing Bugs (`include/milxModel.h`, `src/milxModel.cxx`)
- **Severity**: Critical
- **Files**:
  - [`include/milxModel.h`](../../include/milxModel.h)
  - [`src/milxModel.cxx`](../../src/milxModel.cxx)
- **Fixes**:
  1. `Model::ClearArrays`: Changed `for (int j = 0; j < GetNumberOfArrays(); j++) RemoveArray(0);` (which only removed half the arrays because `GetNumberOfArrays()` decrements on each removal) to `while (mesh->GetPointData()->GetNumberOfArrays() > 0) mesh->GetPointData()->RemoveArray(0);`.
  2. `Model::MaskScalars`: Fixed duplicate check `!model1->GetPointData()->GetScalars()` to check `!model2->GetPointData()->GetScalars()`.
  3. `Model::ScalarThresholdCollection`: Fixed loop to process all items in `collection` starting from item 0 (previously started at `j = 1`, skipping item 0).
  4. `Model::ClipCollection`: Fixed reversed argument order `Clip(belowVal, aboveVal)` (was `Clip(aboveVal, belowVal)`).
  5. `Model::MeanSquaredErrorCollection`: Fixed traversal to accumulate MSE across all models in the collection against `surface1`.
  6. `Model::ProcrustesAlignCollection`: Moved `procrustes->Update()` **before** `procrustes->GetOutput()->GetBlock(0)` to fix a guaranteed null-pointer dereference in VTK 6+.
  7. `Model::IterativeClosestPointsAlignCollection`: Added the reference surface `last` to `alignedCollection`.

### 2.6 Multi-Window Overlays (`Import View From`), Camera Views & Measurement Widgets (`src/Qt/milxQtMain.cpp`, `src/Qt/milxQtImage.cpp`, `src/Qt/milxQtModel.cpp`, `src/Qt/milxQtRenderWindow.cpp`, `src/vtk-ext/vtkImageViewer3.cxx`)
- **Severity**: High
- **Fixes**:
  1. **`Import View From` (`customOperation`)**: Fixed `milxQtMain::windowActionList()` to parent each `QAction` to `win` (`milxQtWindow*`) instead of `currentWindow` (`QMdiSubWindow*`), and added child widget resolution in `milxQtImage::customOperation()` and `milxQtRenderWindow::customOperation()`.
  2. **Camera Positioning in `milxQtModel`**: Fixed `viewToXYPlane()`, `viewToZXPlane()`, and `viewToZYPlane()` to offset the camera position from `modelCentroid` rather than `(0, 0, 0)`, preventing oblique views when `modelCentroid != (0, 0, 0)`. Also fixed `defaultViewMode = 0` (`AXIAL`) in `milxQtMain::readSettings()` and `resetSettings()`, and fixed RGB channel order in `milxQtModel::generateEdges()` and `milxQtMain::display(vtkPolyDataCollection*)`.
  3. **Measurement Widgets & Crosshair Cursor**: Fixed operator precedence bug `(imageBounds[1] - imageBounds[0]) / 10.0` (was `imageBounds[1] - imageBounds[0] / 10`) and used `imageData->GetCenter(imgCenter)` in `milxQtImage::SetupWidgets()` and `vtkImageViewer3::EnableCursor()`.
  4. **Voxelisation Bounds & Interpolation Flags**: Restored `bounds[6]` after padding in `milxQtImage::generateVoxelisedSurface()` so caller arrays (`bigBounds` in `milxQtUnifiedWindow::generateDifference`) are not mutated across iterations, and updated `enableInterpolateDisplay()` / `disableInterpolateDisplay()` to synchronize both `interpolateAct` and `viewer->GetImageActor()` whether called before or after `generateImage()`.

---

## 3. Automated Verification (CTest)

A dedicated biomedical imaging and visualization regression test target (`test_smili_imaging_workflows`) was created in [`tests/test_smili_imaging_workflows/test_smili_imaging_workflows.cpp`](../../tests/test_smili_imaging_workflows/test_smili_imaging_workflows.cpp) and registered in [`tests/CMakeLists.txt`](../../tests/CMakeLists.txt).

### 3.1 Sub-Tests Verified in `test_smili_imaging_workflows`
1. **Test 1**: Physical coordinate preservation (non-zero origin `(-15.5, 24.25, -8.75)`, anisotropic spacing, and oblique rotation direction matrix) across `milx::Image<T>::ResizeImage()` and `milx::Image<T>::ApplyOrientationToVTKImage()`.
2. **Test 2**: 6-component DTI / 4D multi-component `itk::VectorImage<float, 3>` conversion via `milx::Image<VectorImageType>::ConvertITKVectorImageToVTKImage()` without heap corruption.
3. **Test 3**: `vnl_vector<float>` and `vnl_matrix<float>` import via `ImportVectorToImage()` and `ImportMatrixToImage()` verifying deep-copy buffer ownership after caller VNL objects are destroyed.
4. **Test 4**: `milx::Model` (`ClearArrays` removing all 4 arrays, `MaskScalars` null guard, `ProcrustesAlignCollection` on a 3-mesh `vtkPolyDataCollection`, and `milx::DeformableModel::ApplyOrientation` preserving centroid).
5. **Test 5**: `std::vector` / NumPy array bridges in `milxQtModel` (`SetInput` and `getData` on both `double` and `float` meshes) and `milxQtImage` (`setData` and `getFloatData` after temporary vector destruction).
6. **Test 6**: `milxQtImage::generateVoxelisedSurface()` preserving caller `bounds[6]` and `enableInterpolateDisplay()` / `disableInterpolateDisplay()` synchronization.

### 3.2 CTest Execution Output

```text
Test project D:/Antigravity/smili-release/build-vs2022
    Start 1: cli_milxImageApp_help
1/8 Test #1: cli_milxImageApp_help .............   Passed    0.11 sec
    Start 2: cli_milxDICOMApp_help
2/8 Test #2: cli_milxDICOMApp_help .............   Passed    0.11 sec
    Start 3: cli_milxModelApp_help
3/8 Test #3: cli_milxModelApp_help .............   Passed    0.11 sec
    Start 4: cli_milxDeformableModelApp_help
4/8 Test #4: cli_milxDeformableModelApp_help ...   Passed    0.11 sec
    Start 5: cli_milxHausdorffDistance_help
5/8 Test #5: cli_milxHausdorffDistance_help ....   Passed    0.11 sec
    Start 6: test_vtk9_itk5_pipeline
6/8 Test #6: test_vtk9_itk5_pipeline ...........   Passed    0.39 sec
    Start 7: test_qt6_smilx_gui
7/8 Test #7: test_qt6_smilx_gui ................   Passed    0.87 sec
    Start 8: test_smili_imaging_workflows
8/8 Test #8: test_smili_imaging_workflows ......   Passed    0.29 sec

100% tests passed, 0 tests failed out of 8

Total Test time (real) =   2.10 sec
```
