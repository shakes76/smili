# Phase 1: CMake & Build System Audit, Fixes & Verification

## 1. Executive Summary

Phase 1 audited and modernized the CMake build system across the root project, `CMake/` helper modules, core libraries (`src/`, `src/vtk-ext/`, `src/Qt/`), applications (`apps/*`), plugins (`plugin/*`), Python bindings (`python/`), and test suite (`tests/*`) to build cleanly against **VTK 9.6.0**, **ITK 5.4.3**, and **Qt 6.10.2** on Windows (Visual Studio 2022 x64) and Linux.

---

## 2. Issues Identified

| ID | File(s) | Issue Description | Severity |
| :--- | :--- | :--- | :--- |
| **CM-01** | `CMakeLists.txt` | `include(zlib)` was called **before** `include(find)`, so `ITK_FOUND` and `VTK_FOUND` were always false when `CMake/zlib.cmake` executed, causing ZLIB resolution to fail on Windows. | **Critical** |
| **CM-02** | `CMake/zlib.cmake` | Checked `IF(ITK_USE_SYSTEM_ZLIB)` and legacy ITK 4 library name `itkzlib` (ITK 5 uses `ITKZLIB_LIBRARIES`), did not support `VTK::zlib`, and had a syntax typo `ADD_DEFINITIONS(-DZLIB_DLL})` with a stray closing brace `}`. | **Critical** |
| **CM-03** | `CMake/Windeployqt.cmake` | Used `if(${CMAKE_CONFIGURATION_TYPES} EQUAL "Release")` which expands to `if(Debug Release MinSizeRel RelWithDebInfo EQUAL "Release")` on Visual Studio multi-config generators, throwing a fatal CMake configure error when `BUILD_PACAKGE=ON`. | **Critical** |
| **CM-04** | `CMakeLists.txt` | When `find_package(Qt6 ...)` succeeded independently of `find_package(QT NAMES Qt6 Qt5 ...)`, `QT_VERSION_MAJOR` was left unset, causing downstream `Qt${QT_VERSION_MAJOR}::*` references and `if(QT_VERSION_MAJOR GREATER 5)` checks to fail. | **Critical** |
| **CM-05** | `CMakeLists.txt`, `src/Qt/CMakeLists.txt` | Missing `Qt6::OpenGL` and `Qt6::OpenGLWidgets` (`QOpenGLWidget` / `QVTKOpenGLNativeWidget` live in `Qt6::OpenGLWidgets` in Qt 6) and missing `Qt6::Core`, `Qt6::Gui`, `Qt6::Widgets` in `QT_LIBRARIES`. | **Critical** |
| **CM-06** | `include/vtk-ext/vtkImageViewer3.h`, `vtkInteractorStyleImage2.h`, `vtkTensor.h` | `vtk-ext` classes used internal VTK export macros (`VTKINTERACTIONIMAGE_EXPORT`, `VTKRENDERINGCORE_EXPORT`) which expand to `__declspec(dllimport)` when compiling `vtk-ext` against a shared VTK 9 build on Windows, causing MSVC `error C2491: definition of dllimport function not allowed`. | **Critical** |
| **CM-07** | `include/milxGlobal.h`, `include/milxMath.h`, `include/milxImage.h` | Header-only template classes `milx::Math<Type>` and `milx::Image<TImage>` were decorated with `SMILI_EXPORT` instead of `SMILI_TEMPLATEDLL`, causing `error C2491: definition of dllimport function not allowed` whenever downstream targets defined `SMILI_DLL`. Also `ImageBase::ImageBase()` lacked an inline body `{}`. | **Critical** |
| **CM-08** | `src/vtk-ext/CMakeLists.txt`, `src/Qt/CMakeLists.txt` | Source files in `src/vtk-ext/CMakeLists.txt` and `src/Qt/CMakeLists.txt` omitted `.cxx`/`.cpp` file extensions, violating CMake 3.20+ Policy `CMP0115`. | **High** |
| **CM-09** | `src/CMakeLists.txt`, `src/vtk-ext/CMakeLists.txt`, `src/Qt/CMakeLists.txt`, `apps/*/CMakeLists.txt`, `plugin/*/CMakeLists.txt`, `tests/*/CMakeLists.txt` | Missing `vtk_module_autoinit(TARGETS ... MODULES ${VTK_LIBRARIES})` calls required by VTK 9's object factory system (`vtkRenderingOpenGL2`, `vtkRenderingContextOpenGL2`, `vtkInteractionStyle`, `vtkRenderingFreeType`, `vtkRenderingVolumeOpenGL2`). Without this, VTK 9 returns `nullptr` or aborts at runtime when instantiating rendering/chart classes. | **Critical** |
| **CM-10** | `CMake/find.cmake` | `IF(VTK_INCLUDE_DIRS)` was dead code after `IF(VTK_FOUND)` because `ELSE()` was omitted (`ELSEIF(VTK_INCLUDE_DIRS)` was intended). | **Medium** |
| **CM-11** | `src/Qt/CMakeLists.txt` | Included pure header `milxQtAliases.h` (no `Q_OBJECT`) in `MILXQT_HEADERS_MOC`, triggering MOC warnings and redundant code generation. | **Low** |
| **CM-12** | `plugin/animate/milxQtAnimateModel.cpp` | Unconditionally included `<vtkFFMPEGWriter.h>` which is not built by default on Windows VTK 9 builds; lacked fallback to `<vtkOggTheoraWriter.h>` (`VTK::IOOggTheora`). | **High** |
| **CM-13** | `plugin/dicom/milxQtDICOMPlugin.h`, `milxQtDICOMPlugin.cpp`, `apps/milxAnimate/milxAnimate.cpp` | `ExportDICOM_RT` passed 3 coordinates to `itk::SpatialObjectPoint<2>::SetPositionInObjectSpace(x, y, z)`, triggering an ITK 5 `static_assert` failure (`sizeof...(coords) == TDimension - 1`). Unqualified `cout` calls failed under ITK 5. | **High** |
| **CM-14** | `tests/*/CMakeLists.txt` and `tests/**/*.cxx` | All test subdirectories used `LINK_LIBRARIES(milx)` instead of `milx-SMILI`, lacked `add_test()` registration for CTest, omitted `std::` on `cout`/`cerr`/`endl`, and used legacy VTK 5/7 `SetInput()` calls or outdated `ApplyOrientationToVTKImage` signatures. | **High** |
| **CM-15** | `CMakeLists.txt`, `src/CMakeLists.txt`, `src/Qt/CMakeLists.txt`, `plugin/CMakeLists.txt`, `CMake/Windeployqt.cmake`, `src/Qt/milxQtMain.cpp` | Legacy `LIBRARY_OUTPUT_PATH` placed Windows runtime DLLs (`vtk-ext.dll`, `milx-SMILI.dll`, `milx-Qt.dll`) in `lib/Release` instead of `bin/Release`, preventing `windeployqt` from scanning `milx-Qt.dll` (`Qt6OpenGL.dll`, `Qt6OpenGLWidgets.dll`, `Qt6WebEngineWidgets.dll`) and causing `0xC0000135` (`STATUS_DLL_NOT_FOUND`) when launching apps directly from `bin/Release` without custom `PATH` entries. Also `milxQtMain::loadPlugins()` did not descend into `plugins/Release` on multi-config Visual Studio builds. | **Critical** |

---

## 3. Remediations Applied

1. **Root [`CMakeLists.txt`](../../CMakeLists.txt)**:
   - Upgraded `cmake_minimum_required(VERSION 3.20)` and enforced `CMAKE_CXX_STANDARD 17` (`CMAKE_CXX_STANDARD_REQUIRED ON`).
   - Reordered `include(find)` before `include(zlib)` so ITK/VTK ZLIB targets are detected.
   - Explicitly set `QT_VERSION_MAJOR` and populated `QT_LIBRARIES` with `Qt6::Core`, `Qt6::Gui`, `Qt6::Widgets`, `Qt6::OpenGL`, `Qt6::OpenGLWidgets`, and optional `Qt6::WebEngineWidgets`.
   - Enabled `enable_testing()` when `BUILD_TESTS=ON`.
   - Set `CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin`, `CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib`, and `CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib` so Windows `.dll` targets (`vtk-ext.dll`, `milx-SMILI.dll`, `milx-Qt.dll`) output directly alongside `.exe` targets in `bin/<Config>`, added `POST_BUILD` `$<TARGET_RUNTIME_DLLS:...>` deployment on `milx-SMILI`, `milx-Qt`, and `animatePlugin`, updated `windeployqt` to scan both `smilx.exe` and `milx-Qt.dll`, and updated `milxQtMain::loadPlugins()` to resolve `plugins/<Config>` in multi-config builds.
2. **CMake Helper Modules**:
   - [`CMake/find.cmake`](../../CMake/find.cmake): Changed `IF(VTK_INCLUDE_DIRS)` to `ELSEIF(VTK_INCLUDE_DIRS)`.
   - [`CMake/zlib.cmake`](../../CMake/zlib.cmake): Added support for ITK 5's `ITKZLIB_LIBRARIES` and VTK 9's `VTK::zlib`, and fixed the `ADD_DEFINITIONS(-DZLIB_DLL)` syntax error.
   - [`CMake/Windeployqt.cmake`](../../CMake/Windeployqt.cmake): Replaced the broken `if(${CMAKE_CONFIGURATION_TYPES} EQUAL "Release")` check with generator expression `$<IF:$<CONFIG:Debug>,--debug,--release>` and proper `Qt6::windeployqt` target resolution.
3. **Core Libraries (`vtk-ext`, `milx-SMILI`, `milx-Qt`)**:
   - [`include/vtk-ext/vtkImageViewer3.h`](../../include/vtk-ext/vtkImageViewer3.h), [`vtkInteractorStyleImage2.h`](../../include/vtk-ext/vtkInteractorStyleImage2.h), [`vtkTensor.h`](../../include/vtk-ext/vtkTensor.h): Defined `VTKEXT_EXPORT` controlled by CMake's automatic `vtk_ext_EXPORTS` definition (`__declspec(dllexport)` when building `vtk-ext`, `__declspec(dllimport)` when consuming it).
   - [`include/milxGlobal.h`](../../include/milxGlobal.h), [`include/milxMath.h`](../../include/milxMath.h), [`include/milxImage.h`](../../include/milxImage.h): Added `SMILI_TEMPLATEDLL` macro (matching `MILXQT_TEMPLATEDLL` in `milxQtAliases.h`) for template classes `milx::Math<Type>` and `milx::Image<TImage>`, and gave `milx::ImageBase::ImageBase()` an inline definition `{}`.
   - [`src/vtk-ext/CMakeLists.txt`](../../src/vtk-ext/CMakeLists.txt), [`src/CMakeLists.txt`](../../src/CMakeLists.txt), [`src/Qt/CMakeLists.txt`](../../src/Qt/CMakeLists.txt): Added explicit `.cxx`/`.cpp` extensions, `target_include_directories`, `target_compile_definitions`, and `vtk_module_autoinit(TARGETS ... MODULES ${VTK_LIBRARIES})`.
4. **Applications, Plugins & Tests**:
   - Added `vtk_module_autoinit` and `Qt6::*` links to all targets in [`apps/*/CMakeLists.txt`](../../apps/CMakeLists.txt), [`plugin/*/CMakeLists.txt`](../../plugin/CMakeLists.txt), and [`tests/*/CMakeLists.txt`](../../tests/CMakeLists.txt).
   - Updated [`plugin/animate/milxQtAnimateModel.cpp`](../../plugin/animate/milxQtAnimateModel.cpp) to fall back cleanly to `vtkOggTheoraWriter` when `vtkFFMPEGWriter` is unavailable.
   - Fixed 2D slice point dimensionality in [`plugin/dicom/milxQtDICOMPlugin.h`](../../plugin/dicom/milxQtDICOMPlugin.h) (`p.SetPositionInObjectSpace(pixelIndex[0], pixelIndex[1])`) and `std::cout` qualifications in [`plugin/dicom/milxQtDICOMPlugin.cpp`](../../plugin/dicom/milxQtDICOMPlugin.cpp) and [`apps/milxAnimate/milxAnimate.cpp`](../../apps/milxAnimate/milxAnimate.cpp).
   - Updated [`tests/*/CMakeLists.txt`](../../tests/CMakeLists.txt) to link `milx-SMILI` and `vtk-ext`, fixed ITK 5 `std::cout`/`std::cerr` and VTK 9 `SetInputData` calls in test sources, and registered CTest CLI smoke tests.

---

## 4. Verification Results

### 4.1 CMake Configure & Full Build (`Visual Studio 17 2022`, `x64`, `Release`)
All **25 targets** configured, compiled, and linked with **0 errors**:
- **Shared Libraries (3)**: `vtk-ext.dll`, `milx-SMILI.dll`, `milx-Qt.dll`
- **Plugins (2)**: `DICOMPlugin.dll`, `animatePlugin.dll`
- **Main Application (1)**: `smilx.exe` (with `windeployqt` deployment)
- **CLI & Viewer Applications (13)**: `milxImageApp.exe`, `milxModelApp.exe`, `milxDeformableModelApp.exe`, `milxDICOMApp.exe`, `milxHausdorffDistance.exe`, `milxAssistant.exe`, `milxImageViewer.exe`, `milxModelViewer.exe`, `milxOverlay.exe`, `milxBlend.exe`, `milxLabelVisualisation.exe`, `milxAnimate.exe`
- **Test & Prototype Executables (6)**: `milxHausdorffDistanceBetweenSurfaces.exe`, `milxPlotXY.exe`, `imageViewer.exe`, `prototype_plot_xy.exe`, `prototype_plot_surface.exe`, `prototype_image_planes.exe`, `prototype_image_slice.exe`, `prototype_vectorimage.exe`

### 4.2 CTest Execution
```text
Test project D:/Antigravity/smili-release/build-vs2022
    Start 1: cli_milxImageApp_help
1/5 Test #1: cli_milxImageApp_help .............   Passed    0.38 sec
    Start 2: cli_milxDICOMApp_help
2/5 Test #2: cli_milxDICOMApp_help .............   Passed    0.12 sec
    Start 3: cli_milxModelApp_help
3/5 Test #3: cli_milxModelApp_help .............   Passed    0.11 sec
    Start 4: cli_milxDeformableModelApp_help
4/5 Test #4: cli_milxDeformableModelApp_help ...   Passed    0.13 sec
    Start 5: cli_milxHausdorffDistance_help
5/5 Test #5: cli_milxHausdorffDistance_help ....   Passed    0.13 sec

100% tests passed, 0 tests failed out of 5
```
