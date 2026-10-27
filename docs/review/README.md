# SMILI & sMILX: VTK 9.6 / Qt 6.10 / ITK 5.4 Porting Review & Remediation

This directory contains the phase-by-phase technical audit, root-cause analysis, code remediation, and automated verification records for the **Simple Medical Imaging Library Interface (SMILI)** and its flagship multi-window application **`sMILX`** ([Chandra et al., *Computer Methods and Programs in Biomedicine*, 2018](../../papers/CMPB_Chandra2018.pdf)).

## Review Team (Specialist Agents)

1. **Senior VTK/Kitware & ITK Developer (`senior-vtk-developer`)** — VTK 9.6 module system (`VTK::*`), `QVTKOpenGLNativeWidget` & OpenGL 3.2 Core Profile lifecycle, `vtk_module_autoinit`, `vtk-ext` custom classes, ITK 5.4 dynamic threading, and `itkVTKGlue` pipeline synchronization.
2. **Senior Qt 6 GUI & CMake Developer (`senior-qt-developer`)** — Modern CMake (3.20+) target architecture, `Qt6::*` (`Core`, `Gui`, `Widgets`, `OpenGL`, `OpenGLWidgets`, `WebEngineWidgets`), `QMdiArea` + `QVTKOpenGLNativeWidget` window management, Qt 6 signal/slot modernization, and `sMILX` plugin stability.
3. **Medical Imaging Researcher & SMILI Architect (`medical-imaging-researcher`)** — Dual-layer API consistency (`milx::*` ITK/VTK core vs. `milxQt*` Qt GUI wrappers), physical coordinate space and anatomical orientation (Radiological vs. Neurological, direction cosines, origin/spacing), 2D/3D/4D volumes, surface meshes, vector/tensor fields, overlays (`importFrom`), interactive widgets, and NumPy/VNL interoperability.

## Toolchain Environment

| Component | Version | Local Build / Install Path |
| :--- | :--- | :--- |
| **VTK** | 9.6.0 | `D:\Dev\VTK-9.6.0\build-vs2022` |
| **ITK** | 5.4.3 | `D:\Dev\InsightToolkit-5.4.3\build-vs2022` |
| **Qt** | 6.10.2 (MSVC 2022 64-bit) | `D:\Qt\6.10.2\msvc2022_64` |
| **Compiler / Generator** | MSVC 19.44 (`Visual Studio 17 2022`, `x64`) | C++17 |

## Phase Documentation Index

- **[Phase 1: CMake & Build System Audit, Fixes & Verification](phase1-cmake-audit.md)** *(Completed — `b971431`)*
- **[Phase 2: VTK 9.6 & ITK 5.4 Core & Rendering Pipeline Review, Fixes & Verification](phase2-vtk9-itk5-pipeline.md)** *(In Progress)*
- **Phase 3: Qt 6 GUI Architecture & `sMILX` Stability Review, Fixes & Verification** (`phase3-qt6-gui-stability.md`) *(Pending)*
- **Phase 4: Biomedical Imaging & Scientific Visualization Assessment, Fixes & Verification** (`phase4-medical-imaging-workflows.md`) *(Pending)*
- **Phase 5: Final End-to-End Verification, Consolidated Report & Future Roadmap** (`phase5-remediation-roadmap.md`) *(Pending)*
- **Phase 6 (Deferred): QtPython (`PythonQt`) & Shiboken6 (`PySMILI`) Cross-Platform Bindings**
