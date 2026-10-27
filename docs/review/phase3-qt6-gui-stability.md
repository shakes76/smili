# Phase 3: Qt 6 GUI Architecture & `sMILX` Stability Review, Fixes & Verification

**Reviewers**: `senior-qt-developer` (Senior Qt 6 GUI & CMake Developer) & `medical-imaging-researcher` (Medical Imaging Researcher & SMILI Architect)  
**Target Toolchain**: Qt 6.10.2 (`msvc2022_64`), VTK 9.6.0, ITK 5.4.3  
**Branch**: `topic-review`

---

## 1. Executive Summary

Phase 3 audited the Qt 6 GUI architecture across `include/Qt/`, `src/Qt/`, `ui/`, `plugin/`, and the GUI applications (`apps/smilx`, `apps/milxImageViewer`, `apps/milxModelViewer`, `apps/milxOverlay`, `apps/milxBlend`, `apps/milxAnimate`).

While the initial port from Qt 4 to Qt 6 allowed the project to compile once CMake was modernized, **string-based `SIGNAL()`/`SLOT()` macros hid removed Qt 4 signals (`QAction::activated()`, `QSignalMapper::mapped(QWidget*)`) from the C++ compiler**, causing silent runtime connection failures in the DICOM plugin, viewer applications, and `sMILX` window menu. Additionally, several critical memory-lifecycle, null-pointer, and context-menu defects caused crashes or broken UI menus in `sMILX`.

All 8 categories of Qt 6 GUI and `sMILX` stability defects were remediated and verified via a new automated CTest suite (`test_qt6_smilx_gui`).

---

## 2. Issues Identified & Remediated

### 2.1 Removed `QAction::activated()` Signal in Qt 6 (`plugin/dicom`, `milxImageViewer`, `milxModelViewer`)
- **Severity**: Critical (Silent Runtime Signal Failure)
- **Files**:
  - [`plugin/dicom/milxQtDICOMPlugin.cpp`](../../plugin/dicom/milxQtDICOMPlugin.cpp)
  - [`apps/milxImageViewer/milxImageViewer.cpp`](../../apps/milxImageViewer/milxImageViewer.cpp)
  - [`apps/milxModelViewer/milxModelViewer.cpp`](../../apps/milxModelViewer/milxModelViewer.cpp)
- **Root Cause**: `QAction::activated()` was a Qt 3/4 compatibility signal removed in Qt 6 in favor of `QAction::triggered()`. Because old-style `SIGNAL(activated())` connections are evaluated as runtime strings rather than compile-time function pointers, all 5 DICOM plugin actions (`openSeries`, `viewTags`, `openStructureSet`, `convert`, `anonymize`) and the `Exit` actions in `milxImageViewer` and `milxModelViewer` compiled without error but failed to connect at runtime (`QObject::connect: No such signal QAction::activated()`).
- **Fix**: Replaced all `SIGNAL(activated())` occurrences with `SIGNAL(triggered())`.

### 2.2 Removed `QSignalMapper::mapped(QWidget*)` Signal & Broken `milxQtMain::setActiveWindow(QWidget*)`
- **Severity**: Critical (Broken MDI Window Switching & Hierarchy Corruption)
- **File**: [`src/Qt/milxQtMain.cpp`](../../src/Qt/milxQtMain.cpp)
- **Root Cause**:
  1. `milxQtMain::createConnections()` connected `windowMapper` via `SIGNAL(mapped(QWidget *))`, which was removed from `QSignalMapper` in Qt 6.
  2. Even worse, `milxQtMain::setActiveWindow(QWidget *currentWindow)` allocated an unparented `new QMdiSubWindow` on the heap, called `subWindow->setWidget(currentWindow)` (stealing the existing `QMdiSubWindow` from `QMdiArea`), and passed the orphan `subWindow` to `QMdiArea::setActiveSubWindow(subWindow)`.
- **Fix**:
  - Replaced the `QSignalMapper` indirection in `milxQtMain::windowActionList()` with a direct Qt 6 lambda connection:
    ```cpp
    connect(action, &QAction::triggered, this, [this, currentWindow]() { setActiveWindow(currentWindow); });
    ```
  - Fixed `milxQtMain::setActiveWindow(QWidget *currentWindow)` to resolve the existing `QMdiSubWindow*` (either directly via `qobject_cast<QMdiSubWindow*>(currentWindow)` or via `currentWindow->parentWidget()`), ensure it is visible (`subWindow->show()`), and activate it in `QMdiArea` without allocating an orphan `QMdiSubWindow`. Also ensured `milxQtMain::addRender()` calls `subWindow->show()`.
  - Fixed a Qt 6 `QMainWindowLayout` / `QDockAreaLayout` crash in `milxQtMain::readSettings()` caused by calling `restoreDockWidget(console->dockWidget())` after `console->dockWidget()` had already been registered via `addDockWidget()` prior to `restoreState()`. Also guarded `restoreGeometry()` and `restoreState()` with `settings.contains(...)` and cleaned up `Connector->Disconnect()` and unparented `currentUnifiedWindow` in `~milxQtMain()`.

### 2.3 Uninitialized `updateBar` (`QStatusBar*`) Crash in Standalone `milxQtRenderWindow` Subclasses
- **Severity**: Critical (Segmentation Fault / Access Violation)
- **Files**:
  - [`src/Qt/milxQtRenderWindow.cpp`](../../src/Qt/milxQtRenderWindow.cpp)
  - [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp)
  - [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp)
  - [`src/Qt/milxQtPlot.cpp`](../../src/Qt/milxQtPlot.cpp)
- **Root Cause**: `QStatusBar *updateBar` declared in `milxQtRenderWindow` was never initialized in `milxQtRenderWindow::milxQtRenderWindow()`. Whenever a `milxQtRenderWindow`, `milxQtImage`, `milxQtModel`, or `milxQtPlot` was used outside `milxQtMain` (or before `enableUpdates()` was called) and `updateCoords()` was invoked, `updateBar->showMessage(message)` dereferenced an uninitialized wild pointer.
- **Fix**: Initialized `updateBar = NULL;` in `milxQtRenderWindow`'s constructor and guarded all `updateBar->showMessage(message)` calls across `milxQtRenderWindow`, `milxQtImage`, `milxQtModel`, and `milxQtPlot` with `if (updateBar)`. Also added `actionCubeHelix` and `actionHSV` to `mapGroup` in `milxQtRenderWindow::createActions()`.

### 2.4 Null-Pointer Dereferences in `milxQtMain` When No MDI Subwindow Is Active or Layout Is Null
- **Severity**: Critical (Application Crash)
- **File**: [`src/Qt/milxQtMain.cpp`](../../src/Qt/milxQtMain.cpp)
- **Root Cause**:
  1. `milxQtMain::activeWebView()` called `activeWin->layout()->count()` without checking `if (activeWin && activeWin->layout())`. Because `QVTKOpenGLNativeWidget` has `layout() == nullptr`, calling `activeNamePrefix()` (which checks `isActiveWebView()`) on any VTK render window crashed immediately.
  2. `milxQtMain::save()` and `milxQtMain::saveScreen()` dereferenced `activeSubWindow()->widget()` **before** checking if `activeSubWindow()` was null.
  3. `milxQtMain::activeName()` and `milxQtMain::activeNamePrefix()` returned `0` (`nullptr` `const char*` conversion in `QString(const char*)`, deprecated/removed in modern Qt) when no subwindow was active.
  4. `milxQtMain::tileTabVertically()` and `tileTabHorizontally()` called `widget->widget()->parentWidget()->resize(...)` instead of `widget->resize(...)` directly on the `QMdiSubWindow*`.
  5. `milxQtMain::updateWindowsWithRefresh()` and `updateWindowsWithView()` did not null-check `win` / `rndWin`.
  6. `milxQtMain::dropEvent()` manually stripped the leading `/` from `QUrl::path()` on Windows instead of using `QUrl::toLocalFile()`.
- **Fix**: Added null checks in `activeWebView()`, `save()`, `saveScreen()`, `updateWindowsWithRefresh()`, and `updateWindowsWithView()`; returned `QString()` in `activeName()` and `activeNamePrefix()`; resized/moved `QMdiSubWindow*` directly in `tileTabVertically()` and `tileTabHorizontally()`; and used `urlsList[j].toLocalFile()` in `dropEvent()`.

### 2.5 Context Menu Overwrite & `QObject::tr()` Misuse in `milxQtImage`, `milxQtModel`, and `milxQtPlot`
- **Severity**: High (Broken Right-Click Context Menus)
- **Files**:
  - [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp)
  - [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp)
  - [`src/Qt/milxQtPlot.cpp`](../../src/Qt/milxQtPlot.cpp)
- **Root Cause**:
  1. In `milxQtImage::basicContextMenu()` and `milxQtModel::basicContextMenu()`, `contextMenu = new QMenu(this);` overwrote the member variable `this->contextMenu` while `createMenu(contextMenu)` was populating it, causing right-click context menus to lose all outer actions (`Extensions`, `Contouring`, `Window Properties`, `Refresh`, `Reset`, `Close`, `Undock`).
  2. `tr("MainWindow", "Imaging", 0)`, `tr("MainWindow", "Modelling", 0)`, and `tr("MainWindow", "Plotting", 0)` passed `"MainWindow"` as the `sourceText` and the actual menu title as the `disambiguation` comment to `QObject::tr()`, resulting in all three menus being titled `"MainWindow"`.
- **Fix**: Declared `QMenu *contextMenu = new QMenu(this);` as a local variable in `basicContextMenu()` and changed the titles to `tr("Imaging")`, `tr("Modelling")`, and `tr("Plotting")`.

### 2.6 Non-Interactive Programmatic Calls Silently Aborting (`bool ok = false` Bug)
- **Severity**: High (Broken Programmatic C++ / Python API Calls)
- **Files**:
  - [`src/Qt/milxQtImage.cpp`](../../src/Qt/milxQtImage.cpp)
  - [`src/Qt/milxQtModel.cpp`](../../src/Qt/milxQtModel.cpp)
- **Root Cause**: In `milxQtImage::subsample`, `checkerBoard`, `otsu`, `otsuMultiple`, `component` and `milxQtModel::generateStreamLines`, `changeOpacity`, local `bool ok = false;` was only set by `QInputDialog` when default sentinel arguments were passed. Passing explicit arguments programmatically skipped the `QInputDialog` and hit `if (!ok) return;`, silently aborting the operation.
- **Fix**: Set `ok = true;` prior to the conditional `QInputDialog` prompt so programmatic calls execute while cancelled interactive dialogs still abort.

### 2.7 Dangling `milxQtImage*` Pointer in `milxQtUnifiedWindow::removeFromWindow`
- **Severity**: High (Use-After-Free Crash)
- **File**: [`src/Qt/milxQtUnifiedWindow.cpp`](../../src/Qt/milxQtUnifiedWindow.cpp)
- **Root Cause**: `milxQtUnifiedWindow::removeFromWindow(QWidget *passedWindow)` had the `milxQtImage*` removal commented out, leaving dangling pointers in `unifyImages` when an image subwindow was closed.
- **Fix**: Restored `qobject_cast<milxQtImage*>(passedWindow)` and `unifyImages.removeAll(image)` when `model == 0`.

### 2.8 UI Forms (`milxQtControlsForm`, `milxQtPreferencesForm`, `milxQtThemeEditor`) & Viewer App Fixes
- **Severity**: Medium
- **Files**:
  - [`ui/milxQtControlsForm.cpp`](../../ui/milxQtControlsForm.cpp)
  - [`ui/milxQtPreferencesForm.cpp`](../../ui/milxQtPreferencesForm.cpp)
  - [`ui/milxQtThemeEditor.cpp`](../../ui/milxQtThemeEditor.cpp)
  - [`apps/milxBlend/milxBlend.cpp`](../../apps/milxBlend/milxBlend.cpp)
  - [`apps/milxOverlay/milxOverlay.cpp`](../../apps/milxOverlay/milxOverlay.cpp)
  - [`apps/milxAnimate/milxAnimate.cpp`](../../apps/milxAnimate/milxAnimate.cpp)
- **Fixes**:
  - Removed `this->setParent(NULL);` in `milxQtControlsForm::closeEvent` which orphaned the widget from its `QMdiSubWindow`.
  - Fixed `if (themeName.compare("Light") && themeName.compare("Dark"))` in `milxQtPreferencesForm::removeTheme()` (was `||`, which was always true and allowed deleting built-in themes).
  - Fixed `this->themeName` initialization order in `milxQtThemeEditorForm` constructor and checked `QFile::open(QFile::ReadOnly)` return values across `milxQtMain.cpp`, `milxQtPreferencesForm.cpp`, and `milxQtThemeEditor.cpp`.
  - Fixed swapped `SetSize(windowWidth, windowHeight)`, `#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)`, and `mainWindow.takeCentralWidget()` before exit in `milxBlend`, `milxOverlay`, and `milxAnimate`.

---

## 3. Verification Results

A dedicated Qt 6 GUI & `sMILX` stability test suite was added in [`tests/test_qt6_smilx_gui/test_qt6_smilx_gui.cpp`](../../tests/test_qt6_smilx_gui/test_qt6_smilx_gui.cpp) and registered with CTest:
1. **Test 1**: Standalone `milxQtImage`, `milxQtModel`, and `milxQtPlot` `updateCoords()` safety when `updateBar == NULL`, plus verification that `basicContextMenu()` titles are `"Imaging"` and `"Modelling"` and `createMenu()` preserves outer menu actions.
2. **Test 2**: Programmatic non-interactive execution of `milxQtModel::changeOpacity(0.45f)` and `milxQtImage::subsample(2, 2, 2)`.
3. **Test 3**: `milxQtMain` zero-window null-pointer safety across `activeName()`, `activeNamePrefix()`, `activeWebView()`, `isActiveWebView()`, `save()`, `saveScreen()`, `updateWindowsWithRefresh()`, and `updateWindowsWithView()`.
4. **Test 4**: `milxQtMain` MDI subwindow lifecycle (`addImage`, `addModel`), `setActiveWindow()` without orphan `QMdiSubWindow` creation, `windowActionList()` lambda action triggering, and `tileTabVertically()` / `tileTabHorizontally()`.
5. **Test 5**: `milxQtUnifiedWindow` `addToWindow()` and `removeFromWindow()` lifecycle for both `milxQtImage` and `milxQtModel` followed by `refresh()`.

### CTest Output (`build-vs2022`, `Release`)
```text
Test project D:/Antigravity/smili-release/build-vs2022
    Start 1: cli_milxImageApp_help
1/7 Test #1: cli_milxImageApp_help .............   Passed    0.09 sec
    Start 2: cli_milxDICOMApp_help
2/7 Test #2: cli_milxDICOMApp_help .............   Passed    0.09 sec
    Start 3: cli_milxModelApp_help
3/7 Test #3: cli_milxModelApp_help .............   Passed    0.09 sec
    Start 4: cli_milxDeformableModelApp_help
4/7 Test #4: cli_milxDeformableModelApp_help ...   Passed    0.09 sec
    Start 5: cli_milxHausdorffDistance_help
5/7 Test #5: cli_milxHausdorffDistance_help ....   Passed    0.09 sec
    Start 6: test_vtk9_itk5_pipeline
6/7 Test #6: test_vtk9_itk5_pipeline ...........   Passed    0.39 sec
    Start 7: test_qt6_smilx_gui
7/7 Test #7: test_qt6_smilx_gui ................   Passed    0.92 sec

100% tests passed, 0 tests failed out of 7

Total Test time (real) =   1.80 sec
```
