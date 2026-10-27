/*=========================================================================
  Phase 3 Verification Test Suite:
  Qt 6 GUI Architecture & sMILX Stability Verification
=========================================================================*/
#include <iostream>
#include <cmath>
#include <QApplication>
#include <QMenu>
#include <QAction>
#include <QActionGroup>
#include <QMdiSubWindow>

#include "milxQtAliases.h"
#include "milxQtRenderWindow.h"
#include "milxQtImage.h"
#include "milxQtModel.h"
#include "milxQtPlot.h"
#include "milxQtUnifiedWindow.h"
#include "milxQtMain.h"

#include <vtkSphereSource.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>
#include <vtkProperty.h>

#define TEST_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "[FAIL] " << (msg) << " (line " << __LINE__ << ")" << std::endl; \
            return EXIT_FAILURE; \
        } \
    } while (0)

int main(int argc, char *argv[])
{
    milxQtVTKInit();
    QApplication app(argc, argv);

    std::cout << "=== Test 1: Uninitialized updateBar safety & Context Menu integrity ===" << std::endl;
    {
        using ImageType = floatImageType;
        ImageType::Pointer itkImg = ImageType::New();
        ImageType::RegionType region;
        ImageType::IndexType start;
        start.Fill(0);
        ImageType::SizeType size;
        size[0] = 16; size[1] = 16; size[2] = 16;
        region.SetSize(size);
        region.SetIndex(start);
        itkImg->SetRegions(region);
        itkImg->Allocate();
        itkImg->FillBuffer(10.0f);

        QScopedPointer<milxQtImage> qtImg(new milxQtImage);
        qtImg->setName("test_img.nii.gz");
        qtImg->SetInput(itkImg);
        qtImg->generateImage();

        // Call updateCoords without enableUpdates (updateBar == NULL)
        qtImg->updateCoords(qtImg->GetRenderWindow()->GetInteractor());

        // Verify createMenu populates both the inner "Imaging" submenu and outer actions
        QMenu outerImgMenu;
        qtImg->createMenu(&outerImgMenu);
        TEST_ASSERT(outerImgMenu.actions().size() > 3, "milxQtImage::createMenu should populate outer menu actions");
        QMenu *basicImgMenu = outerImgMenu.actions()[0]->menu();
        TEST_ASSERT(basicImgMenu != nullptr, "First item of milxQtImage::createMenu should be basicContextMenu");
        TEST_ASSERT(basicImgMenu->title() == "Imaging", "milxQtImage::basicContextMenu title should be 'Imaging'");

        // Test milxQtModel updateCoords & context menu
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(5.0);
        sphere->Update();

        QScopedPointer<milxQtModel> qtMdl(new milxQtModel);
        qtMdl->setName("test_sphere.vtk");
        qtMdl->SetInput(sphere->GetOutput());
        qtMdl->generateModel();

        qtMdl->updateCoords(qtMdl->GetRenderWindow()->GetInteractor());

        QMenu outerMdlMenu;
        qtMdl->createMenu(&outerMdlMenu);
        TEST_ASSERT(outerMdlMenu.actions().size() > 3, "milxQtModel::createMenu should populate outer menu actions");
        QMenu *basicMdlMenu = outerMdlMenu.actions()[0]->menu();
        TEST_ASSERT(basicMdlMenu != nullptr, "First item of milxQtModel::createMenu should be basicContextMenu");
        TEST_ASSERT(basicMdlMenu->title() == "Modelling", "milxQtModel::basicContextMenu title should be 'Modelling'");

        // Test milxQtPlot updateCoords without updateBar
        QScopedPointer<milxQtPlot> qtPlot(new milxQtPlot);
        qtPlot->updateCoords(qtPlot->GetRenderWindow()->GetInteractor());
    }
    std::cout << "[PASS] Test 1" << std::endl;

    std::cout << "=== Test 2: Programmatic non-interactive API calls (ok = true fix) ===" << std::endl;
    {
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(5.0);
        sphere->Update();

        QScopedPointer<milxQtModel> qtMdl(new milxQtModel);
        qtMdl->setName("opacity_sphere.vtk");
        qtMdl->SetInput(sphere->GetOutput());
        qtMdl->generateModel();

        qtMdl->changeOpacity(0.45f);
        double actualOpacity = qtMdl->GetActor()->GetProperty()->GetOpacity();
        TEST_ASSERT(std::abs(actualOpacity - 0.45) < 1e-3, "milxQtModel::changeOpacity(0.45f) did not apply");

        using ImageType = floatImageType;
        ImageType::Pointer itkImg = ImageType::New();
        ImageType::RegionType region;
        ImageType::IndexType start;
        start.Fill(0);
        ImageType::SizeType size;
        size[0] = 16; size[1] = 16; size[2] = 16;
        region.SetSize(size);
        region.SetIndex(start);
        itkImg->SetRegions(region);
        itkImg->Allocate();
        itkImg->FillBuffer(42.0f);

        QScopedPointer<milxQtImage> qtImg(new milxQtImage);
        qtImg->setName("subsample_img.nii.gz");
        qtImg->SetInput(itkImg);
        qtImg->generateImage();

        qtImg->subsample(2, 2, 2);
        auto subSize = qtImg->GetFloatImage()->GetLargestPossibleRegion().GetSize();
        TEST_ASSERT(subSize[0] == 8 && subSize[1] == 8 && subSize[2] == 8,
                    "milxQtImage::subsample(2, 2, 2) did not downsample 16^3 to 8^3");
    }
    std::cout << "[PASS] Test 2" << std::endl;

    std::cout << "=== Test 3: milxQtMain zero-window null-pointer safety ===" << std::endl;
    {
        milxQtMain mainWin;
        TEST_ASSERT(mainWin.activeName().isEmpty(), "activeName() with no windows should return empty QString");
        TEST_ASSERT(mainWin.activeNamePrefix().isEmpty(), "activeNamePrefix() with no windows should return empty QString");
        TEST_ASSERT(mainWin.activeWebView() == nullptr, "activeWebView() with no windows should return nullptr");
        TEST_ASSERT(!mainWin.isActiveWebView(), "isActiveWebView() with no windows should return false");

        // Must not crash when no subwindow is active
        mainWin.save("should_not_crash.nii.gz");
        mainWin.saveScreen("should_not_crash.png");
        mainWin.updateWindowsWithRefresh();
        mainWin.updateWindowsWithView(0);
        std::cout << "[PASS] Test 3" << std::endl;

        std::cout << "=== Test 4: milxQtMain MDI subwindow lifecycle, windowActionList & tiling ===" << std::endl;
        mainWin.resize(800, 600);
        mainWin.show();
        qApp->processEvents();

        using ImageType = floatImageType;
        ImageType::Pointer itkImg = ImageType::New();
        ImageType::RegionType region;
        ImageType::IndexType start;
        start.Fill(0);
        ImageType::SizeType size;
        size[0] = 16; size[1] = 16; size[2] = 8;
        region.SetSize(size);
        region.SetIndex(start);
        itkImg->SetRegions(region);
        itkImg->Allocate();
        itkImg->FillBuffer(25.0f);

        milxQtImage *imgWin = new milxQtImage;
        imgWin->setName("brain_scan.nii.gz");
        imgWin->SetInput(itkImg);
        imgWin->generateImage();
        mainWin.addImage(imgWin);

        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(10.0);
        sphere->Update();

        milxQtModel *mdlWin = new milxQtModel;
        mdlWin->setName("hippocampus.vtk");
        mdlWin->SetInput(sphere->GetOutput());
        mdlWin->generateModel();
        mainWin.addModel(mdlWin);

        qApp->processEvents();

        TEST_ASSERT(mainWin.getNumberOfWindows() == 2, "Expected 2 MDI subwindows in milxQtMain");
        QList<QMdiSubWindow*> subWins = mainWin.getListOfWindows();
        TEST_ASSERT(subWins.size() == 2, "getListOfWindows() should return 2 subwindows");

        // Activate first subwindow via setActiveWindow(QWidget*)
        mainWin.setActiveWindow(subWins[0]);
        qApp->processEvents();
        TEST_ASSERT(mainWin.getNumberOfWindows() == 2, "setActiveWindow must not create extra orphan QMdiSubWindows");
        TEST_ASSERT(!mainWin.activeName().isEmpty(), "activeName() should not be empty when subwindow is active");
        TEST_ASSERT(!mainWin.activeNamePrefix().isEmpty(), "activeNamePrefix() should not be empty when subwindow is active");
        TEST_ASSERT(!mainWin.isActiveWebView(), "isActiveWebView() on VTK render window must return false without crashing");

        // Test windowActionList with applyMapper = true (uses QAction::triggered lambda)
        QMenu winMenu;
        QActionGroup *grp = mainWin.windowActionList(&winMenu, true, true);
        TEST_ASSERT(grp != nullptr, "windowActionList returned null group");
        TEST_ASSERT(winMenu.actions().size() == 2, "windowActionList should create 2 actions for 2 subwindows");

        // Trigger second window action and verify it activates subWins[1] without adding subwindows
        winMenu.actions()[1]->trigger();
        qApp->processEvents();
        TEST_ASSERT(mainWin.getNumberOfWindows() == 2, "Triggering window action must not increase subwindow count");

        // Test vertical and horizontal tiling
        mainWin.tileTabVertically();
        qApp->processEvents();
        mainWin.tileTabHorizontally();
        qApp->processEvents();

        mainWin.closeTabAllWindows();
        qApp->processEvents();
    }
    std::cout << "[PASS] Test 4" << std::endl;

    std::cout << "=== Test 5: milxQtUnifiedWindow add & removeFromWindow lifecycle ===" << std::endl;
    {
        using ImageType = floatImageType;
        ImageType::Pointer itkImg = ImageType::New();
        ImageType::RegionType region;
        ImageType::IndexType start;
        start.Fill(0);
        ImageType::SizeType size;
        size[0] = 8; size[1] = 8; size[2] = 8;
        region.SetSize(size);
        region.SetIndex(start);
        itkImg->SetRegions(region);
        itkImg->Allocate();
        itkImg->FillBuffer(5.0f);

        QScopedPointer<milxQtImage> qtImg(new milxQtImage);
        qtImg->setName("uni_img.nii.gz");
        qtImg->SetInput(itkImg);
        qtImg->generateImage();

        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(4.0);
        sphere->Update();

        QScopedPointer<milxQtModel> qtMdl(new milxQtModel);
        qtMdl->setName("uni_mdl.vtk");
        qtMdl->SetInput(sphere->GetOutput());
        qtMdl->generateModel();

        QScopedPointer<milxQtUnifiedWindow> uniWin(new milxQtUnifiedWindow);
        uniWin->addToWindow(qtImg.data());
        uniWin->addToWindow(qtMdl.data());

        // Remove both and delete qtImg/qtMdl before calling refresh()
        uniWin->removeFromWindow(qtImg.data());
        uniWin->removeFromWindow(qtMdl.data());
        qtImg.reset();
        qtMdl.reset();

        // Must not access freed qtImg/qtMdl pointers
        uniWin->refresh();
    }
    std::cout << "[PASS] Test 5" << std::endl;

    std::cout << "=== ALL PHASE 3 QT 6 GUI & sMILX TESTS PASSED ===" << std::endl;
    return EXIT_SUCCESS;
}
