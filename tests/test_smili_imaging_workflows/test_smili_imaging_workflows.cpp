/*=========================================================================
  Phase 4 Verification Test Suite:
  Biomedical Imaging & Scientific Visualization Workflows
=========================================================================*/
#include <iostream>
#include <cmath>
#include <vector>
#include <QApplication>
#include <QScopedPointer>
#include <QFileInfo>

#include "milxImage.h"
#include "milxModel.h"
#include "milxDeformableModel.h"
#include "milxFile.h"
#include "milxQtAliases.h"
#include "milxQtImage.h"
#include "milxQtModel.h"
#include "milxQtMain.h"

#include <vtkSphereSource.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataCollection.h>

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

    std::cout << "=== Test 1: Physical Coordinates, Origin/Direction in ResizeImage & ApplyOrientationToVTKImage ===" << std::endl;
    {
        using ImageType = floatImageType;
        ImageType::Pointer img = ImageType::New();
        ImageType::RegionType region;
        ImageType::IndexType start;
        start.Fill(0);
        ImageType::SizeType size;
        size[0] = 16; size[1] = 16; size[2] = 8;
        region.SetSize(size);
        region.SetIndex(start);
        img->SetRegions(region);

        ImageType::PointType origin;
        origin[0] = 12.5; origin[1] = -24.0; origin[2] = 8.0;
        img->SetOrigin(origin);

        ImageType::SpacingType spacing;
        spacing[0] = 1.5; spacing[1] = 1.5; spacing[2] = 2.0;
        img->SetSpacing(spacing);

        ImageType::DirectionType direction;
        direction.SetIdentity();
        direction(0, 0) = -1.0;
        direction(1, 1) = -1.0;
        img->SetDirection(direction);

        img->Allocate();
        img->FillBuffer(100.0f);

        ImageType::SizeType newSize;
        newSize[0] = 8; newSize[1] = 8; newSize[2] = 4;
        ImageType::SpacingType newSpacing;
        newSpacing[0] = 3.0; newSpacing[1] = 3.0; newSpacing[2] = 4.0;

        ImageType::Pointer resized1 = milx::Image<ImageType>::ResizeImage(img, newSize);
        TEST_ASSERT(std::abs(resized1->GetOrigin()[0] - 12.5) < 1e-4 &&
                    std::abs(resized1->GetOrigin()[1] - (-24.0)) < 1e-4 &&
                    std::abs(resized1->GetOrigin()[2] - 8.0) < 1e-4,
                    "ResizeImage(img, size) failed to preserve physical origin");
        TEST_ASSERT(std::abs(resized1->GetDirection()(0, 0) - (-1.0)) < 1e-4 &&
                    std::abs(resized1->GetDirection()(1, 1) - (-1.0)) < 1e-4,
                    "ResizeImage(img, size) failed to preserve direction cosines");

        ImageType::Pointer resized2 = milx::Image<ImageType>::ResizeImage(img, newSize, newSpacing);
        TEST_ASSERT(std::abs(resized2->GetOrigin()[0] - 12.5) < 1e-4 &&
                    std::abs(resized2->GetOrigin()[1] - (-24.0)) < 1e-4,
                    "ResizeImage(img, size, spacing) failed to preserve physical origin");
        TEST_ASSERT(std::abs(resized2->GetDirection()(0, 0) - (-1.0)) < 1e-4,
                    "ResizeImage(img, size, spacing) failed to preserve direction cosines");

        vtkSmartPointer<vtkImageData> vtkImg = vtkSmartPointer<vtkImageData>::New();
        vtkImg->DeepCopy(milx::Image<ImageType>::ConvertITKImageToVTKImage(img));
        vtkSmartPointer<vtkMatrix4x4> uninitMat; // null smart pointer
        vtkSmartPointer<vtkImageData> orientedVtk =
            milx::Image<ImageType>::ApplyOrientationToVTKImage(vtkImg, img, uninitMat, false, true);
        TEST_ASSERT(uninitMat != nullptr, "ApplyOrientationToVTKImage must allocate transformMatrix when null");
        TEST_ASSERT(orientedVtk != nullptr && orientedVtk->GetNumberOfPoints() > 0,
                    "ApplyOrientationToVTKImage returned empty output");
    }
    std::cout << "[PASS] Test 1" << std::endl;

    std::cout << "=== Test 2: 6-Component DTI / 4D Volume ConvertITKVectorImageToVTKImage ===" << std::endl;
    {
        vectorImageType::Pointer vecImg = vectorImageType::New();
        vectorImageType::RegionType region;
        vectorImageType::IndexType start;
        start.Fill(0);
        vectorImageType::SizeType size;
        size[0] = 6; size[1] = 6; size[2] = 4;
        region.SetSize(size);
        region.SetIndex(start);
        vecImg->SetRegions(region);
        vecImg->SetNumberOfComponentsPerPixel(6);
        vecImg->Allocate();

        itk::VariableLengthVector<float> pixel(6);
        for (unsigned int c = 0; c < 6; ++c) {
            pixel[c] = static_cast<float>((c + 1) * 10);
        }
        vecImg->FillBuffer(pixel);

        vtkSmartPointer<vtkImageData> vtkField =
            milx::Image<vectorImageType>::ConvertITKVectorImageToVTKImage(vecImg);
        TEST_ASSERT(vtkField != nullptr, "ConvertITKVectorImageToVTKImage returned null");
        TEST_ASSERT(vtkField->GetNumberOfScalarComponents() == 6,
                    "ConvertITKVectorImageToVTKImage must allocate 6 scalar components, not hardcoded 3");
        for (int c = 0; c < 6; ++c) {
            float val = vtkField->GetScalarComponentAsFloat(2, 2, 1, c);
            TEST_ASSERT(std::abs(val - static_cast<float>((c + 1) * 10)) < 1e-4f,
                        "Component value mismatch in 6-component VTK image");
        }
    }
    std::cout << "[PASS] Test 2" << std::endl;

    std::cout << "=== Test 3: VNL Vector & Matrix Import Buffer Deep-Copy Ownership ===" << std::endl;
    {
        floatImageType::Pointer importedVecImg;
        floatImageType::Pointer importedMatImg;
        {
            vnl_vector<float> tempVec(64, 7.5f);
            floatImageType::SizeType sz;
            sz[0] = 4; sz[1] = 4; sz[2] = 4;
            importedVecImg = milx::Image<floatImageType>::ImportVectorToImage(tempVec, sz);

            vnl_matrix<float> tempMat(8, 8, 3.25f);
            importedMatImg = milx::Image<floatImageType>::ImportMatrixToImage(tempMat);
        } // tempVec and tempMat are destroyed here

        floatImageType::IndexType idx = {{2, 2, 2}};
        TEST_ASSERT(std::abs(importedVecImg->GetPixel(idx) - 7.5f) < 1e-4f,
                    "ImportVectorToImage buffer invalid after vnl_vector destruction");
        floatImageType::IndexType idx2 = {{3, 3, 0}};
        TEST_ASSERT(std::abs(importedMatImg->GetPixel(idx2) - 3.25f) < 1e-4f,
                    "ImportMatrixToImage buffer invalid after vnl_matrix destruction");
    }
    std::cout << "[PASS] Test 3" << std::endl;

    std::cout << "=== Test 4: Surface Mesh & Collection Operations (ClearArrays, MaskScalars, Procrustes, DeformableModel) ===" << std::endl;
    {
        vtkSmartPointer<vtkSphereSource> sphere1 = vtkSmartPointer<vtkSphereSource>::New();
        sphere1->SetRadius(5.0);
        sphere1->SetCenter(10.0, 20.0, 30.0);
        sphere1->Update();

        vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
        poly->DeepCopy(sphere1->GetOutput());

        for (int a = 0; a < 4; ++a) {
            vtkSmartPointer<vtkFloatArray> arr = vtkSmartPointer<vtkFloatArray>::New();
            arr->SetName(("Array_" + std::to_string(a)).c_str());
            arr->SetNumberOfTuples(poly->GetNumberOfPoints());
            arr->FillComponent(0, static_cast<double>(a + 1));
            poly->GetPointData()->AddArray(arr);
        }
        poly->GetPointData()->SetActiveScalars("Array_0");
        TEST_ASSERT(poly->GetPointData()->GetNumberOfArrays() >= 4, "Expected at least 4 arrays before ClearArrays");

        milx::Model::ClearArrays(poly);
        TEST_ASSERT(poly->GetPointData()->GetNumberOfArrays() == 0,
                    "Model::ClearArrays failed to remove all point data arrays");

        // MaskScalars safety when model2 has no scalars
        vtkSmartPointer<vtkPolyData> polyWithScalars = vtkSmartPointer<vtkPolyData>::New();
        polyWithScalars->DeepCopy(sphere1->GetOutput());
        vtkSmartPointer<vtkFloatArray> sArr = vtkSmartPointer<vtkFloatArray>::New();
        sArr->SetName("Vals");
        sArr->SetNumberOfTuples(polyWithScalars->GetNumberOfPoints());
        sArr->FillComponent(0, 5.0);
        polyWithScalars->GetPointData()->SetScalars(sArr);

        // poly has 0 arrays/scalars now; MaskScalars must return safely without crashing
        milx::Model::MaskScalars(polyWithScalars, poly);
        TEST_ASSERT(std::abs(polyWithScalars->GetPointData()->GetScalars()->GetTuple1(0) - 5.0) < 1e-4,
                    "MaskScalars with un-scalared mask should leave model1 untouched");

        // ProcrustesAlignCollection & ScalarThresholdCollection
        vtkSmartPointer<vtkSphereSource> sphere2 = vtkSmartPointer<vtkSphereSource>::New();
        sphere2->SetRadius(5.0);
        sphere2->SetCenter(12.0, 18.0, 29.0);
        sphere2->Update();

        vtkSmartPointer<vtkPolyData> poly2 = vtkSmartPointer<vtkPolyData>::New();
        poly2->DeepCopy(sphere2->GetOutput());
        vtkSmartPointer<vtkFloatArray> sArr2 = vtkSmartPointer<vtkFloatArray>::New();
        sArr2->SetName("Vals");
        sArr2->SetNumberOfTuples(poly2->GetNumberOfPoints());
        sArr2->FillComponent(0, 5.0);
        poly2->GetPointData()->SetScalars(sArr2);

        vtkSmartPointer<vtkPolyDataCollection> col = vtkSmartPointer<vtkPolyDataCollection>::New();
        col->AddItem(polyWithScalars);
        col->AddItem(poly2);

        milx::Model modelProcessor;
        vtkSmartPointer<vtkPolyDataCollection> alignedCol = modelProcessor.ProcrustesAlignCollection(col, true);
        TEST_ASSERT(alignedCol != nullptr && alignedCol->GetNumberOfItems() == 2,
                    "ProcrustesAlignCollection failed to return 2 aligned meshes");

        // DeformableModel::ApplyOrientation with applyOrigin = true centers mesh on image origin
        using ImageType = floatImageType;
        ImageType::Pointer refImg = ImageType::New();
        ImageType::DirectionType dir;
        dir.SetIdentity();
        refImg->SetDirection(dir);
        ImageType::PointType refOrigin;
        refOrigin[0] = 10.0; refOrigin[1] = 20.0; refOrigin[2] = 30.0;
        refImg->SetOrigin(refOrigin);

        milx::DeformableModel defModel;
        defModel.SetInput(sphere1->GetOutput());
        defModel.ApplyOrientation<ImageType>(refImg, true, false);
        double outBounds[6];
        defModel.Result()->GetBounds(outBounds);
        double cx = 0.5 * (outBounds[0] + outBounds[1]);
        double cy = 0.5 * (outBounds[2] + outBounds[3]);
        double cz = 0.5 * (outBounds[4] + outBounds[5]);
        TEST_ASSERT(std::abs(cx - 10.0) < 0.5 && std::abs(cy - 20.0) < 0.5 && std::abs(cz - 30.0) < 0.5,
                    "DeformableModel::ApplyOrientation(applyOrigin=true) failed to place mesh at image origin");
    }
    std::cout << "[PASS] Test 4" << std::endl;

    std::cout << "=== Test 5: std::vector / NumPy Bridge in milxQtModel & milxQtImage ===" << std::endl;
    {
        QScopedPointer<milxQtModel> qtMdl(new milxQtModel);
        {
            std::vector<double> tempPts = {
                1.0, 2.0, 3.0,
                4.0, 5.0, 6.0,
                7.0, 8.0, 9.0,
                10.0, 11.0, 12.0
            };
            qtMdl->SetInput(tempPts);
        } // tempPts destroyed here

        TEST_ASSERT(qtMdl->GetOutput()->GetNumberOfPoints() == 4,
                    "milxQtModel::SetInput(vector<double>) of 12 elements must create 4 3D points, not 12");
        std::vector<double> roundTripPts = qtMdl->getData();
        TEST_ASSERT(roundTripPts.size() == 12, "milxQtModel::getData() should return 12 values");
        TEST_ASSERT(std::abs(roundTripPts[0] - 1.0) < 1e-6 && std::abs(roundTripPts[11] - 12.0) < 1e-6,
                    "milxQtModel::getData() point coordinates corrupted after tempPts destruction");

        // Also test getData() on a standard float-precision vtkSphereSource mesh
        vtkSmartPointer<vtkSphereSource> floatSphere = vtkSmartPointer<vtkSphereSource>::New();
        floatSphere->SetRadius(3.0);
        floatSphere->Update();
        qtMdl->SetInput(floatSphere->GetOutput());
        std::vector<double> floatMeshPts = qtMdl->getData();
        TEST_ASSERT(floatMeshPts.size() == static_cast<size_t>(floatSphere->GetOutput()->GetNumberOfPoints()) * 3,
                    "milxQtModel::getData() failed on VTK_FLOAT mesh");

        // Test milxQtImage::setData(const std::vector<float>&, ...)
        QScopedPointer<milxQtImage> qtImg(new milxQtImage);
        {
            std::vector<float> tempVol(10 * 8 * 6, 42.5f);
            tempVol[0] = 1.25f;
            tempVol[479] = 99.75f;
            qtImg->setData(tempVol, 10, 8, 6, 1.5);
        } // tempVol destroyed here

        std::vector<int> shp = qtImg->shape();
        TEST_ASSERT(shp.size() == 3 && shp[0] == 10 && shp[1] == 8 && shp[2] == 6,
                    "milxQtImage::shape() did not return {10, 8, 6}");
        std::vector<float> outData = qtImg->getData();
        TEST_ASSERT(outData.size() == 480 && std::abs(outData[0] - 1.25f) < 1e-4f && std::abs(outData[479] - 99.75f) < 1e-4f,
                    "milxQtImage::getData() failed after temporary input vector destruction");
        std::vector<float> outFloatData = qtImg->getFloatData();
        TEST_ASSERT(outFloatData.size() == 480,
                    "milxQtImage::getFloatData() returned empty vector after setData(vector<float>)");
    }
    std::cout << "[PASS] Test 5" << std::endl;

    std::cout << "=== Test 6: generateVoxelisedSurface Bounds Integrity & milxQtImage Interpolation Flags ===" << std::endl;
    {
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(5.0);
        sphere->Update();

        double bounds[6] = {-6.0, 6.0, -6.0, 6.0, -6.0, 6.0};
        QScopedPointer<milxQtImage> voxImg(new milxQtImage);
        voxImg->disableInterpolateDisplay();
        TEST_ASSERT(!voxImg->isInterpolated(), "disableInterpolateDisplay() before generateImage() must uncheck interpolateAct");
        voxImg->enableInterpolateDisplay();
        TEST_ASSERT(voxImg->isInterpolated(), "enableInterpolateDisplay() before generateImage() must check interpolateAct");

        double voxSpacing[3] = {1.0, 1.0, 1.0};
        voxImg->generateVoxelisedSurface(sphere->GetOutput(), bounds, voxSpacing);
        TEST_ASSERT(std::abs(bounds[0] - (-6.0)) < 1e-6 && std::abs(bounds[1] - 6.0) < 1e-6,
                    "generateVoxelisedSurface must not mutate caller's bounds array in-place");
    }
    std::cout << "[PASS] Test 6" << std::endl;

    std::cout << "=== Test 7: NIfTI Volume to Surface Mesh Patient-Space Registration (milxQtImage <-> milxQtModel) ===" << std::endl;
    {
        // 1. Synthetic NIfTI-like ITK volume matching AverageNonrigid_R_MR2_x2.nii.gz geometry:
        //    Origin = (21.09375, 8.25, 18.75), Spacing = (0.625, 1.5, 0.625), Size = (184, 74, 226)
        //    Physical patient space bounds:
        //      X in [21.09375, 135.46875], Y in [8.25, 117.75], Z in [18.75, 159.375]
        using ImageType = floatImageType;
        ImageType::Pointer niiImg = ImageType::New();
        ImageType::RegionType region;
        ImageType::IndexType start;
        start.Fill(0);
        ImageType::SizeType size;
        size[0] = 184; size[1] = 74; size[2] = 226;
        region.SetSize(size);
        region.SetIndex(start);
        niiImg->SetRegions(region);

        ImageType::PointType origin;
        origin[0] = 21.09375; origin[1] = 8.25; origin[2] = 18.75;
        niiImg->SetOrigin(origin);

        ImageType::SpacingType spacing;
        spacing[0] = 0.625; spacing[1] = 1.5; spacing[2] = 0.625;
        niiImg->SetSpacing(spacing);

        ImageType::DirectionType direction;
        direction.SetIdentity();
        niiImg->SetDirection(direction);
        niiImg->Allocate();
        niiImg->FillBuffer(50.0f);

        QScopedPointer<milxQtImage> imgWin(new milxQtImage);
        imgWin->SetInput(niiImg);
        imgWin->generateImage(true);

        // In milxQtImage display space (flipY = true across global Y=0), imageData Y bounds must be [-117.75, -8.25]
        double imgDisplayBounds[6];
        imgWin->GetOutput()->GetBounds(imgDisplayBounds);
        std::cout << "  milxQtImage display bounds: ["
                  << imgDisplayBounds[0] << ", " << imgDisplayBounds[1] << "] x ["
                  << imgDisplayBounds[2] << ", " << imgDisplayBounds[3] << "] x ["
                  << imgDisplayBounds[4] << ", " << imgDisplayBounds[5] << "]" << std::endl;
        TEST_ASSERT(std::abs(imgDisplayBounds[0] - 21.09375) < 1e-3 &&
                    std::abs(imgDisplayBounds[1] - 135.46875) < 1e-3,
                    "milxQtImage X bounds mismatch in display space");
        TEST_ASSERT(std::abs(imgDisplayBounds[2] - (-117.75)) < 1e-3 &&
                    std::abs(imgDisplayBounds[3] - (-8.25)) < 1e-3,
                    "milxQtImage Y bounds in display space must be [-117.75, -8.25] so SetUserMatrix(F) maps to [+8.25, +117.75]");
        TEST_ASSERT(std::abs(imgDisplayBounds[4] - 18.75) < 1e-3 &&
                    std::abs(imgDisplayBounds[5] - 159.375) < 1e-3,
                    "milxQtImage Z bounds mismatch in display space");

        // Now place the MRI actor into milxQtModel (via importFrom, which applies SetUserMatrix(transformMatrix^-1))
        QScopedPointer<milxQtModel> mdlWin(new milxQtModel);
        vtkSmartPointer<vtkSphereSource> kneeProxy = vtkSmartPointer<vtkSphereSource>::New();
        kneeProxy->SetCenter(80.4, 53.8, 92.0);
        kneeProxy->SetRadius(10.0);
        kneeProxy->Update();
        mdlWin->SetInput(kneeProxy->GetOutput());
        mdlWin->generateModel();
        mdlWin->importFrom(imgWin.data());

        // Compute world-space Y bounds of the MRI imageData after transformModel = transformMatrix^-1
        vtkSmartPointer<vtkMatrix4x4> transformModel = vtkSmartPointer<vtkMatrix4x4>::New();
        transformModel->DeepCopy(imgWin->getTransformMatrix());
        transformModel->Invert();
        double minPt[4] = {imgDisplayBounds[0], imgDisplayBounds[2], imgDisplayBounds[4], 1.0};
        double maxPt[4] = {imgDisplayBounds[1], imgDisplayBounds[3], imgDisplayBounds[5], 1.0};
        double worldMin[4], worldMax[4];
        transformModel->MultiplyPoint(minPt, worldMin);
        transformModel->MultiplyPoint(maxPt, worldMax);
        double worldYMin = std::min(worldMin[1], worldMax[1]);
        double worldYMax = std::max(worldMin[1], worldMax[1]);
        std::cout << "  Imported MRI actor world Y bounds in milxQtModel: ["
                  << worldYMin << ", " << worldYMax << "]" << std::endl;
        TEST_ASSERT(std::abs(worldYMin - 8.25) < 1e-3 && std::abs(worldYMax - 117.75) < 1e-3,
                    "Imported MRI actor world Y bounds in milxQtModel must match physical patient space [8.25, 117.75]");

        // 2. If the user's actual NIfTI & VTK dataset is available on disk, verify exact registration on it as well
        const QString niiPath = "E:/data/SyngoData_focused/AverageNonrigid_R_MR2_x2.nii.gz";
        const QString vtkPath = "E:/data/SyngoData_focused/AverageNonrigid_R_combined.vtk";
        if (QFileInfo::exists(niiPath) && QFileInfo::exists(vtkPath)) {
            std::cout << "  Verifying with actual dataset: " << niiPath.toStdString()
                      << " and " << vtkPath.toStdString() << std::endl;
            ImageType::Pointer realNii;
            TEST_ASSERT(milx::File::OpenImage<ImageType>(niiPath.toStdString(), realNii),
                        "Failed to open real NIfTI volume");
            vtkSmartPointer<vtkPolyData> realMesh;
            TEST_ASSERT(milx::File::OpenModel(vtkPath.toStdString(), realMesh),
                        "Failed to open real VTK mesh");

            QScopedPointer<milxQtImage> realImgWin(new milxQtImage);
            realImgWin->SetInput(realNii);
            realImgWin->generateImage(true);

            QScopedPointer<milxQtModel> realMdlWin(new milxQtModel);
            realMdlWin->SetInput(realMesh);
            realMdlWin->generateModel();

            double meshBounds[6];
            realMdlWin->GetOutput()->GetBounds(meshBounds);
            double realImgBounds[6];
            realImgWin->GetOutput()->GetBounds(realImgBounds);

            vtkSmartPointer<vtkMatrix4x4> realTransformModel = vtkSmartPointer<vtkMatrix4x4>::New();
            realTransformModel->DeepCopy(realImgWin->getTransformMatrix());
            realTransformModel->Invert();

            double p0[4] = {realImgBounds[0], realImgBounds[2], realImgBounds[4], 1.0};
            double p1[4] = {realImgBounds[1], realImgBounds[3], realImgBounds[5], 1.0};
            double wp0[4], wp1[4];
            realTransformModel->MultiplyPoint(p0, wp0);
            realTransformModel->MultiplyPoint(p1, wp1);

            double mriWorldBounds[6] = {
                std::min(wp0[0], wp1[0]), std::max(wp0[0], wp1[0]),
                std::min(wp0[1], wp1[1]), std::max(wp0[1], wp1[1]),
                std::min(wp0[2], wp1[2]), std::max(wp0[2], wp1[2])
            };
            std::cout << "  Real MRI actor world bounds in milxQtModel: ["
                      << mriWorldBounds[0] << ", " << mriWorldBounds[1] << "] x ["
                      << mriWorldBounds[2] << ", " << mriWorldBounds[3] << "] x ["
                      << mriWorldBounds[4] << ", " << mriWorldBounds[5] << "]" << std::endl;
            std::cout << "  Real VTK mesh bounds in milxQtModel:        ["
                      << meshBounds[0] << ", " << meshBounds[1] << "] x ["
                      << meshBounds[2] << ", " << meshBounds[3] << "] x ["
                      << meshBounds[4] << ", " << meshBounds[5] << "]" << std::endl;

            // The mesh X and Y bounds must lie strictly within the MRI volume's patient-space X and Y bounds
            TEST_ASSERT(meshBounds[0] >= mriWorldBounds[0] && meshBounds[1] <= mriWorldBounds[1],
                        "Real mesh X bounds fall outside MRI actor world X bounds");
            TEST_ASSERT(meshBounds[2] >= mriWorldBounds[2] && meshBounds[3] <= mriWorldBounds[3],
                        "Real mesh Y bounds fall outside MRI actor world Y bounds (misregistration!)");
        }
    }
    std::cout << "[PASS] Test 7" << std::endl;

    std::cout << "=== ALL PHASE 4 MEDICAL IMAGING & VISUALIZATION TESTS PASSED ===" << std::endl;
    return EXIT_SUCCESS;
}
