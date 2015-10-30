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
#ifndef __MILXDEFORMABLEMODEL_H
#define __MILXDEFORMABLEMODEL_H
//ITK
#include <itkImage.h>
#include <itkVTKPolyDataToMesh.h> //itkVTKGlue
#include <itkMeshToVTKPolyData.h> //itkVTKGlue
#include <itkLinearInterpolateImageFunction.h>
#ifdef ITK_USE_REVIEW //Review only members
  #include <itkConformalFlatteningMeshFilter.h>
#endif
//VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
//SMILI
#include <milxGlobal.h>
#include <milxMath.h>
#include <milxModel.h>
#include <milxImage.h>
#include <milxFile.h>

//enums

namespace milx
{
/**
  \class DeformableModel
  \brief Represents a deformable model (i.e. a model with cells and scalar values that has special members for interacting with images).

  It also handles interaction with ITK images and meshes.
*/
class SMILI_EXPORT DeformableModel : public Model
{
public:
  /*!
    \fn DeformableModel::DeformableModel()
    \brief Standard constructor
	*/
  DeformableModel();
  /*!
    \fn DeformableModel::DeformableModel(vtkSmartPointer<vtkPolyData> model)
    \brief Constructor that copies the input model.
	*/
  DeformableModel(vtkSmartPointer<vtkPolyData> model);
  /*!
    \fn DeformableModel::~DeformableModel()
    \brief Standard Destructor
	*/
  virtual ~DeformableModel() {}

  /**
  * \defgroup General Imaging General
  * \brief Members for operating on models with images, requires ITK
  */
  //@{
  /*!
    \fn DeformableModel::ApplyOrientation(itk::SmartPointer<TImage> image, const bool applyOrigin = true, const bool flipY = false)
    \brief Apply the orientation of an image to a model.

    This is useful sometimes when operations such as the marching cubes strips orientation from its result.
    flipY is needed if a conversion is made between ITK and VTK.
    */
  template<typename TImage>
  void ApplyOrientation(itk::SmartPointer<TImage> image, const bool applyOrigin = true, const bool flipY = false);
  /*!
    \fn DeformableModel::VoxeliseAsITKImage(const unsigned char insideValue, double *spacing, double *bounds = NULL, const size_t padVoxels = 1)
    \brief Convert a model to voxelised (imaging) data.

    Spacing is a length 3 array, bounds in a length 6 array (leave NULL to use output of GetBounds())
    padVoxels is the number of pixels to pad the voxelised surface by.
    This is useful sometimes when operations such as distance maps is needed.
    */
  template<typename TImage>
  itk::SmartPointer<TImage> VoxeliseAsITKImage(const unsigned char insideValue, double *spacing, double *bounds = NULL, const size_t padVoxels = 1);
  /*!
    \fn DeformableModel::ConvertVTKPolyDataToITKMesh()
    \brief Convert the current vtkPolyData to an itk::Mesh.
    */
  template<typename MeshType>
  itk::SmartPointer<MeshType> ConvertVTKPolyDataToITKMesh();
  /*!
    \fn DeformableModel::ConvertITKMeshToVTKPolyData(itk::SmartPointer<MeshType> mesh)
    \brief Convert an itk::Mesh to the current vtkPolyData.
    */
  template<typename MeshType>
  void ConvertITKMeshToVTKPolyData(itk::SmartPointer<MeshType> mesh);
#ifdef ITK_USE_REVIEW //Review only members
  /*!
    \fn DeformableModel::Flatten(const size_t flatMode = 1)
    \brief Given a genus zero (no holes) surface, 'flatten' it by conformally mapping the surface to a sphere/plane. If flatMode is 1, a sphere is produced.
    */
  void Flatten(const size_t flatMode = 1);
#endif
  //@}

  /**
  * \defgroup Interactions Imaging Interactions
  * \brief Members for operating on models interacting with images, requires ITK
  */
  //@{
  /*!
    \fn DeformableModel::SurfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<FloatImageType> img, const bool absoluteValues)
    \brief Returns the scalars of the mesh (float values per vertex) according to the values found within the image at corresponding pixel locations for each vertex.

    Takes into account of image orientation etc.
	*/
  template<typename TImage>
  static vtkSmartPointer<vtkFloatArray> SurfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<TImage> img, const bool absoluteValues);
  /**
   * \fn DeformableModel::MarkSurfaceInsideImage(vtkSmartPointer<vtkPolyData> &surface, itk::SmartPointer<TImage> img, vtkFloatArray *weights, const float outsideValue = 0.0)
   * \brief Mark existing surface with all parts outside the image as 'outsideValue' value. Inside will remain the same scalars if present, else will be marked as 1.0.
   * Uses include weights for the missing data etc.
   */
  template<typename TImage>
  static void MarkSurfaceInsideImage(vtkSmartPointer<vtkPolyData> &surface, itk::SmartPointer<TImage> img, vtkFloatArray *weights, const float outsideValue = 0.0);
  //@}

  //Collection operations
  /**
  * \defgroup DeformableModel_Collection Deformable Model Collection Operations
  * \brief Members for operating of a collection of deformable models (i.e. model batching)
  */
  //@{
  /**
   * \fn DeformableModel::ApplyOrientationCollection(vtkSmartPointer<vtkPolyDataCollection> collection, itk::SmartPointer<TImage> refImage, const bool applyOrigin = true, const bool flipY = false)
   * \brief Apply orientation/direction of a reference image to a collection of models.
   */
  template<typename TImage>
  void ApplyOrientationCollection(vtkSmartPointer<vtkPolyDataCollection> collection, itk::SmartPointer<TImage> refImage, const bool applyOrigin = true, const bool flipY = false);
  /**
   * \fn DeformableModel::VoxeliseCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType spacing, std::vector< typename itk::SmartPointer<TImage> > &images)
   * \brief Voxelise a collection of models with isotropic spacing given.
   * A vector of image pointers in set to 'images'
   */
  template<typename TImage>
  void VoxeliseCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType spacing, std::vector< typename itk::SmartPointer<TImage> > &images);
#ifdef ITK_USE_REVIEW //Review only members
  /**
   * \fn DeformableModel::FlattenCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t flatMode)
   * \brief Flatten a collection of models. See Flatten().
   */
  void FlattenCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t flatMode);
#endif
  //@}

protected:

private:

};

template<typename TImage>
itk::SmartPointer<TImage> DeformableModel::VoxeliseAsITKImage(const unsigned char insideValue, double *spacing, double *bounds, const size_t padVoxels)
{
  vtkSmartPointer<vtkImageData> img = Model::Voxelise(insideValue, spacing, bounds, padVoxels);

  return milx::Image<TImage>::ConvertVTKImageToITKImage(img);
}

template<typename MeshType>
itk::SmartPointer<MeshType> DeformableModel::ConvertVTKPolyDataToITKMesh()
{
  typedef itk::VTKPolyDataToMesh<MeshType>  MeshFilterType;
  typename MeshFilterType::Pointer filter = MeshFilterType::New();
  filter->SetInput(CurrentModel);
  try
  {
    std::cout << "Converting to ITK Mesh" << std::endl;
    filter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed ITK Conversion");
    PrintError(ex.GetDescription());
  }

  return filter->GetOutput();
}

template<typename MeshType>
void DeformableModel::ConvertITKMeshToVTKPolyData(itk::SmartPointer<MeshType> mesh)
{
  typedef itk::MeshToVTKPolyData<MeshType>  MeshFilterType;
  typename MeshFilterType::Pointer filter = MeshFilterType::New();
  filter->SetInput(mesh);
  try
  {
    std::cout << "Converting to VTK PolyData" << std::endl;
    filter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed VTK Conversion");
    PrintError(ex.GetDescription());
  }
  vtkSmartPointer<vtkPolyData> convMesh = filter->GetOutput();

  CurrentModel->DeepCopy(convMesh);
  std::cout << "Number of Points in VTK PolyData: " <<  CurrentModel->GetNumberOfPoints() << std::endl;
  std::cout << "Number of Cells in VTK PolyData: " <<  CurrentModel->GetNumberOfCells() << std::endl;
}

#ifdef ITK_USE_REVIEW //Review only members
void DeformableModel::Flatten(const size_t flatMode)
{
  //Triangulate
  Triangulate();

  //convert to ITK
  typedef itk::Mesh<vtkFloatingPointType, 3> MeshType;
  itk::SmartPointer<MeshType> mesh = ConvertVTKPolyDataToITKMesh<MeshType>();
  std::cout << "Number of Points in ITK Mesh: " <<  mesh->GetNumberOfPoints() << std::endl;

  //flatten
  typedef itk::ConformalFlatteningMeshFilter<MeshType, MeshType>  FlatFilterType;
  FlatFilterType::Pointer filter = FlatFilterType::New();
  filter->SetInput(mesh);
  //~ filter->SetPolarCellIdentifier(0);
  if(flatMode == 1)
      filter->MapToPlane();
  else
      filter->MapToSphere();
  try
  {
    std::cout << "Trying to flatten..." << std::endl;
    filter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Flattening");
    PrintError(ex.GetDescription());
  }
  itk::SmartPointer<MeshType> flatMesh = filter->GetOutput();

  //convert back to VTK
  ConvertITKMeshToVTKPolyData<MeshType>(flatMesh);
}
#endif

template<typename TImage>
void DeformableModel::ApplyOrientation(itk::SmartPointer<TImage> image, const bool applyOrigin, const bool flipY)
{
  ///Get orientation and origin
  typename TImage::DirectionType direction = image->GetDirection();
  typename TImage::PointType origin = image->GetOrigin();

  coordinate centroid = Math<double>::Centroid(CurrentModel->GetPoints());

  vtkSmartPointer<vtkMatrix4x4> flipMatrix = vtkSmartPointer<vtkMatrix4x4>::New(); //start with identity matrix
  flipMatrix->Identity();
  //Flip the image for VTK coordinate system
  if(flipY)
    flipMatrix->SetElement(1,1,-1); //flip

  std::cout << "Direction to be applied:\n" << direction << std::endl;
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New(); //start with identity matrix
    matrix->Identity();
    for (int i = 0; i < 3; i ++)
      for (int k = 0; k < 3; k ++)
        matrix->SetElement(i,k, direction(i,k));

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Identity();
    transform->PostMultiply();
    transform->Concatenate(flipMatrix); //flip
    transform->Translate(-centroid[0], -centroid[1], -centroid[2]); //remove centroid of mesh
    transform->Concatenate(matrix); //direction
  if(applyOrigin)
    transform->Translate(origin.GetDataPointer()); //add image origin displacement

  Model::SetTransform(transform);
}

template<typename TImage>
vtkSmartPointer<vtkFloatArray> DeformableModel::SurfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<TImage> img, const bool absoluteValues)
{
  const int numberOfPoints = surface->GetNumberOfPoints();

  typedef itk::Point<double, 3> InputImagePointType;
  typedef itk::ContinuousIndex<double, 3 > ContinuousIndexType;

  typedef itk::LinearInterpolateImageFunction<TImage, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(img);

  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
  scalars->SetName("Distance");
  scalars->SetNumberOfTuples(numberOfPoints);
  scalars->SetNumberOfComponents(1);
//  scalars->FillComponent(0, 0.0);

  InputImagePointType point;
  ContinuousIndexType index;
  for(int i = 0; i < numberOfPoints; i++)
  {
    double position[3];
    surface->GetPoint(i, position);

    point[0] = position[0];
    point[1] = position[1];
    point[2] = position[2];

    img->TransformPhysicalPointToContinuousIndex(point, index);

    // If inside, mark
    if(interpolator->IsInsideBuffer(index))
    {
      double valueFound = 0.0;
      if(absoluteValues)
        valueFound = fabs(interpolator->EvaluateAtContinuousIndex(index));
      else
        valueFound = interpolator->EvaluateAtContinuousIndex(index);
      scalars->SetTuple1(i, valueFound);
    }
  }

  return scalars;
}

template<typename TImage>
void DeformableModel::MarkSurfaceInsideImage(vtkSmartPointer<vtkPolyData> &surface, itk::SmartPointer<TImage> img, vtkFloatArray *weights, const float outsideValue)
{
  const int numberOfPoints = surface->GetNumberOfPoints();
  vtkSmartPointer<vtkFloatArray> surfaceWeights = vtkFloatArray::SafeDownCast(surface->GetPointData()->GetScalars());
  bool hasScalars = false;
  if(surfaceWeights)
    hasScalars = true;

  typedef itk::Point<double, 3> InputImagePointType;
  typedef itk::ContinuousIndex<double, 3> ContinuousIndexType;

  weights->SetName("Weights");
  weights->SetNumberOfTuples(numberOfPoints);
  weights->SetNumberOfComponents(1);
  weights->FillComponent(0, outsideValue);

  typedef itk::LinearInterpolateImageFunction<TImage, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(img);

//  std::cout << "Number of Points " << numberOfPoints << std::endl;
//  std::cout << "Image Spacing " << this->GetInput()->GetSpacing() << std::endl;
  InputImagePointType point;
  ContinuousIndexType index;
  for(int i = 0; i < numberOfPoints; i++)
  {
    double position[3];
    surface->GetPoint(i, position);

    point[0] = position[0];
    point[1] = position[1];
    point[2] = position[2];

    img->TransformPhysicalPointToContinuousIndex(point, index);

    // If inside, mark
    if(interpolator->IsInsideBuffer(index))
    {
      if(hasScalars)
        weights->SetTuple1(i, surfaceWeights->GetTuple1(i));
      else
        weights->SetTuple1(i, 1.0);
    }
  }

  surface->GetPointData()->SetScalars(weights);
}

//Collection ops
template<typename TImage>
void DeformableModel::ApplyOrientationCollection(vtkSmartPointer<vtkPolyDataCollection> collection, itk::SmartPointer<TImage> refImage, const bool applyOrigin, const bool flipY)
{
  const size_t n = collection->GetNumberOfItems();
  InternalInPlaceOperation = true;

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
    Model::SetInput(mesh);
    ApplyOrientation<TImage>(refImage, applyOrigin, flipY);
    mesh->DeepCopy(Model::Result());
  }

  InternalInPlaceOperation = false;
}

template<typename TImage>
void DeformableModel::VoxeliseCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType spacing, std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = collection->GetNumberOfItems();
  double isoSpacing[3];
  InternalInPlaceOperation = true;

  isoSpacing[0] = spacing;
  isoSpacing[1] = spacing;
  isoSpacing[2] = spacing;

  images.clear();
  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
    {
      vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
      Model::SetInput(mesh);
      itk::SmartPointer<TImage> img = VoxeliseAsITKImage<TImage>(1, isoSpacing); //!< scale the model inplace
      images.push_back(img);
    }

  InternalInPlaceOperation = false;
}

#ifdef ITK_USE_REVIEW //Review only members
void DeformableModel::FlattenCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t flatMode)
{
  const size_t n = collection->GetNumberOfItems();
  InternalInPlaceOperation = true;

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
    {
      vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
      Model::SetInput(mesh);
      Flatten(flatMode);
      mesh->DeepCopy(Model::Result());
    }

  InternalInPlaceOperation = false;
}
#endif

} //end namespace milx

#endif //__MILXDEFORMABLEMODEL_H
