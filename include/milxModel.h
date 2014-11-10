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
#ifndef __MILXMODEL_H
#define __MILXMODEL_H

//VTK
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkImageData.h> //voxelise
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyDataCollection.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkTransform.h>
#include <vtkTransformCollection.h>
#include <vtkFloatArray.h>
#include <vtkMutableUndirectedGraph.h>

#include <milxGlobal.h>

//enums
enum GlyphType {Circle = 0, Triangle, Cross, Square, Diamond};

namespace milx
{
/**
  \class Model
  \brief Represents a model (i.e. a model with cells and scalar values) and their common operations. Also allows batch operations on collection of models.

  The class members are divided into two groups, atmoic and collection members. Atomic operations operate on one model at a time, whereas collection members use the atomic
  operations in a loop to batch that operation over multiple models.

  Unlike the milx::Image class, this class tracks the current state of the data using pointers which can be retrieved using the Result() (or GetOutput()) and PreviousResult() members.
  This class is used extensively throughout the milxQtModel class and in milxSurfaceDiagnosticApp application.
  See the usage examples given below:

  Processing Pipeline
  \code
  milx::Model model;
	  model.SetInput(surface);
	  model.Triangulate();
	  model.LaplacianSmoothing(250); //Kill some of the steps caused by MCs
	  model.WindowedSincSmoothing(20); //Preserve features
	  model.Clean();
	  model.Update();
  \endcode

  ICP Registration
  \code
  ...
  vtkSmartPointer<vtkPolyData> mesh = m_RawSurfaces->GetNextItem();

  milx::Model initialisation(m_AtlasSurface);
  initialisation.MatchInformation(mesh, true); //Get correct approx scale
  initialisation.IterativeClosestPointsRegistration(mesh, false); //Doesn't handle different FoV's
  \endcode

  Landmark based registration
  \code
  milx::Model model(m_AtlasSurface);
  model.LandmarkBasedRegistration(mesh, true); ///requires correspondence
  mesh->DeepCopy(model.GetOutput());
  \endcode

  Voxelise model
  \code
  milx::Model model1(surface1);
  vtkSmartPointer<vtkImageData> voxelisedModel1 = model1.Voxelise(255, labelSpacing.GetDataPointer(), bigBounds);
  \endcode

  Multi-use example
  \code
	//Read labels and generate iso surface
	//Clip mesh to ensure same FoV
	vtkSmartPointer<vtkFloatArray> weights = vtkSmartPointer<vtkFloatArray>::New();
	milx::Model model(surface);
	model.MarkSurfaceInsideImage<FloatImageType>(surface, distanceMap, weights); //sets weights as scalars in here
	model.Clip(1.0, 1.0);

	//Save distances as scalars on mesh
	const bool absValues = true;
	vtkSmartPointer<vtkFloatArray> scalars = model.SurfaceScalarsFromImage<FloatImageType>(model.GetOutput(), distanceMap, absValues);
	scalars->SetName("Surface Distance Errors");
	model.SetScalars(scalars);

	//Restore correspondence to copy scalars over to full (unclipped) mesh, since the mesh is clipped to ensure same FoV
	cout << "Regions of Surface Outside the Image are marked with -1." << endl;
	weights->FillComponent(0, -1.0);
	for(int j = 0; j < model.GetOutput()->GetNumberOfPoints(); j ++)
	{
	  vtkIdType index = surface->FindPoint(model.GetOutput()->GetPoint(j));
	  weights->SetTuple1(index, scalars->GetTuple1(j));
	}
  \endcode

*/
class SMILI_EXPORT Model
{
public:
  /*!
    \fn Model::Model()
    \brief Standard constructor
	*/
  Model();
  /*!
    \fn Model::Model(vtkSmartPointer<vtkPolyData> model)
    \brief Constructor that copies the input model.
	*/
  Model(vtkSmartPointer<vtkPolyData> model);
  /*!
    \fn Model::~Model()
    \brief Standard Destructor
	*/
  virtual ~Model() {}

  /**
  * \defgroup IO Class IO Members
  * \brief Members for model input/output
  */
  //@{
  /*!
    \fn Model::AddInput(vtkSmartPointer<vtkPolyData> model)
    \brief Assigns and coninually appends meshes provided to the class.
	*/
  void AddInput(vtkSmartPointer<vtkPolyData> model);
  /*!
    \fn Model::AddArray(vtkSmartPointer<vtkDataArray> array)
    \brief Appends an array to the model.
	*/
  inline void AddArray(vtkSmartPointer<vtkDataArray> array)
  { CurrentModel->GetPointData()->AddArray(array);  }
  /*!
    \fn Model::SetInput(vtkSmartPointer<vtkPolyData> model)
    \brief Assigns the input model to class
	*/
  void SetInput(vtkSmartPointer<vtkPolyData> model);
  /*!
    \fn Model::SetInputData(vtkSmartPointer<vtkPolyData> model)
    \brief Assigns the input model to class. Same as SetInput and meant for VTK version > 5.
	*/
  inline void SetInputData(vtkSmartPointer<vtkPolyData> model)
  { SetInput(model); }
  /*!
    \fn Model::SetTransform(vtkSmartPointer<vtkTransform> newTransform)
    \brief Transforms the model by given transform
	*/
  void SetTransform(vtkSmartPointer<vtkAbstractTransform> newTransform);
  /*!
    \fn Model::SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph)
    \brief Converts a graph to polydata/mesh ready for display.
	*/
  void SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph);
  /*!
    \fn Model::SetPoints(vtkSmartPointer<vtkPoints> modelPoints)
    \brief Assign points to the model
	*/
  void SetPoints(vtkSmartPointer<vtkPoints> modelPoints);
  /*!
    \fn Model::SetPolys(vtkSmartPointer<vtkCellArray> modelPolys)
    \brief Assign polygons to the model via the array given.
	*/
  void SetPolys(vtkSmartPointer<vtkCellArray> modelPolys);
  /*!
    \fn Model::SetVectors(vtkSmartPointer<vtkDataArray> modelVectors)
    \brief Sets the vector field of the mesh (vectors per vertex of the surface)
	*/
  void SetVectors(vtkSmartPointer<vtkDataArray> modelVectors);
  /*!
    \fn Model::SetScalars(vtkSmartPointer<vtkDataArray> modelScalars)
    \brief Sets the scalar field of the mesh (values per vertex of the surface)
	*/
  void SetScalars(vtkSmartPointer<vtkDataArray> modelScalars);
  /*!
    \fn Model::SetActiveScalars(vtkSmartPointer<vtkDataArray> modelScalars)
    \brief Sets the array of name as the scalar field of the mesh (values per vertex of the surface)
	*/
  inline void SetActiveScalars(std::string nameOfArray)
  { CurrentModel->GetPointData()->SetActiveScalars(nameOfArray.c_str());  }

  /*!
    \fn Model::Result()
    \brief Returns the current model, i.e. the result of the latest operation.

  Could be NULL, the user must check.
	*/
  inline vtkSmartPointer<vtkPolyData>& Result()
  { return CurrentModel;  }
  /*!
    \fn Model::PreviousResult()
    \brief Returns the previous model, i.e. the result of the penultimate operation.

  Could be NULL, the user must check.
	*/
  inline vtkSmartPointer<vtkPolyData>& PreviousResult()
  { return PreviousModel;  }
  /*!
    \fn Model::GetInput()
    \brief Returns the input model, i.e. the very most initial model.

  Could be NULL, the user must check.
	*/
  inline vtkSmartPointer<vtkPolyData>& GetInput()
  { return InputModel;  }
  /*!
    \fn Model::GetOutput()
    \brief Returns the current model, i.e. the result of the latest operation ITK/VTK style.
	*/
  inline vtkSmartPointer<vtkPolyData>& GetOutput()
  { return Result();  }
#if VTK_MAJOR_VERSION <=5
  /*!
    \fn Model::GetProducerPort()
    \brief Returns the current model, i.e. the result of the latest operation ITK/VTK style.
	*/
  inline vtkAlgorithmOutput* GetProducerPort()
  { return Result()->GetProducerPort();  }
#endif
  /*!
      \fn Model::GetNumberOfPoints()
      \brief Returns the total number of points of the model currently held.
  */
  inline vtkIdType GetNumberOfPoints()
  { return Result()->GetNumberOfPoints(); }
  /*!
      \fn Model::GetNumberOfArrays()
      \brief Returns the total number of point-based arrays in the model.
  */
  inline vtkIdType GetNumberOfArrays()
  { return Result()->GetPointData()->GetNumberOfArrays(); }
  /*!
    \fn Model::GetTransform()
    \brief Returns the current net transform applied to the model, i.e. concatenations of all (rigid/linear) transforms.
	*/
  inline vtkSmartPointer<vtkTransform>& GetTransform()
  { return CurrentTransform;  }
  /*!
    \fn Model::GetPoints()
    \brief Returns the points of the model.
	*/
  inline vtkSmartPointer<vtkPoints> GetPoints()
  { return Result()->GetPoints();  }
  /*!
    \fn Model::GetScalars()
    \brief Returns the active scalars of the model.
	*/
  inline vtkSmartPointer<vtkDataArray> GetScalars()
  { return Result()->GetPointData()->GetScalars();  }
  /*!
    \fn Model::GetScalarRange(double *range)
    \brief Returns the scalar range of the model.
    */
  inline void GetScalarRange(double *range)
  { if(GetScalars()) Result()->GetScalarRange(range);  }
  /*!
    \fn Model::GetVectors()
    \brief Returns the active vectors of the model.
	*/
  inline vtkSmartPointer<vtkDataArray> GetVectors()
  { return Result()->GetPointData()->GetVectors();  }
  /*!
    \fn Model::GetNormals()
    \brief Returns the normal vectors of the model.
	*/
  inline vtkSmartPointer<vtkDataArray> GetNormals()
  { return Result()->GetPointData()->GetNormals();  }
  /*!
    \fn Model::GetTensors()
    \brief Returns the active tensors (3x3 matrix per vertex) of the model.
	*/
  inline vtkSmartPointer<vtkDataArray> GetTensors()
  { return Result()->GetPointData()->GetTensors();  }
  /*!
    \fn Model::Update()
    \brief Updates the model to be the most current.
	*/
  inline void Update()
  {
  #if VTK_MAJOR_VERSION > 5
    CurrentModel->Modified();
  #else
    CurrentModel->Update();
  #endif // VTK_MAJOR_VERSION
  }
  //@}

  /**
  * \defgroup Atomic Atomic Operations Members
  * \brief Members for atomic (one model at a time) operations
  */
  //@{
  /*!
    \fn Model::Clean()
    \brief Removes duplicate points etc. from model.
	*/
  void Clean();
  /*!
    \fn Model::Triangulate()
    \brief Ensures the polygons in the model are triangulated. Useful when some filters required triangulation before use.
	*/
  void Triangulate();
  /*!
    \fn Model::MatchInformation(vtkSmartPointer<vtkPolyData> matchMesh, const bool rescale)
    \brief Matches the centroid and centroid scale of the model to the one given.
	*/
  void MatchInformation(vtkSmartPointer<vtkPolyData> matchMesh, const bool rescale);
  /*!
    \fn Model::IterativeClosestPointsRegistration(vtkSmartPointer<vtkPolyData> fixedMesh, const bool similarity = false, const int maxPoints = 1024)
    \brief Computes the ICP registration of the current mesh to the one provided.
	*/
  void IterativeClosestPointsRegistration(vtkSmartPointer<vtkPolyData> fixedMesh, const bool similarity = false, const int maxPoints = 1024);
  /*!
    \fn Model::LandmarkBasedRegistration(vtkSmartPointer<vtkPolyData> fixedMesh, const bool similarity = false)
    \brief Computes the landmark transform registration of the current mesh to the one provided. Requires point correspondence and same number of points.
	*/
  void LandmarkBasedRegistration(vtkSmartPointer<vtkPolyData> fixedMesh, const bool similarity = false);
  /*!
    \fn Model::Decimate(const float reduceFactor)
    \brief Decimate the mesh (in terms of points) using the Decimate PRO. Generally use QuadricDecimate() instead.
	*/
  void Decimate(const float reduceFactor);
  /*!
    \fn Model::QuadricDecimate(const float reduceFactor)
    \brief Decimate the mesh (in terms of points) using the Quadric algorithm.
	*/
  void QuadricDecimate(const float reduceFactor);
  /*!
    \fn Model::ClusterDecimate()
    \brief Decimate the mesh (in terms of points) using the Quadric Clustering algorithm.

    This is a one-step mesh simiplification operation. See QuadricDecimate() for similar approach with more control.
	*/
  void ClusterDecimate();
  /*!
    \fn Model::LaplacianSmoothing(const int iterations)
    \brief Smooth the mesh (in terms of points) using the Laplacian algorithm. It is unstable, use WindowedSincSmoothing() instead unless massive non-edge preserving smoothing is required.
	*/
  void LaplacianSmoothing(const int iterations);
  /*!
    \fn Model::WindowedSincSmoothing(const int iterations, const float passBand = 0.1)
    \brief Smooth the mesh (in terms of points) using the Taubin (1995) optimal filter.
	*/
  void WindowedSincSmoothing(const int iterations, const float passBand = 0.1);
  /*!
    \fn Model::Curvature(const bool meanCurvature = true)
    \brief Compute the curvature (using Laplace-Beltrami operator) of the mesh.
	*/
  void Curvature(const bool meanCurvature = true);
  /*!
    \fn Model::Gradient()
    \brief Compute the gradient (using Laplace-Beltrami operator) of the scalars on the mesh.
	*/
  void Gradient();
  /*!
    \fn Model::DistanceMap(vtkSmartPointer<vtkPoints> loopPoints)
    \brief Compute the distance map of the mesh points given a loop (as points) assumed to be on the mesh surface.
	*/
  void DistanceMap(vtkSmartPointer<vtkPoints> loopPoints);
  /*!
    \fn Model::GaussianWeightScalars(float stddev, float clampValue)
    \brief Evaluate the Gaussian function with the given standard deviation for all values of the scalars on the mesh.

    Clamp value is the lowest allowed value, so that the scalars will never have a value below this value.
        */
  void GaussianWeightScalars(float stddev, float clampValue);
  /*!
    \fn Model::Rotate(const bool xAxis, const bool yAxis, const bool zAxis, const float angle, const coordinate centre)
    \brief Rotate the mesh along the selected axes from the centre by angle.
	*/
  void Rotate(const bool xAxis, const bool yAxis, const bool zAxis, const float angle, const coordinate centre);
  /*!
    \fn Model::Flip(const bool xAxis, const bool yAxis, const bool zAxis)
    \brief Flip the mesh along the selected axes.
	*/
  void Flip(const bool xAxis, const bool yAxis, const bool zAxis);
   /*!
    \fn Model::Clip(const coordinateType belowValue, const coordinateType aboveValue)
    \brief Clip the mesh based on a scalar threshold of the mesh scalars.
	*/
  void Clip(const coordinateType belowValue, const coordinateType aboveValue);
  /*!
    \fn Model::Threshold(const coordinateType belowVal, const coordinateType aboveVal, const coordinateType outsideVal)
    \brief Threshold the scalars on the mesh based on a levels provided.
	*/
  void Threshold(const coordinateType belowVal, const coordinateType aboveVal, const coordinateType outsideVal);
  /*!
    \fn Model::Threshold(const coordinateType belowVal, const coordinateType aboveVal, const coordinateType insideVal, const coordinateType outsideVal)
    \brief Threshold the scalars on the mesh based on a levels provided. Inside value means a binary threshold is effectively done.
	*/
  void Threshold(const coordinateType belowVal, const coordinateType aboveVal, const coordinateType insideVal, const coordinateType outsideVal);
  /*!
    \fn Model::Mask(vtkSmartPointer<vtkPolyData> maskMesh)
    \brief Mask the scalars on the mesh based on values of maskMesh being >= 1.
	*/
  void Mask(vtkSmartPointer<vtkPolyData> maskMesh);
  /** Iso surface uses marching cubes to generate surface. The minValue is essentially a threshold below exclusive of its value.*/
  void IsoSurface(vtkSmartPointer<vtkImageData> img, double minValue);
  /** Voxelises mode. Returns NULL is failed.*/
  vtkSmartPointer<vtkImageData> Voxelise(const unsigned char insideValue, double *spacing, double *bounds = NULL, const size_t padVoxels = 1);
  /*!
    \fn Model::DelaunayGraph3D()
    \brief Compute the Delaunay 3D graph of the points in the mesh.
	*/
  void DelaunayGraph3D();
  /*!
    \fn Model::Delaunay2DTriangulation()
    \brief Compute the Delaunay 2D triangluation of the points in the mesh.
	*/
  void Delaunay2DTriangulation();
  /*!
    \fn Model::DelaunayTriangulation()
    \brief Compute the Delaunay 3D triangluation of the points in the mesh.
	*/
  void DelaunayTriangulation();
  /*!
    \fn Model::Append(vtkSmartPointer<vtkPolyData> model)
    \brief Appends the model provided to the member into the current model
	*/
  void Append(vtkSmartPointer<vtkPolyData> model);
  /*!
    \fn Model::Scale(vtkSmartPointer<vtkPolyData> model, const coordinateType scale, const bool useCentroid = true)
    \brief Scales the coordinates of each point in the model provided by scale

  useCentroid enables if the scaling is to be done wrt the centroid of the model.
  Note: The scale is just applied to the model provided and is not copied internally like the filters!
  Note: Does not keep track of Previous Result.
	*/
  static void Scale(vtkSmartPointer<vtkPolyData> model, const coordinateType scale, const bool useCentroid = true);
  /*!
    \fn Model::Scale(vtkSmartPointer<vtkPolyData> model, const coordinateType scale, const bool useCentroid = true)
    \brief Scales the coordinates of each point in the internal model provided by scale

  useCentroid enables if the scaling is to be done wrt the centroid of the model.
  Note: Does not keep track of Previous Result.
	*/
  inline void Scale(const coordinateType scale, const bool useCentroid = true)
  { Scale(CurrentModel, scale, useCentroid); }
  //@}

  /**
  * \defgroup Generate Generate Display Members
  * \brief Members for generating certain types of things from the models
  */
  //@{
  /*!
    \fn Model::GenerateVertices()
    \brief Generate (only) vertices for display from point data. Good for when only points in mesh.
	*/
  void GenerateVertices();
  /*!
    \fn Model::GenerateVertexScalars()
    \brief Generate vertex scalars for display from point data. This colours the mesh by point ids.
	*/
  void GenerateVertexScalars();
  /*!
    \fn Model::GenerateVerticesAs2DGlyphs(const GlyphType glyphType)
    \brief Generate (only) vertices (as 2D glyphs like crosses etc.) for display from point data. Good for when only points in mesh.
	*/
  void GenerateVerticesAs2DGlyphs(const GlyphType glyphType);
  /*!
    \fn Model::GenerateTubes()
    \brief Generate tubes along edges of the mesh.
	*/
  void GenerateTubes();
  /*!
    \fn Model::GeneratePointModel(const double newScale)
    \brief Generate spheres of scaling for each point of the mesh.
	*/
  void GeneratePointModel(const double newScale);
  /*!
    \fn Model::GenerateSpheres(const double newScale)
    \brief Generate spheres of scaling for each point of the mesh. Same as GeneratePointModel().
	*/
  inline void GenerateSpheres(const double newScale)
  { GeneratePointModel(newScale); }
  /*!
    \fn Model::GenerateSampledPoints(const double distance)
    \brief Generate sampled points (only points) for the mesh at given spacing on triangles.
	*/
  void GenerateSampledPoints(const double distance);
  /*!
    \fn Model::GenerateVectorField(const double newScale = 0.0, const bool useNormals = false)
    \brief Generate vector field for vectors present (as an array) in the mesh at given scaling.

    If useNormals is true, then normals will be shown. Set scale to zero (default) if auto scaling needed.
	*/
  void GenerateVectorField(const double newScale = 0.0, const bool useNormals = false);
  /*!
    \fn Model::GenerateTensorField(const double newScale = 0.0)
    \brief Generate tensor field for tensors present (as an array) in the mesh at given scaling.

    Set scale to zero (default) if auto scaling needed.
	*/
  void GenerateTensorField(const double newScale = 0.0);
  /*!
    \fn Model::GenerateHedgehog(const double newScale = 0.0)
    \brief Generate a line field for vectors present (as an array) in the mesh at given scaling.

    Set scale to zero (default) if auto scaling needed.
	*/
  void GenerateHedgehog(const double newScale = 0.0);
  /*!
    \fn Model::GenerateNormals(const int pointNormals = 0)
    \brief Generate normal vectors for each point on the mesh. If pointNormals is (0,1,2), then normals are for (points,cells,both).
	*/
  void GenerateNormals(const int pointNormals = 0);
  /*!
    \fn Model::GenerateKMeansClustering(int numberOfClusters)
    \brief Generate k-means clustering for the point positions in the mesh.
	*/
  void GenerateKMeansClustering(int numberOfClusters);
  /*!
    \fn Model::GenerateQuantisedPoints(float quantiseFactor)
    \brief Generate points of the mesh to have integral cordinates.
	*/
  void GenerateQuantisedPoints(float quantiseFactor);
  /*!
    \fn Model::GenerateCappedBoundaries()
    \brief Generate capped boundaries for open meshes. This closes off open ends of a mesh.
	*/
  void GenerateCappedBoundaries();
  /*!
    \fn Model::GenerateRegions()
    \brief Generate connected regions for the mesh labelled by scalar values.
	*/
  void GenerateRegions();
  /*!
    \fn Model::GenerateElevation(double x, double y, double z)
    \brief Generate scalar field for the mesh proportional to the elevation (dot product) with a given vector/direction.
	*/
  void GenerateElevation(double x, double y, double z);
  /*!
    \fn Model::GenerateElevation()
    \brief Generate scalar field for the mesh proportional to the elevation in the z-direction.
	*/
  void GenerateElevation();
  /*!
    \fn Model::GenerateReebGraph()
    \brief Generate Reeb Graph from scalar field on mesh (uses elevation if no field present).
	*/
  //~ void GenerateReebGraph();
  //@}

  /**
  * \defgroup Scalar Scalar operations
  * \brief Members for operating of model scalars/arrays, modifies internal polydata
  */
  //@{
  /*!
    \fn Model::ResetScalars()
    \brief Resets the scalars to zero of the current model

  Note: Does not keep track of Previous Result.
	*/
  void ResetScalars();
  /*!
    \fn Model::FillScalars(const coordinateType value)
    \brief Sets all the scalars of the current model to value

  Note: Does not keep track of Previous Result.
	*/
  void FillScalars(const coordinateType value);
  /*!
    \fn Model::ClearScalars(vtkSmartPointer<vtkPolyData> &mesh)
    \brief Removes all the scalars from the current model

  Note: Does not keep track of Previous Result. Use RemoveScalars() instead.
	*/
  static void ClearScalars(vtkSmartPointer<vtkPolyData> &mesh);
  /*!
    \fn Model::ClearArrays(vtkSmartPointer<vtkPolyData> &mesh)
    \brief Removes all the arrays from the current model

  Note: Does not keep track of Previous Result. Use RemoveArrays() instead.
	*/
  static inline void ClearArrays(vtkSmartPointer<vtkPolyData> &mesh)
  {
  #if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION < 8)
      mesh->GetPointData()->SetActiveScalars(NULL);
  #else
    ///Strip existing scalars/vectors
    for(int j = 0; j < mesh->GetPointData()->GetNumberOfArrays(); j ++)
      mesh->GetPointData()->RemoveArray(0);
  #endif
  }
  /*!
    \fn Model::RemoveScalars()
    \brief Removes all the scalars from the current model. Same as ClearScalars()
	*/
  void RemoveScalars();
  /*!
    \fn Model::RemoveArrays()
    \brief Removes all the arrays from the current model. Same as ClearArrayss()
	*/
  void RemoveArrays();
  /*!
    \fn Model::ScalarDifference(vtkSmartPointer<vtkPolyData> model1, vtkSmartPointer<vtkPolyData> model2, vtkSmartPointer<vtkPolyData> result)
    \brief Calculates differences within scalars of two meshes assuming the meshes have the same number of points.

  Note: Does not keep track of Previous Result.
	*/
	static void ScalarDifference(vtkSmartPointer<vtkPolyData> model1, vtkSmartPointer<vtkPolyData> model2, vtkSmartPointer<vtkPolyData> result);
	/*!
    \fn Model::ThresholdScalars(vtkSmartPointer<vtkPolyData> model, const coordinateType aboveVal, const coordinateType belowVal, const coordinateType outsideVal)
    \brief Thresholds scalar values of mesh in place.

  Note: Does not keep track of Previous Result.
	*/
	static void ThresholdScalars(vtkSmartPointer<vtkPolyData> model, const coordinateType aboveVal, const coordinateType belowVal, const coordinateType outsideVal);
	/*!
    \fn Model::BinaryThresholdScalars(vtkSmartPointer<vtkPolyData> model, const coordinateType aboveVal, const belowVal, coordinateType insideVal = 1.0, const coordinateType outsideVal = 0.0)
    \brief Thresholds scalar values of mesh in place. Scalars outside the given range is set to outsideVal and scalars inside the range set to insideVal.

  Note: Does not keep track of Previous Result. Setting insideVal = outsideVal retains scalar values within the range.
	*/
	static void BinaryThresholdScalars(vtkSmartPointer<vtkPolyData> model, const coordinateType aboveVal, const coordinateType belowVal, const coordinateType insideVal = 1.0, const coordinateType outsideVal = 0.0);
	/*!
    \fn Model::MaskScalars(vtkSmartPointer<vtkPolyData> model1, vtkSmartPointer<vtkPolyData> model2)
    \brief Masks the scalars in model1 with scalars in model2. Only values corresponding to >= 1 in model2 are kept in model1.
	*/
	static void MaskScalars(vtkSmartPointer<vtkPolyData> model1, vtkSmartPointer<vtkPolyData> model2);
  //@}

  /**
  * \defgroup Collection Collection operations
  * \brief Members for operating of a collection of models
  */
  //@{
  /*!
    \fn Model::ProcrustesAlign(vtkSmartPointer<vtkPolyData> source, vtkSmartPointer<vtkPolyData> target, bool rigid = false)
    \brief Aligns the source mesh to target mesh assuming the meshes have the same number of points and correspondence. The transform is applied in-place to the source mesh.

  The transformation matrix is returned, it is NULL if failure was encountered.
  The meshes are aligned in a least squares sense. The member returns the aligned mesh and is internally allocated.
  The scalars are preserved from the source. The rigid argument makes the alignment rigid body, by default similarity is used instead.
	*/
  static vtkSmartPointer<vtkMatrix4x4> ProcrustesAlign(vtkSmartPointer<vtkPolyData> source, vtkSmartPointer<vtkPolyData> target, bool rigid = false);

  //Batch (collection) operations
  /*!
    \fn Model::ConcatenateCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Concatenates all the models in the collection together into one model.

  This member calls Append() over the members of the collection to give the concatenated result.
  The order of the models in the collection will effect the order of the point IDs of the
  resultant model.
	*/
  void ConcatenateCollection(vtkSmartPointer<vtkPolyDataCollection> collection);
  /*!
    \fn Model::ScaleCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType scale)
    \brief Scales the coordinates of each point in all the models in the collection.

  The collection is modified in place.
	*/
  void ScaleCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType scale);
  /*!
    \fn Model::SmoothCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t iterations)
    \brief Smoothes all points in all the models in the collection using the Windowed Sinc algorithm.

  The collection is modified in place.
	*/
  void SmoothCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t iterations);
  /*!
    \fn Model::LaplacianCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t iterations)
    \brief Smoothes all points in all the models in the collection using the Laplacian algorithm.

  The collection is modified in place.
    */
  void LaplacianCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t iterations);
  /*!
    \fn Model::FlipCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const bool xAxis, const bool yAxis, const bool zAxis)
    \brief Flips all points in all the models in the collection using the axis provided.

  The collection is modified in place.
    */
  void FlipCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const bool xAxis, const bool yAxis, const bool zAxis);
  /*!
    \fn Model::DecimateCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType factor)
    \brief Decimates all points in all the models in the collection using the Quadric algorithm algorithm.

  The collection is modified in place.
	*/
  void DecimateCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType factor);
  /*!
    \fn Model::SplitCollection(vtkSmartPointer<vtkPolyDataCollection> collection, vtkSmartPointer<vtkPolyDataCollection> components, std::vector< vtkSmartPointer<vtkPolyDataCollection> > &splitCollections)
    \brief Splits each of the models in the collection into a collection of models.

  The order of the models in the components needs to match the order of the models in the collection.
	*/
  void SplitCollection(vtkSmartPointer<vtkPolyDataCollection> collection, vtkSmartPointer<vtkPolyDataCollection> components, std::vector< vtkSmartPointer<vtkPolyDataCollection> > &splitCollections);
  /*!
    \fn Model::ScalarDifferenceCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Computes the difference in scalars cummulatively over the collection wrt the first mesh assuming the meshes have the same number of points.

  The result is effectively a sum of differences across the collection with respect to the first model.
	*/
  void ScalarDifferenceCollection(vtkSmartPointer<vtkPolyDataCollection> collection);
  void ScalarDifferenceCollection(vtkSmartPointer<vtkPolyDataCollection> collection1, vtkSmartPointer<vtkPolyDataCollection> collection2, vtkSmartPointer<vtkPolyDataCollection> &results);
  /*!
    \fn Model::ScalarThresholdCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Computes the Threshold of scalar values of meshes in place over the collection wrt the first mesh assuming the meshes have the same number of points.

  The collection is modified in place.
	*/
  void ScalarThresholdCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType aboveVal, const coordinateType belowVal);
  /*!
    \fn Model::ScalarStatisticsCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Computes the statistics (min, max, variance etc.) of the scalars over the collection assuming the meshes have the same number of points.

    Use Result() to get result, which is populated by replicating the first mesh with arrays for each statistic.
	*/
  void ScalarStatisticsCollection(vtkSmartPointer<vtkPolyDataCollection> collection);
  /*!
    \fn Model::ScalarRemoveCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Removes the scalars over the collection assuming the meshes have the same number of points.

    The collection is modified in place.
	*/
  void ScalarRemoveCollection(vtkSmartPointer<vtkPolyDataCollection> collection);
  /*!
    \fn Model::ScalarCopyCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Copies the scalars over the collection (from the first mesh) assuming the meshes have the same number of points.

    The collection is modified in place.
	*/
  void ScalarCopyCollection(vtkSmartPointer<vtkPolyDataCollection> collection);
  /*!
    \fn Model::MeanSquaredErrorCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Computes the MSE of all models in collection accummulatively.
	*/
  double MeanSquaredErrorCollection(vtkSmartPointer<vtkPolyDataCollection> collection);
  /*!
    \fn Model::ProcrustesAlignCollection(vtkSmartPointer<vtkPolyDataCollection> collection, bool rigid = false)
    \brief Aligns the collection to the mean mesh (computed internally) of the collection assuming the meshes have the same number of points and correspondence.

  The entire collection is then registered to one another in a least squares sense.
  The member returns a collection of aligned meshes and the mesh internally stored is the mean mesh.
  The mean mesh scalars are the mean of the scalars of all the meshes.
  The rigid argument makes the alignment rigid body, by default similarity is used instead.
	*/
  vtkSmartPointer<vtkPolyDataCollection> ProcrustesAlignCollection(vtkSmartPointer<vtkPolyDataCollection> collection, bool rigid = false);
  /*!
    \fn Model::IterativeClosestPointsAlignCollection(vtkSmartPointer<vtkPolyDataCollection> collection, bool rigid = false)
    \brief Aligns the collection to the mean mesh (computed internally) of the collection without needing correspondence.

  The entire collection is then registered to the first mesh in a least squares sense.
  The member returns a collection of aligned meshes and the mesh internally stored is the mean mesh.
  The rigid argument makes the alignment rigid body, by default similarity is used instead.
	*/
  vtkSmartPointer<vtkPolyDataCollection> IterativeClosestPointsAlignCollection(vtkSmartPointer<vtkPolyDataCollection> collection, bool rigid = false, vtkSmartPointer<vtkTransformCollection> tCollection=0);
  //@}

protected:
  vtkSmartPointer<vtkPolyData> InputModel; //!< Holds the initial model in the pipeline
  vtkSmartPointer<vtkPolyData> CurrentModel; //!< Holds the current model in the pipeline
  vtkSmartPointer<vtkPolyData> PreviousModel; //!< Holds the previous model in the pipeline
  vtkSmartPointer<vtkTransform> CurrentTransform; //!< transform applied to model

  bool InternalInPlaceOperation; //!< Used by the collection members to assign via pointers rather than deep copys

  inline bool IsCurrentModel()
  {
    if(!CurrentModel)
    {
      PrintError("Input Model not set. Ignoring Operation.");
      return false;
    }

    return true;
  }
  inline void UpdateModelState(vtkSmartPointer<vtkPolyDataAlgorithm> filter)
  {
    CurrentModel = filter->GetOutput();
    PreviousModel = filter->GetPolyDataInput(0);
    if(!InputModel)
      InputModel = PreviousModel;
  }
  inline void UpdateModelStateFromModel(vtkSmartPointer<vtkPolyData> mesh)
  {
    PreviousModel = CurrentModel;
    CurrentModel = mesh;
    if(!InputModel)
      InputModel = PreviousModel;
  }

private:

};

} //end namespace milx

#endif //__MILXMODEL_H
