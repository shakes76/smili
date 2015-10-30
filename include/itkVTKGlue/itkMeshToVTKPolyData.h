#ifndef __MeshToVTKPolyData_h__
#define __MeshToVTKPolyData_h__

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include "itkPoint.h"
#include "itkObject.h"


namespace itk
{
  
/** 
  \class MeshToVTKPolyData
  \brief 
    \warning
  \sa 
  */

template <class TMesh >
class MeshToVTKPolyData : public Object
{

 public:

  /** Standard class typedefs. */
  typedef MeshToVTKPolyData       Self;
  typedef Object             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshToVTKPolyData, Object);

  typedef TMesh TriangleMeshType;
  typedef typename TriangleMeshType::MeshTraits TriangleMeshTraits;
  typedef typename TriangleMeshType::PointType                       PointType;
  typedef typename TriangleMeshType::PointsContainer                 InputPointsContainer;
  typedef typename InputPointsContainer::Pointer            InputPointsContainerPointer;
  typedef typename InputPointsContainer::Iterator           InputPointsContainerIterator;
  typedef typename TriangleMeshType::CellType                        CellType; 
  
  typedef typename TriangleMeshType::CellsContainerPointer           CellsContainerPointer;
  typedef typename TriangleMeshType::CellsContainerIterator          CellsContainerIterator;

  /**
  The SetInput method provides pointer to the vtkPolyData
  */
  void SetInput(TriangleMeshType * mesh);
  TriangleMeshType * GetInput();

  vtkSmartPointer<vtkPolyData> GetOutput();

  void Update();

 private:
  MeshToVTKPolyData( void );
  virtual ~MeshToVTKPolyData( void );
  MeshToVTKPolyData(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename TriangleMeshType::Pointer m_itkTriangleMesh;

  vtkSmartPointer<vtkPoints>  m_Points;
  vtkSmartPointer<vtkPolyData> m_PolyData;
  vtkSmartPointer<vtkCellArray> m_Polys;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshToVTKPolyData.txx"
#endif

#endif
