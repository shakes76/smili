/*=========================================================================
  Copyright: (c) CSIRO, Australia

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.
=========================================================================*/

#include "vtkAreaSimplificationMetric.h"

//STL
#include <map>
//VTK
#include <vtkTriangle.h>

vtkStandardNewMacro(vtkAreaSimplificationMetric);

double vtkAreaSimplificationMetric::ComputeMetric(vtkDataSet *mesh,
  vtkDataArray *scalarField, vtkIdType startCriticalPoint,
  vtkAbstractArray* vertexList, vtkIdType endCriticalPoint)
{
  // In this example, the metric algorithm just evaluates the area of the
  // surface region corresponding to the arc of the Reeb graph passed as an
  // argument.
  // As a result, the arcs corresponding to small surface regions (below the
  // threshold specified to the simplificatin filter) will be
  // simplified in priority in the surface simplification algorithm.

  double  fieldLowerBound = scalarField->GetComponent(startCriticalPoint,0),
          fieldUpperBound = scalarField->GetComponent(endCriticalPoint,0);

  double  cumulativeArea = 0;

  std::map<vtkIdType, bool> visitedTriangles;

  for(int i = 0; i < vertexList->GetNumberOfTuples(); i++)
    {
    int vId = vertexList->GetVariantValue(i).ToInt();
    vtkIdList *starTriangleList = vtkIdList::New();

    mesh->GetPointCells(vId, starTriangleList);

    for(int j = 0; j < starTriangleList->GetNumberOfIds(); j++)
      {
      vtkIdType tId = starTriangleList->GetId(j);
      vtkTriangle *t = vtkTriangle::SafeDownCast(mesh->GetCell(tId));
      std::map<vtkIdType, bool>::iterator tIt = visitedTriangles.find(tId);
      if(tIt == visitedTriangles.end())
        {
        if((scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
          <= fieldUpperBound)
          &&(scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
          <= fieldUpperBound)
          &&(scalarField->GetComponent(t->GetPointIds()->GetId(2), 0)
          <= fieldUpperBound)
          &&(scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
          >= fieldLowerBound)
          &&(scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
          >= fieldLowerBound)
          &&(scalarField->GetComponent(t->GetPointIds()->GetId(2), 0)
          >= fieldLowerBound))
          {
          // the triangle fully maps inside the arc function interval
          cumulativeArea += t->ComputeArea();
          }
        visitedTriangles[tId] = true;
        }
      }

    starTriangleList->Delete();
    }

  return cumulativeArea/(this->UpperBound - this->LowerBound);
}
