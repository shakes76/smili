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
#ifndef __MILXMATH_H
#define __MILXMATH_H

#ifndef VTK_ONLY
  //VNL
  #include <vnl/vnl_error.h>
  #include <vnl/vnl_matrix.h>
  #include <vnl/vnl_vector_fixed.h>
  #include <vnl/algo/vnl_matrix_inverse.h>
#endif
#ifndef ITK_ONLY
  //VTK
  #include <vtkMath.h>
  #include <vtkPoints.h>
#endif
//SMILI
#include <milxGlobal.h>

namespace milx
{

/**
  \class Math
  \brief A general math object for computing common operations such a centroids etc.
*/
template<class Type = double>
class SMILI_EXPORT Math
{
public:
#ifndef VTK_ONLY
  /**
   * \defgroup Stats Statistical Methods
   * \brief Common statistical functions
   */
  //@{
  //! \brief Compute the centroid (or the mean) of a data vector
  static Type Centroid(const vnl_vector<Type> &data);
  //! \brief Compute the centroid (or the mean vector) of a data matrix
  ///
  /// Assumes that the rows contain the observations and the
  /// columns contain the variables
  static vnl_vector<Type> Centroid(const vnl_matrix<Type> &data);
  //! \brief Compute the centroid size or scale for a series of points
  ///
  /// norm allows the scale to be average scale per point
  static Type CentroidSize(const vnl_matrix<Type> &data, const vnl_vector<Type> &centroid, bool norm = false);
#ifndef ITK_ONLY
  //! \brief Compute the centroid (or the mean vector) of a data matrix
  ///
  /// Assumes that the data matrix rows contain the observations and the
  /// columns contain the variables
  /// Overloaded for points, which are common in vtkPolyData
  static vnl_vector<Type> Centroid(vtkPoints *points);
  //! \brief Compute the centroid size or scale for a series of points
  ///
  /// Overloaded for points, which are common in vtkPolyData
  /// norm allows the scale to be average scale per point
  static Type CentroidSize(vtkPoints *points, const vnl_vector<Type> &centroid, bool norm = false);
#endif
#endif
#ifdef VTK_ONLY
  //! \brief Compute the centroid (or the mean vector) of a data matrix
  ///
  /// Assumes that the data matrix rows contain the observations and the
  /// columns contain the variables
  /// Overloaded for points, which are common in vtkPolyData
  static Type* Centroid(vtkPoints *points);
  //! \brief Compute the centroid size or scale for a series of points
  ///
  /// Overloaded for points, which are common in vtkPolyData
  /// norm allows the scale to be average scale per point
  static Type CentroidSize(vtkPoints *points, const Type* centroid, bool norm = false);
#endif

#ifndef VTK_ONLY
  //! \brief Compute the covariance matrix from a vector. EXPERIMENTAL
  ///
  /// Use with caution, this member computes the inter-component covariance matrix
  /// for the given vector
  static vnl_matrix<Type> CovarianceMatrix(vnl_vector<Type> sourceVector);
  //! \brief Compute the covariance matrix between two vectors. EXPERIMENTAL
  ///
  /// Use with caution, this member computes the inter-component covariance matrix
  /// for the given vectors
  static vnl_matrix<Type> CovarianceMatrix(vnl_vector<Type> sourceVector, vnl_vector<Type> targetVector);
  //! \brief Compute the covariance matrix from a datamatrix
  ///
  /// Assumes that the data matrix rows contain the observations and the
  /// columns contain the variables. Centroid is removed from each row.
  /// The result is a nxn matrix where n is the number of variables
  /// The diagonal should contain the variances of the variables
  static vnl_matrix<Type> CovarianceMatrix(const vnl_matrix<Type> &data, const vnl_vector<Type> &centroid);
  //! \brief Compute the covariance matrix from a datamatrix
  ///
  /// Assumes that the data matrix rows contain the observations and the
  /// columns contain the variables. Centroid is assumed to have removed already.
  /// The result is a nxn matrix where n is the number of variables
  /// The diagonal should contain the variances of the variables
  static vnl_matrix<Type> CovarianceMatrix(const vnl_matrix<Type> &data);
#ifndef ITK_ONLY
  //! \brief Compute the covariance matrix from a datamatrix
  ///
  /// Assumes that the data matrix rows contain the observations and the
  /// columns contain the variables. Centroid is removed from each row.
  /// The result is a nxn matrix where n is the number of variables
  /// The diagonal should contain the variances of the variables
  /// Overloaded for points, which are common in vtkPolyData
  static vnl_matrix<Type> CovarianceMatrix(vtkPoints *points, const vnl_vector<Type> &centroid);
#endif
#endif

#ifndef VTK_ONLY
  //! \brief Compute the Mahalanobis distance between two vectors with an unknown covariance matrix. EXPERIMENTAL
  ///
  /// Mahalanobis distance takes into account the covariance between vectors in
  /// estimating the separation between them. Returns squared distance.
  /// This version assumes you do not know the covariance matrix, so an approximate one
  /// is constructed for matching purposes
  /// Vectors should be of the same length
  static Type MahalanobisDistance(const vnl_vector<Type> &source, const vnl_vector<Type> &target);
  //! \brief Compute the Mahalanobis distance of a vector to a population. EXPERIMENTAL
  ///
  /// Mahalanobis distance takes into account the covariance between the population and the vector in
  /// estimating the separation between them, so its more robust to outlyers. Returns squared distance.
  /// This version assumes you know the inverse covariance matrix of a population.
  /// The result of this member depends solely on the order of your covariance matrix
  static Type MahalanobisDistance(const vnl_vector<Type> &target, const vnl_vector<Type> &mean, const vnl_matrix<Type> &invCovMatrix);
  //! \brief Compute the Mahalanobis distance of a vector to a population.  EXPERIMENTAL
  ///
  /// Mahalanobis distance takes into account the covariance between the population and the vector in
  /// estimating the separation between them, so its more robust to outlyers. Returns squared distance.
  /// This version assumes you know the covariance matrix of a population and the inverse
  /// covariance matrix is assigned to invCovMatrix
  /// The result of this member depends solely on the order of your covariance matrix
  static Type MahalanobisDistance(const vnl_vector<Type> &target, const vnl_vector<Type> &mean, const vnl_matrix<Type> &covMatrix, vnl_matrix<Type> &invCovMatrix);
  //@}
#endif

  /**
   * \defgroup Measures Measure Methods
   * \brief Common measurement functions such as MSE etc.
   */
  //@{
  //! \brief Compute the MSE of the given points sets.
  ///
  /// The MSE is computed across all points between the two point sets.
  static Type MeanSquaredError(vtkPoints *sourcePoints, vtkPoints *targetPoints);
  //@}

protected:

private:
  Math() {};
  ~Math() {};

};

#ifndef VTK_ONLY
template<class Type>
Type Math<Type>::Centroid(const vnl_vector<Type> &data)
{
  const vnl_size_t n = data.size();
  Type centroid = 0.0;

  for (vnl_size_t j = 0; j < n; j ++)
    centroid += data(j); ///Sum
  centroid /= n; ///then average

  return centroid;
}

template<class Type>
vnl_vector<Type> Math<Type>::Centroid(const vnl_matrix<Type> &data)
{
  const vnl_size_t n = data.rows();
  const vnl_size_t m = data.cols();
  vnl_vector<Type> centroid(m, 0.0);

  for (vnl_size_t j = 0; j < n; j ++)
  {
    centroid += data.get_row(j); ///Sum
  }
  centroid /= n; ///then average

  return centroid;
}

template<class Type>
Type Math<Type>::CentroidSize(const vnl_matrix<Type> &data, const vnl_vector<Type> &centroid, bool norm)
{
  const vnl_size_t n = data.rows();
  Type newScale = 0;

  for (int j = 0; j < n; j ++)
  {
    vnl_vector<Type> rowVector( data.get_row(j) );

    ///Accumulate the squared distance from centroid
    newScale += vtkMath::Distance2BetweenPoints(rowVector.data_block(), centroid.data_block());
  }
  
  if(norm)
    newScale /= n;

  return sqrt(newScale); ///Norm and return
}

#ifndef ITK_ONLY
template<class Type>
vnl_vector<Type> Math<Type>::Centroid(vtkPoints *points)
{
  const vtkIdType n = points->GetNumberOfPoints();
  vnl_vector_fixed<Type, 3> centroid(0.0);

  for (vtkIdType j = 0; j < n; j ++)
  {
    vnl_vector_fixed<Type, 3> location( points->GetPoint(j) );

    centroid += location; ///Sum
  }
  centroid /= n; ///then average

  return centroid.as_vector();
}

template<class Type>
Type Math<Type>::CentroidSize(vtkPoints *points, const vnl_vector<Type> &centroid, bool norm)
{
  Type newScale = 0;

  for (vtkIdType j = 0; j < points->GetNumberOfPoints(); j ++)
  {
    vnl_vector_fixed<Type, 3> location( points->GetPoint(j) );

    ///Accumulate the squared distance from centroid
    newScale += vtkMath::Distance2BetweenPoints(location.data_block(), centroid.data_block());
  }
  
  if(norm)
    newScale /= points->GetNumberOfPoints();

  return sqrt(newScale); ///Norm and return
}
#endif
#endif

#ifdef VTK_ONLY //Code that does not depend on VNL
template<class Type>
Type* Math<Type>::Centroid(vtkPoints *points)
{
  const vtkIdType n = points->GetNumberOfPoints();
  //~ Type centroid[3] = {0.0, 0.0, 0.0};
  Type *centroid = new Type[3];

  for (vtkIdType j = 0; j < n; j ++)
  {
    Type location[3];
    
    points->GetPoint(j, location);

    centroid[0] += location[0]/n; ///Sum
    centroid[1] += location[1]/n; ///Sum
    centroid[2] += location[2]/n; ///Sum
  }

  return centroid;
}

template<class Type>
Type Math<Type>::CentroidSize(vtkPoints *points, const Type* centroid, bool norm)
{
  Type newScale = 0;

  for (vtkIdType j = 0; j < points->GetNumberOfPoints(); j ++)
  {
    Type location[3];
    
    points->GetPoint(j, location);

    ///Accumulate the squared distance from centroid
    newScale += vtkMath::Distance2BetweenPoints(location, centroid);
  }
  
  if(norm)
    newScale /= points->GetNumberOfPoints();

  return sqrt(newScale); ///Norm and return
}
#endif

#ifndef VTK_ONLY
template<class Type>
vnl_matrix<Type> Math<Type>::CovarianceMatrix(vnl_vector<Type> sourceVector)
{
  const vnl_size_t n = sourceVector.size();

  sourceVector -= sourceVector.mean();

  vnl_matrix<Type> covarianceMatrix(n, n, 0.0);
  covarianceMatrix += outer_product(sourceVector, sourceVector); ///take outer product and add to matrix

  return covarianceMatrix;
}

template<class Type>
vnl_matrix<Type> Math<Type>::CovarianceMatrix(vnl_vector<Type> sourceVector, vnl_vector<Type> targetVector)
{
  const vnl_size_t n = sourceVector.size();
  const vnl_size_t m = targetVector.size();

  if (m != n)
  {
    vnl_error_vector_dimension ("Math<Type>::CovarianceMatrix", m, n);
  }

  sourceVector -= sourceVector.mean();
  targetVector -= targetVector.mean();

  vnl_matrix<Type> covarianceMatrix(n, n, 0.0);
  covarianceMatrix += outer_product(sourceVector, targetVector); ///take outer product and add to matrix

  return covarianceMatrix;
}

template<class Type>
vnl_matrix<Type> Math<Type>::CovarianceMatrix(const vnl_matrix<Type> &data, const vnl_vector<Type> &centroid)
{
  const vnl_size_t n = data.rows();
  const vnl_size_t m = data.cols();
  vnl_size_t norm = 1;

  if (n > 1)
    norm = n-1;

  if (m != centroid.size())
  {
    vnl_error_vector_dimension ("Math<Type>::CovarianceMatrix", m, centroid.size());
  }

  vnl_matrix<Type> covarianceMatrix(m, m, 0.0);
  for (vnl_size_t j = 0; j < n; j ++)
  {
    vnl_vector<Type> rowVector(data.get_row(j));

    rowVector -= centroid; ///Subtract centroid

    covarianceMatrix += outer_product(rowVector, rowVector) / norm; ///take outer product and add to matrix
  }

  return covarianceMatrix;
}

template<class Type>
vnl_matrix<Type> Math<Type>::CovarianceMatrix(const vnl_matrix<Type> &data)
{
  const vnl_size_t n = data.rows();
  const vnl_size_t m = data.cols();
  vnl_size_t norm = 1;

  if (n > 1)
    norm = n-1;

  vnl_matrix<Type> covarianceMatrix(m, m, 0.0);
  for (vnl_size_t j = 0; j < n; j ++)
  {
    vnl_vector<Type> rowVector(data.get_row(j));

    covarianceMatrix += outer_product(rowVector, rowVector) / norm; ///take outer product and add to matrix
  }

  return covarianceMatrix;
}

#ifndef ITK_ONLY
template<class Type>
vnl_matrix<Type> Math<Type>::CovarianceMatrix(vtkPoints *points, const vnl_vector<Type> &centroid)
{
  const vnl_size_t n = points->GetNumberOfPoints();
  vnl_size_t norm = 1;

  if (n > 1)
    norm = n-1;

  if (centroid.size() != 3)
  {
    vnl_error_vector_dimension ("Math<Type>::CovarianceMatrix", 3, centroid.size());
  }

  vnl_matrix<Type> covarianceMatrix(3, 3, 0.0);
  for (vnl_size_t j = 0; j < n; j ++)
  {
    vnl_vector_fixed<Type, 3> location( points->GetPoint(j) );

    location -= centroid; ///Subtract centroid

    covarianceMatrix += outer_product(location.as_vector(), location.as_vector()) / norm; ///take outer product and add to matrix
  }

  return covarianceMatrix;
}
#endif
#endif

#ifndef VTK_ONLY
template<class Type>
Type Math<Type>::MahalanobisDistance(const vnl_vector<Type> &source, const vnl_vector<Type> &target)
{
  const vnl_size_t n = source.size();
  const vnl_size_t m = target.size();

  if (n != m)
  {
    //Throw exception
    vnl_error_vector_dimension ("Math<Type>::MahalanobisDistance", n, m);
  }

  ///Type 1 Computation (row wise)
  /*vnl_matrix<Type> dataMatrix(2, n, 0.0);
  dataMatrix.set_row(0, target);
  dataMatrix.set_row(1, source);

  //	cout << "Compute Centroid" << endl;
  const vnl_vector<Type> dataCentroid = milx::Math<Type>::Centroid(dataMatrix);
    cout << "Centroid: " << dataCentroid << endl;
  //	cout << "Compute Covariance Matrix" << endl;
  vnl_matrix<Type> dataCovariance = milx::Math<Type>::CovarianceMatrix(dataMatrix, dataCentroid);
  cout << "Covariance: " << dataCovariance.rows() << "x" << dataCovariance.cols() << endl;
    cout << "Covariance: \n" << dataCovariance << endl;*/

  ///Type 3 (column wise)
  vnl_matrix<Type> dataMatrix(n, 2, 0.0);
  dataMatrix.set_column(0, target);
  dataMatrix.set_column(1, source);

  const vnl_vector<Type> dataCentroid = Math<Type>::Centroid(dataMatrix);
  //cout << "Centroid: " << dataCentroid << endl;
  vnl_matrix<Type> dataCovariance = Math<Type>::CovarianceMatrix(dataMatrix, dataCentroid);
  //cout << "Cov: " << dataCovariance.rows() << "x" << dataCovariance.cols() << endl;
  //cout << "Covariance: \n" << dataCovariance << endl;

  ///Type 4 (Cov Mat is 1x1)
  /*vnl_matrix<Type> dataMatrix(m, 1, 0.0);
  dataMatrix.set_column(0, target);

  //	const Type centroid = milx::Math<Type>::Centroid(target);
  //    const vnl_vector<Type> dataCentroid(m, centroid);
    const vnl_vector<Type> dataCentroid = milx::Math<Type>::Centroid(dataMatrix);
    cout << "Centroid: " << dataCentroid << endl;
  vnl_matrix<Type> dataCovariance = milx::Math<Type>::CovarianceMatrix(dataMatrix, dataCentroid);
  cout << "Cov: " << dataCovariance.rows() << "x" << dataCovariance.cols() << endl;
    cout << "Covariance: \n" << dataCovariance << endl;*/

  ///Type 5
  /*vnl_matrix<Type> dataMatrix(1, m, 0.0);
  dataMatrix.set_row(0, target);

  //    const vnl_vector<Type> dataCentroid = milx::Math<Type>::Centroid(dataMatrix);
  //    const Type centroid = milx::Math<Type>::Centroid(target);
  //    const vnl_vector<Type> dataCentroid(m, centroid);
  const vnl_vector<Type> dataCentroid(m, 0.0);
  cout << "Centroid: " << source << endl;
  vnl_matrix<Type> dataCovariance = milx::Math<Type>::CovarianceMatrix(dataMatrix, dataCentroid);
  cout << "Cov: " << dataCovariance.rows() << "x" << dataCovariance.cols() << endl;
  dataCovariance.set_identity();
  cout << "Covariance: \n" << dataCovariance << endl;*/

  ///Type 6
//    vnl_matrix<Type> dataCovariance = milx::Math<Type>::CovarianceMatrix(source, target);

  ///Type 7
//    vnl_matrix<Type> dataCovariance = milx::Math<Type>::CovarianceMatrix(target);

//	vnl_matrix<Type> dataDifference(n, 1, 0.0);
  vnl_matrix<Type> dataDifference = dataMatrix;
//    const vnl_vector<Type> sourceZeroMean = source - source.mean();
//    const vnl_vector<Type> targetZeroMean = target - target.mean();
//    dataDifference.set_column(0, targetZeroMean - sourceZeroMean);
//    dataDifference.set_column(0, target - source);

  //cout << "Compute Inverse Covariance Matrix" << endl;
  vnl_matrix<Type> dataInvCovariance = vnl_matrix_inverse<Type>(dataCovariance);
  //cout << "Inv Covariance: \n" << dataInvCovariance << endl;

//	const vnl_matrix<Type> distance = dataDifference.transpose()*dataInvCovariance*dataDifference;
  const vnl_matrix<Type> distance = dataDifference*dataInvCovariance*dataDifference.transpose();

  //cout << "Distance Size: " << distance.rows() << "x" << distance.cols() << endl;
  //cout << "Distance: " << distance << endl;
  return 1.0/distance(0,0); //cause of covariance, maximum correlation is returned so invert
}

template<class Type>
Type Math<Type>::MahalanobisDistance(const vnl_vector<Type> &target, const vnl_vector<Type> &mean, const vnl_matrix<Type> &invCovMatrix)
{
  const vnl_size_t n = target.size();
  const vnl_size_t m = mean.size();

  if (n != m)
  {
    //Throw exception
    vnl_error_vector_dimension ("Math<Type>::MahalanobisDistance", n, m);
  }

  vnl_matrix<Type> dataDifference(n, 1, 0.0);
  dataDifference.set_column(0, target - mean);

  const vnl_matrix<Type> distance = dataDifference.transpose()*invCovMatrix*dataDifference;

  return distance(0,0);
}

template<class Type>
Type Math<Type>::MahalanobisDistance(const vnl_vector<Type> &target, const vnl_vector<Type> &mean, const vnl_matrix<Type> &covMatrix, vnl_matrix<Type> &invCovMatrix)
{
  vnl_matrix<Type> dataInvCovariance = vnl_matrix_inverse<Type>(covMatrix);

  invCovMatrix = dataInvCovariance;

  return MahalanobisDistance(target, mean, invCovMatrix);
}

#ifndef ITK_ONLY
template<class Type>
Type Math<Type>::MeanSquaredError(vtkPoints *sourcePoints, vtkPoints *targetPoints)
{
  const vnl_size_t numberOfPoints = sourcePoints->GetNumberOfPoints();
  Type mse = 0.0;

  for (vnl_size_t v = 0; v < numberOfPoints; v ++) ///For each point
  {
    coordinate p( sourcePoints->GetPoint(v) );
    coordinate p2( targetPoints->GetPoint(v) );

    ///Accumulate distance squared for each point across all sets
    mse += vtkMath::Distance2BetweenPoints(p.data_block(), p2.data_block());
  }
  mse /= numberOfPoints;

  return mse;
}
#endif
#endif

#ifdef VTK_ONLY //Code that does not depend on VNL
template<class Type>
Type Math<Type>::MeanSquaredError(vtkPoints *sourcePoints, vtkPoints *targetPoints)
{
  const size_t numberOfPoints = sourcePoints->GetNumberOfPoints();
  Type mse = 0.0;

  for (size_t v = 0; v < numberOfPoints; v ++) ///For each point
  {
    Type p[3], p2[3];
    sourcePoints->GetPoint(v, p);
    targetPoints->GetPoint(v, p2);

    ///Accumulate distance squared for each point across all sets
    mse += vtkMath::Distance2BetweenPoints(p, p2)/numberOfPoints;
  }

  return mse;
}
#endif

} //end namespace milx

#endif //__MILXMATH_H
