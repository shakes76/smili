/*=========================================================================
  Program: MILX MSK
  Module: WPCA Class
  Author: Shekhar Chandra
  Modified by:
  Language: C++
  Created: Tue 19/04/2011

  Copyright: (c) 2009-11 CSIRO, Australia.

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/
/**
  \class vtkWeightedPCAAnalysisFilter
  \brief Performs weighted principal component analysis (WPCA) of a set of aligned pointsets

  The EM algorithm is used to construct a subspace that spans the weighted PCA space.
  The weights here are spatial weights and are to be assigned per point.
  The eigenvalues and eigenvectors of the spanned subspace is then computed using SVD.

  See the following journal for more details:
  Skočaj, D.; Leonardis, A. & Bischof, H.
  Weighted and robust learning of subspace representations
  Pattern Recogn., Elsevier Science Inc., 2007, 40, 1556-1569

  Note about weighted mean and unweighted mean: The weighted mean and the unweighted mean
  are the same if each dimension is equally weighted.

  Note about the initial guess: By default, the standard PCA is not computed first where atleast one
  iteration of the EM algorithm would be computed. The system is initialised with random values.
  The algorithm continues until tolerance is reached or the maximum number of iterations is exceeded.
  Change this behaviour using RandomGuessOff().

  Usage Example:
  \code
  vtkSmartPointer<vtkWeightedPCAAnalysisFilter> wpca = vtkSmartPointer<vtkWeightedPCAAnalysisFilter>::New();
      wpca->SetNumberOfInputs(4);
      wpca->SetConvergenceTolerance(1e-2); //1e-3
      wpca->SetMaximumIterations(20); //100
      wpca->SetInput(0, procrustes->GetOutput(0), weights);
      wpca->SetInput(1, procrustes->GetOutput(1), lesserWeights);
      wpca->SetInput(2, procrustes->GetOutput(2), lesserWeights);
      wpca->SetInput(3, procrustes->GetOutput(3), lesserWeights);
      ...
      wpca->Update();
  \endcode
*/

#ifndef __vtkWeightedPCAAnalysisFilter_h
#define __vtkWeightedPCAAnalysisFilter_h

#include "vtkPointSetAlgorithm.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "milxWin32Header.h"

class vtkPointSet;

class MILX_EXPORT vtkWeightedPCAAnalysisFilter : public vtkPointSetAlgorithm
{
public:
  vtkTypeMacro(vtkWeightedPCAAnalysisFilter,vtkPointSetAlgorithm);

  // Description:
  //! Prints information about the state of the filter.
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  //! Creates with similarity transform.
  static vtkWeightedPCAAnalysisFilter *New();

  // Description:
  //! Get/Set the eigenvalues sorted in descending order
  vtkSetObjectMacro(EigenValues, vtkFloatArray);
  vtkGetObjectMacro(EigenValues, vtkFloatArray);
  inline vtkFloatArray* GetEvals()
  { return this->GetEigenValues(); }

  // Description:
  //! Get/Set the eigenvectors sorted in descending order
  vtkSetMacro(EigenVectors, vnl_matrix<double>);
  vtkGetMacro(EigenVectors, vnl_matrix<double>);

  // Description:
  //! Get/Set the eigenvectors sorted in descending order
  vtkSetMacro(InverseWeightedEigenVectors, vnl_matrix<double>);
  vtkGetMacro(InverseWeightedEigenVectors, vnl_matrix<double>);

  // Description:
  //! Get/Set the data matrix used for the WPCA
  vtkSetMacro(DataMatrix, vnl_matrix<double>);
  vtkGetMacro(DataMatrix, vnl_matrix<double>);

  // Description:
  //! Get/Set the weighted mean shape
  vtkSetMacro(MeanShape, vnl_vector<double>);
  vtkGetMacro(MeanShape, vnl_vector<double>);

  // Description:
  //! Get/Set the weights for each point in the training set
  vtkSetMacro(Weights, vnl_matrix<double>);
  vtkGetMacro(Weights, vnl_matrix<double>);

  // Description:
  //! Get/Set the precomputed inverse weights for each point in the training set
  vtkSetMacro(InverseWeights, vnl_matrix<double>);
  vtkGetMacro(InverseWeights, vnl_matrix<double>);

  // Description:
  //! Get/Set the reconstructed training sets used to compute the eigen-decomposition.
  vtkGetMacro(Reconstruction, vnl_matrix<double>);

  // Description:
  //! Stop when difference between actual and reconstructed results are less than this value (default is 1e-3)
  vtkSetMacro(ConvergenceTolerance, float);
  vtkGetMacro(ConvergenceTolerance, float);

  // Description:
  //! Set/Get the precision of the modes to use in recovering missing data
  vtkSetMacro(Precision, float);
  vtkGetMacro(Precision, float);

  // Description:
  //! Set/Get the sum of the eigenvalues found
  vtkSetMacro(SumOfEigenValues, double);
  vtkGetMacro(SumOfEigenValues, double);

  // Description:
  //! Stop if this number of iterations is reached (default is 3)
  vtkSetMacro(MaximumIterations, size_t);
  vtkGetMacro(MaximumIterations, size_t);

  // Description:
  //! Number of iterations to use for the iterative reconstruction algorithm for recovering missing data (default is MaximumIterations)
  //! After Generating the model, it is the number of iterations that EM algorithm converged with.
  vtkSetMacro(ReconstructionIterations, size_t);
  vtkGetMacro(ReconstructionIterations, size_t);

  //!Get/Set whether to use the iterative reconstruction EM algorithm for missing data or to use the more general weighted PCA
  //!Default On.
  vtkSetMacro(MissingMode, bool);
  vtkGetMacro(MissingMode, bool);
  vtkBooleanMacro(MissingMode, bool);

  //!Get/Set whether eigenvalues and eigenvectors are already known or pre-loaded, if so, don't bother doing EM algorithm
  //!Use stored result.
  vtkSetMacro(PreLoaded, bool);
  vtkGetMacro(PreLoaded, bool);
  vtkBooleanMacro(PreLoaded, bool);

  //!Get/Set whether to use random values as the initialisation to the EM algorithm when building the Weighted PCA solution
  //!Else use the standard PCA as initialisation. Applies to missing Mode OFF only and Default On
  vtkSetMacro(RandomGuess, bool);
  vtkGetMacro(RandomGuess, bool);
  vtkBooleanMacro(RandomGuess, bool);

  // Description:
  //! Returns the final reconstruction error found for the current eign-decomposition
  vtkGetMacro(ReconstructionError, float);

  // Description:
  //! Specify how many pointsets are going to be given as input.
  void SetNumberOfInputs(int n);

  // Description:
  //! Specify the input pointset with index idx and weights.
  //! Weights are assigned to each point as a float so its size should
  //! be equal to the number of points of the point set
  //! Call SetNumberOfInputs before calling this function.
  void SetInput(int idx, vtkPointSet* p); //Assumes weights assigned elsewhere
  void SetInput(int idx, vtkPointSet* p, vtkFloatArray *weights);
  void SetInput(int idx, vtkDataObject* input, vtkFloatArray *weights);

  // Description:
  //! Retrieve the input with index idx (usually only used for pipeline
  //! tracing).
  vtkPointSet* GetInput(int idx);

  /**
    \fn vtkWeightedPCAAnalysisFilter::ParameteriseShape(vtkPointSet* shape, vtkFloatArray *weights, const int k)
    \brief Converts the weighted shape into a parameterised version according to the components held internally.

    Weighted Shape is replaced with a parameterised version.
  */
  void ParameteriseShape(vtkPointSet* shape, vtkFloatArray *weights, const int k);

  // Description:
  //! Fills the shape with:
  //
  // mean + b[0] * sqrt(eigenvalue[0]) * eigenvector[0]
  //      + b[1] * sqrt(eigenvalue[1]) * eigenvector[1]
  // ...
  //      + b[sizeb-1] * sqrt(eigenvalue[bsize-1]) * eigenvector[bsize-1]
  //
  // here b are the parameters expressed in standard deviations
  // bsize is the number of parameters in the b vector
  // This function assumes that shape is allready allocated
  // with the right size, it just moves the points.
  void GetParameterisedShape(vtkFloatArray *b, vtkPointSet* shape);
  void GetParameterisedShape(vtkFloatArray *b, vtkPointSet* shape, vtkFloatArray *weights);

  // Description:
  //! Return the bsize parameters b that best model the given weighted shape
  //! The appropriate algorithm is used given Missing Mode.
  // (in standard deviations).
  // That is that the given shape will be approximated by:
  //
  // shape ~ mean + b[0] * sqrt(eigenvalue[0]) * eigenvector[0]
  //              + b[1] * sqrt(eigenvalue[1]) * eigenvector[1]
  //         ...
  //              + b[bsize-1] * sqrt(eigenvalue[bsize-1]) * eigenvector[bsize-1]
  void GetShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize);
  void GetMissingShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize);
  void GetStandardShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize);
  void GetWeightedShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize);

  // Description:
  //! Retrieve how many modes are necessary to model the given proportion of the variation.
  //! proportion should be between 0 and 1
  int GetModesRequiredFor(double proportion);

  //! Save or Load the Weighted PCA solution from/to file.
  bool Save(const char * filename);
  bool Load(const char * filename);

  //! Member for applying the PCA of matrices provided using SVD
  void ApplyStandardPCA(const vnl_matrix<double> &data, vnl_matrix<double> &eigenVecs, vnl_vector<double> &eigenVals);

protected:
  vtkWeightedPCAAnalysisFilter();
  ~vtkWeightedPCAAnalysisFilter();

  // Description:
  //! Usual data generation method.
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  //! Request data for when missing mode is on
  virtual int RequestMissingData();
  //! Request data for when standard mode
  virtual int RequestStandardData();
  //! Request data for when missing mode is off
  virtual int RequestWeightedData();
  //! Request data for when missing mode is off and weights are temporal (same per training set)
  virtual int RequestTemporallyWeightedData();
  virtual int FillInputPortInformation(int port, vtkInformation *info);
  /**
    \fn vtkWeightedPCAAnalysisFilter::ProjectShape(vtkPointSet *shape, int bsize, vtkFloatArray *weights = 0)
    \brief Project the weighted shape into the eigenspace so its represented by the components.
    \return Return the components representing the shape as a sx1 matrix.
  */
  virtual vnl_vector<double> ProjectShape(vtkPointSet *shape, int bsize, vtkFloatArray *weights = 0);
  /**
    \fn vtkWeightedPCAAnalysisFilter::BackProjectShape(vnl_vector<double> components, vtkPointSet *shape, const int k, bool unweightedShape = false)
    \brief Re-Project the components (using k components) into real space so its represented as a shape.

    The shape is replaced with the new shape.
  */
  virtual void BackProjectShape(vnl_vector<double> components, vtkPointSet *shape, const int k);
  virtual void BackProjectWeightedShape(vnl_vector<double> components, vtkPointSet *shape, vtkFloatArray *weights, const int k);
  //! Convert input basis {v} to orthonormal basis {e}. Utilised for missing mode off
  // Each basis vector is a column of v, and likewise the orthonormal bases are returned as columns of e
  void Gram_Schmidt_Orth(const vnl_matrix<double>& v, vnl_matrix<double>& e);

protected:
  vtkWeightedPCAAnalysisFilter(const vtkWeightedPCAAnalysisFilter&);  // Not implemented.
  void operator=(const vtkWeightedPCAAnalysisFilter&);  // Not implemented.

  float ConvergenceTolerance; // Stop when difference between actual and reconstructed results are less than this
  float Precision; // Precision of the missing data iteratively reconstruction
  size_t MaximumIterations; // Stop if this number of iterations is reached
  size_t ReconstructionIterations; // Stop if this number of iterations is reached
  bool MissingMode; //!< Use the iterative reconstruction algorithm for missing data?
  bool PreLoaded; //!< Solution is preloaded? Such as from file.
  bool RandomGuess; //!< Use random guess in the EM algorithm to build the WPCA solution?
  float ReconstructionError; //!< Final Reconstruction error of the eigen-decomposition
  double SumOfEigenValues;

  // Eigenvalues
  vtkFloatArray *EigenValues;

  // Matrix where each column is an eigenvector
  vnl_matrix<double> EigenVectors;
  vnl_matrix<double> InverseWeightedEigenVectors;

  // The mean shape in a vector
  vnl_vector<double> MeanShape;
  vnl_vector<double> TotalMissing;
  vnl_matrix<double> DataMatrix;
  vnl_matrix<double> Coefficients; //!< Reconstruction coefficients for training set
  vnl_matrix<double> Reconstruction; //!< reconstruction of the data matrix used in missing mode on
  vnl_matrix<double> Weights; //!< Weights for the analysis
  vnl_matrix<double> WeightMatrix;
  vnl_matrix<double> InverseWeights; //!< 1/Weights for the analysis

  //! Member for computing the reconstruction from the eigen-decomposition of the data matrix.
  double Reconstruct(vnl_matrix<double> &recon, const vnl_matrix<double> &data, const vnl_matrix<double> &eigenVecs, const vnl_matrix<double> &coefficients, const vnl_matrix<double> &means);
  void printProgBar( int percent );
};

#endif // header guard, vtkWeightedPCAAnalysisFilter
