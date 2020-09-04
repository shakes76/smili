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
#include <time.h>

#include "vtkWeightedPCAAnalysisFilter.h"
#include "vtkSmartPointer.h"
#include "vtkExecutive.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPolyData.h"
#include "vtkMath.h"
#include "vtkFloatArray.h"
#include "vtkCommand.h"

#include <vcl_vector.h>
#include "vnl/vnl_random.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_qr.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

vtkStandardNewMacro(vtkWeightedPCAAnalysisFilter);

//----------------------------------------------------------------------------
// protected
vtkWeightedPCAAnalysisFilter::vtkWeightedPCAAnalysisFilter()
{
  this->EigenValues = vtkFloatArray::New();
  ConvergenceTolerance = 1e-9;
  Precision = 0.95;
  MaximumIterations = 3;
  ReconstructionIterations = MaximumIterations;
  MissingMode = true;
  PreLoaded = false;
  RandomGuess = true;
  ReconstructionError = 1e9;
  SumOfEigenValues = 0.0;
}

//----------------------------------------------------------------------------
// protected
vtkWeightedPCAAnalysisFilter::~vtkWeightedPCAAnalysisFilter()
{
  if (this->EigenValues)
  {
    this->EigenValues->Delete();
  }
}

//----------------------------------------------------------------------------
// protected
int vtkWeightedPCAAnalysisFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  ///get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  ///get the input and output
  vtkPointSet *input = vtkPointSet::SafeDownCast(
                         inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet *output = vtkPointSet::SafeDownCast(
                          outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDebugMacro(<<"Execute()");

  const int N_SETS = this->GetNumberOfInputConnections(0);

  vtkInformation *tmpInfo;
  vtkPointSet *tmpInput;
  ///copy the inputs across
  output->DeepCopy(input);
  for (int i = 1; i < N_SETS; i ++)
  {
    tmpInfo = inputVector[0]->GetInformationObject(i);
    tmpInput = 0;
    if (tmpInfo)
    {
      tmpInput = vtkPointSet::SafeDownCast(tmpInfo->Get(vtkDataObject::DATA_OBJECT()));
    }
    this->GetOutput(i)->DeepCopy(tmpInput);
  }

  ///the number of points is determined by the first input (they must all be the same)
  const int N_POINTS = input->GetNumberOfPoints();

  vtkDebugMacro(<<"N_POINTS is " <<N_POINTS);

  if (N_POINTS == 0)
  {
    vtkErrorMacro(<<"No points!");
    return 0;
  }

  ///all the inputs must have the same number of points to consider executing
  for (int i = 1; i < N_SETS; i ++)
  {
    tmpInfo = inputVector[0]->GetInformationObject(i);
    tmpInput = 0;
    if (tmpInfo)
    {
      tmpInput =
        vtkPointSet::SafeDownCast(tmpInfo->Get(vtkDataObject::DATA_OBJECT()));
    }
    else
    {
      continue;
    }
    if (tmpInput->GetNumberOfPoints() != N_POINTS)
    {
      vtkErrorMacro(<<"The inputs have different numbers of points!");
      return 0;
    }
  }

  ///Check if preloaded, if so sone and return
  if (MissingMode && PreLoaded)
  {
    vtkDebugMacro(<<"Robust PCA is Preloaded. Aborting Generation.");
    cerr << "Robust PCA is Preloaded. Aborting Generation." << endl;
    Coefficients = InverseWeightedEigenVectors*DataMatrix;
    Reconstruction = EigenVectors*Coefficients;
    return 1;
  }

  ///Number of shapes
  const int s = N_SETS;

  ///Number of points in a shape
  const int n = N_POINTS;

  ///Observation Matrix [number of points * 3 X number of shapes]
  DataMatrix.clear();
  DataMatrix.set_size(3*n, s); //Matrix X
  DataMatrix.fill(0.0);

  //~ cerr << "Creating " << 3*n << "x" << s << " Data Matrix" << endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < s; j++)
    {
      tmpInfo = inputVector[0]->GetInformationObject(j);
      tmpInput = 0;
      if (tmpInfo)
      {
        tmpInput = vtkPointSet::SafeDownCast(
                     tmpInfo->Get(vtkDataObject::DATA_OBJECT()));
      }
      else
      {
        continue;
      }
      double p[3];
      tmpInput->GetPoint(i, p);
      DataMatrix(i*3,   j) = p[0];
      DataMatrix(i*3+1, j) = p[1];
      DataMatrix(i*3+2, j) = p[2];
    }
  }

  ///Check if preloaded, if so sone and return
  if (!MissingMode && PreLoaded)
  {
    vtkDebugMacro(<<"Weighted PCA is Preloaded. Aborting Generation.");
    cerr << "Weighted PCA is Preloaded. Aborting Generation." << endl;
//    InverseWeightedEigenVectors = vnl_matrix_inverse<double>(EigenVectors).pinverse();
    Coefficients = InverseWeightedEigenVectors*DataMatrix;
    Reconstruction = EigenVectors*Coefficients;
    return 1;
  }

  if(MissingMode)
    return RequestMissingData();
//    return RequestStandardData();

//  return RequestWeightedData();
  return RequestTemporallyWeightedData();
}

int vtkWeightedPCAAnalysisFilter::RequestMissingData()
{
  const size_t n = this->GetOutput(0)->GetNumberOfPoints();
  const size_t s = this->GetNumberOfInputConnections(0);

  ///The mean shape is also calculated
  MeanShape.clear();
  MeanShape.set_size(3*n);
  MeanShape.fill(0);
  TotalMissing.clear();
  TotalMissing.set_size(3*n);
  TotalMissing.fill(0);

  for(size_t i = 0; i < s; i++)
    {
    MeanShape += element_product(DataMatrix.get_column(i), Weights.get_column(i));
    TotalMissing += Weights.get_column(i);
    }
  MeanShape = element_quotient(MeanShape, TotalMissing);

  const vnl_vector<double> onesCol(s, 1.0);
  const vnl_matrix<double> MeansMatrix = outer_product(MeanShape, onesCol);

  DataMatrix -= MeansMatrix;

  ///Form Weights Matrix
  WeightMatrix.clear();
  WeightMatrix.set_size(3*n, s);
  WeightMatrix.fill(0);
  for(unsigned i = 0; i < s; i++)
    WeightMatrix.set_column(i, Weights.get_column(i));

  //Check if there is any missing data
  bool noMissingData = false;
//  if(WeightMatrix.min_value() == 1.0)
//  {
//    cout << "No Missing Data in Training Set. Not using iterative reconstruction for model." << endl;
//    noMissingData = true;
//  }

  ///Initialisation from PCA of non-missing pixels
  DataMatrix = element_product(DataMatrix, WeightMatrix);

  ///PCA
  vnl_vector<double> eigenValues(s, 1.0);
  this->EigenValues->SetNumberOfValues(s);
  this->EigenValues->FillComponent(0, 0.0);
  this->EigenVectors.clear();
  this->EigenVectors.set_size(3*n, s); //matrix U
  this->EigenVectors.fill(0.0);
  Coefficients.clear();
  Coefficients.set_size(s, s); //matrix A
  Coefficients.fill(1.0);
  InverseWeightedEigenVectors.clear();
  Reconstruction.clear();

  ApplyStandardPCA(DataMatrix, this->EigenVectors, eigenValues);
  for (size_t j = 0; j < s; j++)
    EigenValues->SetValue(j, eigenValues[j]); //!< Extract eigenvalues after normalising by s-1

  //Net error output csv file
  ofstream outFile("WPCA_Net_Error.csv", ios::app);

  vtkDebugMacro(<< "Robust PCA with tolerance " << ConvergenceTolerance << " and max iterations " << MaximumIterations);
  double error = 1e9, previousError = 0.0;
  unsigned iteration = 0;
  while(error > ConvergenceTolerance && iteration < MaximumIterations && !noMissingData)
  {
    std::cerr << "Iteration " << iteration << " ------ " << std::endl;

    ///Solve RPCA
    vnl_matrix<double> WeightedEigenVectors = element_product(WeightMatrix, this->EigenVectors);
//    WeightedEigenVectors.normalize_columns();

    ///E-Step
    //~ std::cerr << "Computing E-Step" << std::endl;
    InverseWeightedEigenVectors = vnl_matrix_inverse<double>(WeightedEigenVectors).pinverse();
    Coefficients = InverseWeightedEigenVectors*DataMatrix;

    ///M-step
    //~ std::cerr << "Computing M-Step" << std::endl;
    ///Reconstruct
    double currentError = Reconstruct(Reconstruction, DataMatrix, this->EigenVectors, Coefficients, MeansMatrix);

    if(iteration == 0)
      error = fabs(currentError);
    else
      error = fabs( previousError - currentError );

    previousError = fabs(currentError);
    std::cerr << "Net error: " << error << std::endl;
    outFile << iteration << ", " << currentError << ", " << error << std::endl;

    DataMatrix = element_product(DataMatrix, WeightMatrix);

    //Replace missing pixels with reconstructed values
    for(size_t j = 0; j < Reconstruction.rows(); j ++)
      for(size_t k = 0; k < Reconstruction.cols(); k ++)
      {
        if(WeightMatrix(j,k) == 0) //missing pixels
          DataMatrix(j,k) = Reconstruction(j,k) - MeansMatrix(j,k);
      }

    iteration ++;
    ReconstructionIterations = iteration;

    ///PCA
    ApplyStandardPCA(DataMatrix, this->EigenVectors, eigenValues);
    for (size_t j = 0; j < s; j++)
      EigenValues->SetValue(j, eigenValues[j]); //!< Extract eigenvalues

    InvokeEvent(vtkCommand::ProgressEvent);
  }
  MaximumIterations = ReconstructionIterations; //remember this number for next time
//  cout << "Iterations Used: " << ReconstructionIterations << endl;

  if(noMissingData)
  {
    InverseWeightedEigenVectors = EigenVectors.transpose(); //!< Used for projection, naming convention broken here to avoid branch of code for missing mode
    Coefficients = InverseWeightedEigenVectors*DataMatrix;
    Reconstruct(Reconstruction, DataMatrix, this->EigenVectors, Coefficients, MeansMatrix);
  }

  SumOfEigenValues = 0.0;
  for (int i = 0; i < this->EigenValues->GetNumberOfTuples(); i ++)
  {
    SumOfEigenValues += this->EigenValues->GetValue(i);
  }
  cout << "Sum of Eigenvalues: " << SumOfEigenValues << endl;
  cout << "Robust PCA Complete" << endl;
  outFile.close();

  return 1;
}

int vtkWeightedPCAAnalysisFilter::RequestStandardData()
{
  const size_t n = this->GetOutput(0)->GetNumberOfPoints();
  const size_t s = this->GetNumberOfInputConnections(0);

  ///The mean shape is also calculated
  MeanShape.clear();
  MeanShape.set_size(3*n);
  MeanShape.fill(0);
  TotalMissing.clear();
  TotalMissing.set_size(3*n);
  TotalMissing.fill(0);

  for(size_t i = 0; i < s; i++)
    {
    MeanShape += element_product(DataMatrix.get_column(i), Weights.get_column(i));
    TotalMissing += Weights.get_column(i);
    }
  MeanShape = element_quotient(MeanShape, TotalMissing);

  const vnl_vector<double> onesCol(s, 1.0);
  const vnl_matrix<double> MeansMatrix = outer_product(MeanShape, onesCol);

  DataMatrix -= MeansMatrix;

  ///PCA
  vnl_vector<double> eigenValues(s, 1.0);
  this->EigenValues->SetNumberOfValues(s);
  this->EigenValues->FillComponent(0, 0.0);
  this->EigenVectors.clear();
  this->EigenVectors.set_size(3*n, s); //matrix U
  this->EigenVectors.fill(0.0);
  Coefficients.clear();
  Coefficients.set_size(s, s); //matrix A
  Coefficients.fill(1.0);
  InverseWeightedEigenVectors.clear();
  Reconstruction.clear();

  ApplyStandardPCA(DataMatrix, this->EigenVectors, eigenValues);
  for (size_t j = 0; j < s; j++)
    EigenValues->SetValue(j, eigenValues[j]); //!< Extract eigenvalues after normalising by s-1

  InverseWeightedEigenVectors = EigenVectors.transpose(); //!< Used for projection, naming convention broken here to avoid branch of code for missing mode
  Coefficients = InverseWeightedEigenVectors*DataMatrix;
  ///Reconstruct
  double currentError = Reconstruct(Reconstruction, DataMatrix, this->EigenVectors, Coefficients, MeansMatrix);
  std::cerr << "Error: " << currentError << std::endl;

  SumOfEigenValues = 0.0;
  for (int i = 0; i < this->EigenValues->GetNumberOfTuples(); i ++)
  {
    SumOfEigenValues += this->EigenValues->GetValue(i);
  }
  cout << "Sum of Eigenvalues: " << SumOfEigenValues << endl;
  cout << "Standard PCA Complete" << endl;

  return 1;
}

int vtkWeightedPCAAnalysisFilter::RequestWeightedData()
{
  const size_t n = this->GetOutput(0)->GetNumberOfPoints();
  const size_t s = this->GetNumberOfInputConnections(0);

  ///The mean shape is also calculated
  MeanShape.clear();
  MeanShape.set_size(3*n);
  MeanShape.fill(0.0);
  TotalMissing.clear();
  TotalMissing.set_size(3*n);
  TotalMissing.fill(0.0);

//  cerr << "Weights: " << Weights << endl;
  ///Zero mean data
  cerr << "Computing Mean" << endl;
  cerr << "Data Size: " << DataMatrix.rows() << "x" << DataMatrix.cols() << endl;
  for(size_t i = 0; i < s; i++)
    {
    MeanShape += element_product(DataMatrix.get_column(i), Weights.get_column(i));
    TotalMissing += Weights.get_column(i);
    }
  MeanShape = element_quotient(MeanShape, TotalMissing);

  ///Remove the mean
  const vnl_vector<double> onesRow(s, 1.0);
  const vnl_matrix<double> MeansMatrix = outer_product(MeanShape, onesRow);

  DataMatrix -= MeansMatrix;

  ///Setup EigenVectors
  EigenVectors.clear();
  EigenVectors.set_size(3*n, s); //matrix U
  EigenVectors.fill(0.0);
  this->EigenValues->SetNumberOfValues(s);
  ///Setup Coefficients
  Coefficients.clear();
  Coefficients.set_size(s, s); //matrix A
  Coefficients.fill(1.0);
  Reconstruction.clear();
  InverseWeightedEigenVectors.clear();

  ///Initialise PRNG
  vnl_random randGen( time(NULL) );
  if (RandomGuess)
  {
    cerr << "Using Random Initial guess" << endl;
    /// Initialise eigenvectors to random values
    for (size_t j = 0; j < 3*n; j ++)
      for (size_t k = 0; k < s; k ++)
        EigenVectors(j, k) = randGen.drand64();
  }
  else
  {
    cerr << "Using PCA as Initial guess" << endl;

    ///Use PCA as initial guess
    vnl_vector<double> eigenValues(s, 1.0);
    ApplyStandardPCA(DataMatrix, this->EigenVectors, eigenValues);

    /// Initialise eigenvectors by adding random values (normally distributed) to avoid
    /// being caught within the PCA (unweighted) solution
    /*for (size_t j = 0; j < 3*n; j ++)
    {
      float vecMean = EigenVectors.get_row(j).mean();
      for (size_t k = 0; k < s; k ++)
        EigenVectors(j, k) += vecMean*randGen.normal64();
    }*/

    ///Guess coefficients
    Coefficients = vnl_matrix_inverse<double>(EigenVectors).pinverse()*DataMatrix;
    cerr << "PCA Initial guess eigenvalues: \n" << eigenValues << endl;
  }

  ///Compute Weighted PCA
  ///Loop until convergence
  float difference = 1;
  size_t iterations = 0;

//  vtkDebugMacro(<<"Iterating to EM Solution");
  cerr << "Iterating to EM Solution" << endl;
  ///Notation 'o' means element-wise (Hadamard) product, '_' is a subscript, ^-1 means the pseudo inverse
//  vnl_vector<double> meanVector(EigenVectors.rows(), 0.0);
  const vnl_matrix<double> weightsSqrt = Weights.apply(sqrt);
//  vnl_matrix<double> zeroMeanData(3*n, s, 0.0), meanw(3*n, s, 0.0);
  vnl_matrix<double> bestEigenVectors, bestCoefficients;
  ///Compute (w_.j o x_.j)
  const vnl_matrix<double> weightedData = element_product(weightsSqrt, DataMatrix); //do this only once

  ///Compute error
  Reconstruction = EigenVectors*Coefficients + MeansMatrix; ///Reconstruct all shapes (mean centered)
  ///Compute difference of steps for convergence testing
  const vnl_matrix<double> pcaError = DataMatrix - Reconstruction; ///Reconstruction error
  difference = pcaError.frobenius_norm();
  cerr << "Initial Guess Reconstruction MSE of " << difference << endl;

  while (iterations < MaximumIterations)
  {
    /*
    ///Compute error
    meanw = outer_product(MeanShape, onesRow); //works
    //      outFile << "\nMean W: \n" << meanw << endl;
    reconstruction = EigenVectors*Coefficients; ///Reconstruct all shapes (mean centered)
    //      outFile << "\nrecon: \n" << reconstruction << endl;
    zeroMeanData = DataMatrix - meanw; ///zero mean
    //      outFile << "\nzero mean: \n" << zeroMeanData << endl;
    ///Compute difference of steps for convergence testing
    const vnl_matrix<double> errorZeroMean = zeroMeanData - reconstruction; ///Reconstruction error
    error = DataMatrix - reconstruction; ///Reconstruction error
    //      outFile << "\ndiff: \n" << errorZeroMean << endl;
    difference = errorZeroMean.frobenius_norm();
    cout << "Computed Iteration " << iterations << " with Reconstruction MSE of " << difference << endl;

    ///Compute updated weighted mean
    SubtractMeanColumn(error, Weights, MeanShape, UnweightedMeanShape);
    meanw = outer_product(MeanShape, onesRow); //works
    zeroMeanData = DataMatrix - meanw; ///update zero mean data

    //    const vnl_matrix<double> oldCoefficients = Coefficients;
    //    const vnl_matrix<double> oldEigenVectors = EigenVectors;

    ///E-Step
    //    vtkDebugMacro(<<"Computing E-Step");
    cout << "Computing E-Step" << endl;
    for (int j = 0; j < s; j ++)
    {
      ///Compose matrix
      const vnl_matrix<double> weights = outer_product(weightsSqrt.get_column(j), onesRow);
      const vnl_matrix<double> weightedEigenVectors = element_product(weights, EigenVectors);
    //      const vnl_matrix_inverse<double> eigenInverse(weightedEigenVectors.transpose()*weightedEigenVectors);
      const vnl_matrix_inverse<double> eigenInverse(EigenVectors.transpose()*weightedEigenVectors); //so to get sxs matrix to invert
    //      const vnl_matrix_inverse<double> eigenInverse(weightedEigenVectors); //so to get sxs matrix to invert
    //        cerr << "eigenInverse Size: " << eigenInverse.inverse().rows() << "x" << eigenInverse.inverse().cols() << endl;
      const vnl_matrix<double> vectors = eigenInverse.inverse()*weightedEigenVectors.transpose();
    //      const vnl_matrix<double> vectors = eigenInverse.pinverse();
      const vnl_vector<double> weightedData = element_product(weightsSqrt.get_column(j), zeroMeanData.get_column(j));

      ///Compute a
      Coefficients.set_column(j, vectors*weightedData);

      printProgBar( static_cast<float>(j)/s*100 );
    }
    printProgBar( 100 );

    ///M-step
    //    vtkDebugMacro(<<"Computing M-Step");
    cout << "Computing M-Step" << endl;
    for (int j = 0; j < 3*n; j ++)
    {
      ///Compose matrix
      const vnl_matrix<double> weights = outer_product(onesRow, weightsSqrt.get_row(j));
      const vnl_matrix<double> weightedCoeffs = element_product(weights, Coefficients);
    //      const vnl_matrix_inverse<double> coeffInverse(weightedCoeffs.transpose()*weightedCoeffs);
      const vnl_matrix_inverse<double> coeffInverse(Coefficients.transpose()*weightedCoeffs);
    //      const vnl_matrix_inverse<double> coeffInverse(weightedCoeffs);
      const vnl_matrix<double> coeffs = coeffInverse.inverse()*weightedCoeffs.transpose();
    //      const vnl_matrix<double> coeffs = coeffInverse.pinverse();
      const vnl_vector<double> weightedData = element_product(weightsSqrt.get_row(j), zeroMeanData.get_row(j));

      ///Compute u
      EigenVectors.set_row(j, weightedData*coeffs);

      printProgBar( static_cast<float>(j)/(3*n)*100 );
    }
    printProgBar( 100 );
    */
    /*
    ///E-Step
//    vtkDebugMacro(<<"Computing E-Step");
    cerr << "Computing E-Step" << endl;
    for (size_t j = 0; j < s; j ++)
    {
      ///Compose matrix
      ///Notation 'o' means element-wise (Hadamard) product, '_' is a subscript, ^-1 means the pseudo inverse
      const vnl_matrix<double> weights = outer_product(weightsSqrt.get_column(j), onesRow);
      ///Compute (w_.j o U_.j)^-1
      const vnl_matrix<double> weightedEigenVectors = vnl_matrix_inverse<double>(element_product(weights, EigenVectors)).pinverse();//Pseudo inverse, result is a sx3n matrix
      ///Compute (w_.j o x_.j)
      const vnl_vector<double> weightedData = element_product(weightsSqrt.get_column(j), DataMatrix.get_column(j));

      ///Compute a
      Coefficients.set_column(j, weightedEigenVectors*weightedData);

      printProgBar( static_cast<float>(j)/s*100 );
    }
    printProgBar( 100 );

    ///M-step
//    vtkDebugMacro(<<"Computing M-Step");
    cerr << "Computing M-Step" << endl;
    for (size_t j = 0; j < 3*n; j ++)
    {
      ///Compose matrix
      const vnl_matrix<double> weights = outer_product(onesRow, weightsSqrt.get_row(j));
      ///Compute (w_i. o x_i.)
      const vnl_vector<double> weightedData = element_product(weightsSqrt.get_row(j), DataMatrix.get_row(j));
      ///Compute (w_i. o A_i.)^-1
      const vnl_matrix<double> weightedCoefficients = vnl_matrix_inverse<double>(element_product(weights, Coefficients)).pinverse(); //Pseudo inverse

      ///Compute u
      EigenVectors.set_row(j, weightedData*weightedCoefficients);

      printProgBar( static_cast<float>(j)/(3*n)*100 );
    }
    printProgBar( 100 );
    */
    ///E-Step
//    vtkDebugMacro(<<"Computing E-Step");
    cerr << "Computing E-Step" << endl;
    for (size_t j = 0; j < s; j ++)
    {
      ///Notation 'o' means element-wise (Hadamard) product, '_' is a subscript, ^-1 means the pseudo inverse
      const vnl_matrix<double> weights = outer_product(weightsSqrt.get_column(j), onesRow);
      ///Compute (w_.j o U_.j)^-1
      const vnl_matrix<double> weightedEigenVectors = vnl_matrix_inverse<double>(element_product(weights, EigenVectors)).pinverse();//Pseudo inverse, result is a sx3n matrix

      ///Compute a
      Coefficients.set_column(j, weightedEigenVectors*weightedData.get_column(j));

      printProgBar( static_cast<float>(j)/s*100 );
    }
    printProgBar( 100 );

    ///M-step
//    vtkDebugMacro(<<"Computing M-Step");
    cerr << "Computing M-Step" << endl;
    for (size_t j = 0; j < 3*n; j ++)
    {
      ///Compose matrix
      const vnl_matrix<double> weights = outer_product(onesRow, weightsSqrt.get_row(j));
      ///Compute (w_i. o A_i.)^-1
      const vnl_matrix<double> weightedCoefficients = vnl_matrix_inverse<double>(element_product(weights, Coefficients)).pinverse(); //Pseudo inverse

      ///Compute u
      EigenVectors.set_row(j, weightedData.get_row(j)*weightedCoefficients);

      printProgBar( static_cast<float>(j)/(3*n)*100 );
    }
    printProgBar( 100 );

    ///Compute error
    Reconstruction = element_product(weightsSqrt, EigenVectors*Coefficients) + MeansMatrix; ///Reconstruct all shapes (mean centered)
    ///Compute difference of steps for convergence testing
    const vnl_matrix<double> error = weightedData - Reconstruction; ///Reconstruction error
    difference = error.frobenius_norm();
    cerr << "Computed Iteration " << iterations << " with Reconstruction MSE of " << difference << endl;

    if(difference <= ReconstructionError)
    {
      cout << "Lower Reconstruction Error found at " << iterations << ". Saving result. " << endl;
      ReconstructionError = difference;
      bestEigenVectors = EigenVectors;
      bestCoefficients = Coefficients;
    }

    if (difference < ConvergenceTolerance && iterations > 0) //If PCA guess at least do one iteration
      break;

    iterations ++;
  }
//  cout << "Converged at: " << iterations << " with MSE of " << difference << endl;
//  cout << "Coefficients: " << endl;
//  cout << Coefficients << endl;

  EigenVectors = bestEigenVectors;
  Coefficients = bestCoefficients;
  cerr << "Converged with MSE of " << ReconstructionError << endl;

  ///Contents of Eigenvectors span principal subspace and do not actually present eigenvectors and eigenvalues
  ///Orthonormalise this matrix, project the data using this subspace and compute actual eigenvectors and values using SVD.
//    ///Compute QR decomposition to get orthonormal U
//    vnl_qr<double> QR(EigenVectors); //!< Use Q matrix is orthogonal and its columns are normalised
//    vnl_matrix<double> orthoEigenVectors = QR.Q().extract(3*n, s); //!< Extract First s columns of Q (3nx3n matrix) are whats needed
  ///Compute Gram-Schmidt Orthogonalisation to get orthonormal U
  ///\note Gram-Schmidt may not be as good as SVD or QR Decomposition (see above code for QR) but it can handle a very large number of points
//  vnl_matrix<double> orthoEigenVectors(3*n, s, 0.0);
//  Gram_Schmidt_Orth(EigenVectors, orthoEigenVectors); //!< Orthogonalise U using Gram-Schmidt
//  cout << "orthoEigenVectors: " << orthoEigenVectors.rows() << "x" << orthoEigenVectors.cols()  << ", magnitude: " << orthoEigenVectors.get_column(0).magnitude() << endl;

//  const vnl_matrix<double> projectedCov = vnl_matrix_inverse<double>(orthoEigenVectors).pinverse()*DataMatrix; //!< Project data into subspace, Z = U^T X
  const vnl_matrix<double> projectedCov = vnl_matrix_inverse<double>(EigenVectors).pinverse()*DataMatrix; //!< Project data into subspace, Z = U^T X
  vnl_vector<double> eigenVals;
  vnl_matrix<double> tmpEigenVectors;
  ApplyStandardPCA(projectedCov, tmpEigenVectors, eigenVals); ///PCA on A' to get U' and eigenvalues.
//  EigenVectors = orthoEigenVectors*tmpEigenVectors; ///Rotate U for U'
  InverseWeightedEigenVectors = EigenVectors.transpose(); //!< Used for projection, naming convention broken here to avoid branch of code for missing mode
  cout << "Eigenvectors: " << EigenVectors.rows() << "x" << EigenVectors.cols() << ", magnitude: " << EigenVectors.get_column(0).magnitude() << endl;
  cout << "Eigenvalues: " << endl;
  cout << eigenVals << endl;
  //cout << "Eigenvectors: " << endl;
  //cout << EigenVectors << endl;

  // Copy data to output structures
  ///Save eigenvalues into float array
  SumOfEigenValues = 0.0;
  for (size_t j = 0; j < s; j++)
  {
    EigenValues->SetValue(j, eigenVals[j]); //!< Extract eigenvalues after normalising by s-1
    SumOfEigenValues += eigenVals[j];

    for (size_t i = 0; i < n; i++)
    {
      double x = EigenVectors(i*3  , j);
      double y = EigenVectors(i*3+1, j);
      double z = EigenVectors(i*3+2, j);

      this->GetOutput(j)->GetPoints()->SetPoint(i, x, y, z);
    }
  }
  cout << "Sum of Eigenvalues: " << SumOfEigenValues << endl;
  cout << "Weight PCA Complete" << endl;

  return 1;
}

int vtkWeightedPCAAnalysisFilter::RequestTemporallyWeightedData()
{
  const size_t n = this->GetOutput(0)->GetNumberOfPoints();
  const size_t s = this->GetNumberOfInputConnections(0);

  ///The mean shape is also calculated
  MeanShape.clear();
  MeanShape.set_size(3*n);
  MeanShape.fill(0.0);
  TotalMissing.clear();
  TotalMissing.set_size(3*n);
  TotalMissing.fill(0.0);

//  cerr << "Weights: " << Weights << endl;
  ///Zero mean data
  cerr << "Computing Mean" << endl;
  cerr << "Data Size: " << DataMatrix.rows() << "x" << DataMatrix.cols() << endl;
  for(size_t i = 0; i < s; i++)
    {
    MeanShape += element_product(DataMatrix.get_column(i), Weights.get_column(i));
    TotalMissing += Weights.get_column(i);
    }
  MeanShape = element_quotient(MeanShape, TotalMissing);

  ///Remove the mean
  const vnl_vector<double> onesRow(s, 1.0);
  const vnl_matrix<double> MeansMatrix = outer_product(MeanShape, onesRow);

  DataMatrix -= MeansMatrix;
  DataMatrix = element_product(Weights, DataMatrix); //Assuming here that the weights are the same for each training set

  ///Setup EigenVectors
  EigenVectors.clear();
  EigenVectors.set_size(3*n, s); //matrix U
  EigenVectors.fill(0.0);
  this->EigenValues->SetNumberOfValues(s);
  ///Setup Coefficients
  Coefficients.clear();
  Coefficients.set_size(s, s); //matrix A
  Coefficients.fill(1.0);
  Reconstruction.clear();
  InverseWeightedEigenVectors.clear();

  vnl_vector<double> eigenValues(s, 1.0);
  ApplyStandardPCA(DataMatrix, this->EigenVectors, eigenValues);
  InverseWeightedEigenVectors = EigenVectors.transpose(); //!< Used for projection, naming convention broken here to avoid branch of code for missing mode

  // Copy data to output structures
  ///Save eigenvalues into float array
  SumOfEigenValues = 0.0;
  for (size_t j = 0; j < s; j++)
  {
    EigenValues->SetValue(j, eigenValues[j]); //!< Extract eigenvalues
    SumOfEigenValues += eigenValues[j];

    for (size_t i = 0; i < n; i++)
    {
      double x = EigenVectors(i*3  , j);
      double y = EigenVectors(i*3+1, j);
      double z = EigenVectors(i*3+2, j);

      this->GetOutput(j)->GetPoints()->SetPoint(i, x, y, z);
    }
  }
  cout << "Sum of Eigenvalues: " << SumOfEigenValues << endl;
  cout << "Spatially Weighted PCA Complete" << endl;
  return 1;
}

void vtkWeightedPCAAnalysisFilter::ApplyStandardPCA(const vnl_matrix<double> &data, vnl_matrix<double> &eigenVecs, vnl_vector<double> &eigenVals)
{
  const double norm = 1.0/(data.cols()-1);
  const vnl_matrix<double> T = (data.transpose()*data)*norm; //D^T.D is smaller so more efficient

  //SVD
  vnl_svd<double> svd(T); //!< Form Projected Covariance matrix and compute SVD, ZZ^T
  svd.zero_out_absolute(); ///Zero out values below 1e-8 but greater than zero

  ///pinverse unnecessary?
//  eigenVecs = data*vnl_matrix_inverse<double>(svd.U()).pinverse().transpose(); //!< Extract eigenvectors from U, noting U = V^T since covariance matrix is real and symmetric
  eigenVecs = data*svd.U(); //!< Extract eigenvectors from U, noting U = V^T since covariance matrix is real and symmetric
  eigenVecs.normalize_columns();

  eigenVals = svd.W().diagonal();
}

double vtkWeightedPCAAnalysisFilter::Reconstruct(vnl_matrix<double> &recon, const vnl_matrix<double> &data, const vnl_matrix<double> &eigenVecs, const vnl_matrix<double> &coefficients, const vnl_matrix<double> &means)
{
  const int n = this->GetOutput(0)->GetNumberOfPoints();
  const int s = this->GetNumberOfInputConnections(0);
  int numberOfModes = this->GetModesRequiredFor(Precision);
  std::cerr << "Modes (" << Precision << " precision) in Reconstruction: " << numberOfModes << std::endl;

  vnl_matrix<double> truncatedVecs(3*n, s, 0.0);
  for(int j = 0; j < numberOfModes; j ++)
    truncatedVecs.set_column(j, eigenVecs.get_column(j));

  recon = truncatedVecs*coefficients + means; ///Reconstruct all shapes
//  recon = eigenVecs*coefficients + means; ///Reconstruct all shapes
  ///Compute difference of steps for convergence testing
  const vnl_matrix<double> Error = data - recon; ///Reconstruction error
  double error = Error.rms();
  std::cerr << "Reconstruction MSE of " << error << std::endl;

  return error;
}

void vtkWeightedPCAAnalysisFilter::ParameteriseShape(vtkPointSet* shape, vtkFloatArray *weights, const int k)
{
  vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
    b->SetNumberOfComponents(1);
    b->SetNumberOfTuples(k);

  //Project
  GetShapeParameters(shape, weights, b, k);

  //Backproject, missing mode already has backprojected shape so skip
  if(!MissingMode)
    GetParameterisedShape(b, shape);
}

//----------------------------------------------------------------------------
// public
void vtkWeightedPCAAnalysisFilter::GetParameterisedShape(vtkFloatArray *b, vtkPointSet* shape)
{
//  cerr << "GetParameterisedShape " << endl;
  const int bsize = b->GetNumberOfTuples();
  const int n = this->GetOutput(0)->GetNumberOfPoints();

  if (shape->GetNumberOfPoints() != n)
  {
    vtkErrorMacro(<<"Input shape does not have the correct number of points");
    return;
  }

  vnl_vector<double> coeffs(bsize, 0.0); //Column Vector

  for (int i = 0; i < bsize; i ++)
    coeffs(i) = b->GetValue(i)*sqrt(EigenValues->GetValue(i));
//    coeffs(i) = b->GetValue(i);

  BackProjectShape(coeffs, shape, bsize);
}

void vtkWeightedPCAAnalysisFilter::GetParameterisedShape(vtkFloatArray *b, vtkPointSet* shape, vtkFloatArray *weights)
{
//  cerr << "GetParameterisedShape " << endl;
  const int bsize = b->GetNumberOfTuples();
  const int n = this->GetOutput(0)->GetNumberOfPoints();

  if (shape->GetNumberOfPoints() != n)
  {
    vtkErrorMacro(<<"Input shape does not have the correct number of points");
    return;
  }

  vnl_vector<double> coeffs(bsize, 0.0); //Column Vector

  for (int i = 0; i < bsize; i ++)
    coeffs(i) = b->GetValue(i)*sqrt(EigenValues->GetValue(i));
//    coeffs(i) = b->GetValue(i);

  BackProjectWeightedShape(coeffs, shape, weights, bsize);
}

//----------------------------------------------------------------------------
// public
void vtkWeightedPCAAnalysisFilter::GetShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize)
{
  if(MissingMode)
    GetMissingShapeParameters(shape, weights, b, bsize);
//    GetStandardShapeParameters(shape, weights, b, bsize);
  else
    GetWeightedShapeParameters(shape, weights, b, bsize);
}

void vtkWeightedPCAAnalysisFilter::GetMissingShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize)
{
//  cerr << "Reconstructing with " << ReconstructionIterations << " iterations" << endl;
  ///Use iterative reconstruction to fill in missing parts and return final parameters.
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(shape->GetPoints()); //Save current shape

  /*double range[2];
  weights->GetRange(range);
  if(range[0] == 1.0) //minimum value is 1.0 then no missing data
  {
    cout << "No missing data. Not reconstructing iteratively." << endl;
    vnl_vector<double> projectedShape = ProjectShape(shape, bsize); //Unweighted since no missing pixels
    for (int j = 0; j < bsize; j++)
      b->SetValue(j, projectedShape[j]/sqrt(EigenValues->GetValue(j))); //!< Extract eigenvalues after normalising by s-1
//      b->SetValue(j, projectedShape[j]); //!< Extract eigenvalues after normalising by s-1
    BackProjectShape(projectedShape, shape, bsize); ///shape in stddevs by convention
    return; //if so, stop
  }*/

  ///Project first
  vnl_vector<double> projectedShape = ProjectShape(shape, bsize, weights); ///Result of projection is a sx1 matrix
  for (int j = 0; j < bsize; j++)
    b->SetValue(j, projectedShape[j]/sqrt(EigenValues->GetValue(j))); //!< Extract eigenvalues after normalising by s-1
//    b->SetValue(j, projectedShape[j]); //!< Extract eigenvalues after normalising by s-1
  BackProjectShape(projectedShape, shape, bsize); ///shape in stddevs by convention

  ///Replace recovered missing pixels and project again
  ///(iterative reconstruction)
  for(unsigned k = 0; k < ReconstructionIterations; k ++)
  {
    for(int j = 0; j < points->GetNumberOfPoints(); j ++)
    {
      if(weights->GetValue(j) == 0)
      {
        double point[3];
        shape->GetPoint(j, point);
        points->SetPoint(j, point); //Replace missing points with reconstructed ones
      }
    }
    shape->SetPoints(points);
    projectedShape = ProjectShape(shape, bsize); //Unweighted since missing pixels replaced
    for (int j = 0; j < bsize; j++)
      b->SetValue(j, projectedShape[j]/sqrt(EigenValues->GetValue(j))); //!< Extract eigenvalues after normalising by s-1
//      b->SetValue(j, projectedShape[j]); //!< Extract eigenvalues after normalising by s-1
    BackProjectShape(projectedShape, shape, bsize); ///shape in stddevs by convention
  }
}

void vtkWeightedPCAAnalysisFilter::GetStandardShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize)
{
  ///Use iterative reconstruction to fill in missing parts and return final parameters.
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(shape->GetPoints()); //Save current shape

  ///Project first
  vnl_vector<double> projectedShape = ProjectShape(shape, bsize, weights); ///Result of projection is a sx1 matrix
  for (int j = 0; j < bsize; j++)
    b->SetValue(j, projectedShape[j]/sqrt(EigenValues->GetValue(j))); //!< Extract eigenvalues after normalising by s-1
//    b->SetValue(j, projectedShape[j]); //!< Extract eigenvalues after normalising by s-1
  BackProjectShape(projectedShape, shape, bsize); ///shape in stddevs by convention
}

void vtkWeightedPCAAnalysisFilter::GetWeightedShapeParameters(vtkPointSet *shape, vtkFloatArray *weights, vtkFloatArray *b, int bsize)
{
//  cerr << "Weighted Shape Parameters" << endl;
  b->SetNumberOfValues(bsize);

  const size_t n = this->GetOutput(0)->GetNumberOfPoints();
  const vnl_vector<double> onesRow(b->GetNumberOfTuples(), 1.0);
  vnl_vector<double> weightsSqrt(3*n, 1.0);
  for(vtkIdType j = 0; j < weights->GetNumberOfTuples(); j ++)
  {
    weightsSqrt(3*j) = sqrt(weights->GetTuple1(j));
    weightsSqrt(3*j+1) = sqrt(weights->GetTuple1(j));
    weightsSqrt(3*j+2) = sqrt(weights->GetTuple1(j));
  }

  const vnl_matrix<double> weightsMatrix = outer_product(weightsSqrt, onesRow);
  const vnl_matrix<double> trimmedEigenVectors = EigenVectors.extract(3*n, b->GetNumberOfTuples());
  ///Compute (w_.j o U_.j)^-1
  const vnl_matrix<double> weightedEigenVectors = vnl_matrix_inverse<double>(element_product(weightsMatrix, trimmedEigenVectors)).pinverse();//Pseudo inverse, result is a sx3n matrix

  vnl_vector<double> weightedShapevec(n*3, 0.0); //Column vector
  ///Copy shape and subtract mean shape
  for(size_t i = 0; i < n; i++)
  {
    double p[3];
    shape->GetPoint(i, p);

    weightedShapevec(i*3) = (p[0] - MeanShape[i*3])*weightsSqrt(3*i);
    weightedShapevec(i*3+1) = (p[1] - MeanShape[i*3+1])*weightsSqrt(3*i+1);
    weightedShapevec(i*3+2) = (p[2] - MeanShape[i*3+2])*weightsSqrt(3*i+2);
  }

  ///Compute a
  vnl_vector<double> coefficients = weightedEigenVectors*weightedShapevec;

  for(vtkIdType j = 0; j < coefficients.size(); j ++)
    b->SetTuple1(j, coefficients(j)/sqrt(EigenValues->GetValue(j)));
}

bool vtkWeightedPCAAnalysisFilter::Load(const char * filename)
{
  cerr << "Loading Robust SSM with FileName " << filename << endl;

  typedef int headerType; //long long is 64-bit always
  typedef vtkFloatingPointType writeType; //vtkFloatingPointType is normally a double

  std::ifstream fin(filename, std::ios::binary);
  if (!fin.is_open())
  {
    cerr << " Robust SSM3D File " << filename << " does not exist " << endl;
    return false;
  }

  int headerSize = 3;
  headerType *header = new headerType[headerSize];
  fin.read((char*)header, headerSize*sizeof(headerType));

  headerType totalShapes = header[0];
  headerType totalAlignedPoints = header[0]*header[1]*header[2];
  headerType totalMeanPoints = header[1]*header[2];
  headerType totalWeightValues = header[0]*header[1]*header[2];
  headerType totalEigenvalues = totalShapes;
  headerType totalEigenvectors = totalAlignedPoints;

  headerType dataSize = totalAlignedPoints + totalMeanPoints + totalWeightValues + totalEigenvalues + totalEigenvectors;

  //std::cout << "Amount of data to read is " << dataSize << std::endl;
  cerr << "Number of shapes is " << header[0] << endl;
  cerr << "Number of Dimensions " << header[1] << endl;
  cerr << "Number of Landmarks " << header[2] << endl;

  this->SetNumberOfInputs(totalShapes);

  if ((header[0] <=0) || (header[1] <= 0) || (header[2] <= 0))
  {
    delete [] header;
    fin.close();
    cerr << "Not a valid WPCA file: " << filename << endl;
    return false;
  }

  // Check file size to ensure it's valid using header info
  headerType begin;
  begin = fin.tellg();
  fin.seekg (0, std::ios::end);

  fin.seekg(begin);
  writeType *data = new writeType[dataSize];
  fin.read((char*)data, dataSize*sizeof(writeType));
  if (!fin)
  {
    delete [] header;
    delete [] data;
    fin.close();
    cerr << "Failed reading WPCA data from file " << filename << endl;
    return false;
  }

  int alignedCount = 0;
  // For each shape get PointSet
  cerr << "Loading Shapes" << endl;
  for (int idx = 0; idx < totalShapes; idx ++)
  {
    vtkPolyData * aligned = vtkPolyData::New();
    vtkPoints *points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
    {
      // Store each point in array before adding to pointset
      writeType *value = new writeType[header[1]];
      for (int k = 0; k < header[1]; k++)
      {
        value[k] = data[alignedCount];
        alignedCount++;
      }
      points->InsertNextPoint(value);
      delete [] value;
    }
    aligned->SetPoints(points);
    points->Delete();

    this->SetNthInputConnection(0, idx, aligned ? aligned->GetProducerPort() : 0);
  }

  MeanShape.set_size(totalMeanPoints);
  int meanCount = 0;
  for (int j = 0; j < totalMeanPoints; j++)
  {
    MeanShape[j] = data[totalAlignedPoints + meanCount];
    meanCount ++;
  }

  //std::cout << "Number of points in shape " << i << " is " << shape->GetNumberOfPoints() << std::endl;

  //Get Weights, Eigenvalues and vectors.
  //Get Weights
  vnl_matrix<double> weights(static_cast<double>(totalWeightValues)/totalShapes, totalShapes, 0.0);
  double *weightsPtr = weights.data_block();
  for (int j = 0; j < totalWeightValues; j ++)
  {
    weightsPtr[j] = data[totalAlignedPoints + totalMeanPoints + j];
    //~ cout << m_Weights->GetValue(j) << ", ";
  }
  //~ cout << endl;
  this->SetWeights(weights);

  //Get Eigenvalues
  vtkFloatArray *eigenVals = vtkFloatArray::New();
  eigenVals->SetNumberOfTuples(totalEigenvalues);
  eigenVals->SetNumberOfComponents(1);
  for (int j = 0; j < eigenVals->GetNumberOfTuples(); j ++)
  {
    eigenVals->SetValue(j, data[totalAlignedPoints + totalMeanPoints + totalWeightValues + j]);
    //~ cout << eigenVals->GetValue(j) << ", ";
  }
  //~ cout << endl;
  this->SetEigenValues(eigenVals);
  cerr << "Loaded Eigenvalues: " << this->GetEigenValues()->GetNumberOfTuples() << endl;

  //Get Eigenvectors
  vnl_matrix<double> vectors(static_cast<double>(totalEigenvectors)/totalShapes, totalShapes, 0.0);
  double *vectorsPtr = vectors.data_block();
//  cerr << "Eigenvectors: " << endl;
  for (int j = 0; j < totalEigenvectors; j ++)
  {
    vectorsPtr[j] = data[totalAlignedPoints + totalMeanPoints + totalWeightValues + totalEigenvalues + j];
//    cerr << vectorsPtr[j] << ", ";
  }
//  cerr << endl;
  this->SetEigenVectors(vectors);
  cerr << "Loaded Eigenvectors: " << this->GetEigenVectors().rows() << "x" << this->GetEigenVectors().cols() << endl;

  //Set the weights in the WPCA object
  PreLoaded = true;
  Update();
  cerr << "Updated PCA" << endl;

  delete [] data;
  delete [] header;
  cerr << "Finished Reading WPCA File" << endl;
  return true;
}

bool vtkWeightedPCAAnalysisFilter::Save(const char * filename)
{
  cerr << "About to Save WPCA with filename " << filename << endl;
  /// We want to save the current model to file
  /// File format is header, aligned (point) data,
  typedef int headerType; //long long is 64-bit always
  typedef vtkFloatingPointType writeType; //vtkFloatingPointType is normally a double

  const size_t n = this->GetOutput(0)->GetNumberOfPoints();
  headerType sz = this->GetNumberOfInputConnections(0);
  if (sz > 0)
  {
    std::ofstream fout(filename, std::ios::binary);
    if (!fout.is_open())
    {
      cerr << "Cannot create WPCA File " << filename;
      return false;
    }

    int headerSize = 3;
    headerType * header = new headerType[headerSize];
    header[0] = sz;
    header[1] = 3;
    header[2] = n;

    headerType totalShapes = header[0];
    headerType totalAlignedPoints = header[0]*header[1]*header[2];
    headerType totalMeanPoints = header[1]*header[2];
    headerType totalWeightValues = header[0]*header[1]*header[2];
    headerType totalEigenvalues = totalShapes;
    headerType totalEigenvectors = totalAlignedPoints;

    //std::cout << "Header is " << header[0] << " " << header[1] << " " << header[2] << std::endl;
    headerType sizeData = totalAlignedPoints + totalMeanPoints + totalWeightValues + totalEigenvalues + totalEigenvectors;
    writeType *writer = new writeType [sizeData];
    int count = 0;

    ///Extract aligned points
    cerr << "Extracting Aligned Shapes." << endl;
    for (int i = 0; i < sz; i++)
    {
      vtkPoints * points = this->GetOutput(i)->GetPoints();
      for (int j = 0; j < points->GetNumberOfPoints(); j++)
      {
        double value[3];

        points->GetPoint(j, value);
        writer[count] = value[0];
        count++;
        writer[count] = value[1];
        count++;
        writer[count] = value[2];
        count++;
      }
    }

    cerr << "Extracting Mean Shape." << endl;
    for (size_t j = 0; j < n; j++)
    {
      writer[count] = MeanShape[3*j];
      count++;
      writer[count] = MeanShape[3*j+1];
      count++;
      writer[count] = MeanShape[3*j+2];
      count++;
    }

    ///Extract Weights
    cerr << "Extracting Weights." << endl;
    for (size_t j = 0; j < Weights.rows(); j ++)
    {
      for(size_t k = 0; k < Weights.cols(); k ++)
      {
        writer[count] = Weights(j,k);
        count++;
      }
    }

    ///Extract Eigenvalues
    cerr << "Extracting Eigenvalues." << endl;
    for (int j = 0; j < EigenValues->GetNumberOfTuples(); j++)
    {
      writer[count] = EigenValues->GetValue(j);
      count++;
    }

    ///Extract Eigenvectors
    cerr << "Extracting Eigenvectors." << endl;
    for (size_t j = 0; j < EigenVectors.rows(); j ++)
    {
      for(size_t k = 0; k < EigenVectors.cols(); k ++)
      {
        writer[count] = EigenVectors(j,k);
        count++;
      }
    }

    // Write data to file
    cerr << "Writing Data to File." << endl;
    fout.write((char *)(header),headerSize*sizeof(headerType));
    fout.write((char *)(writer),sizeData*sizeof(writeType));
    fout.close();
    delete [] header;
    delete [] writer;
    cerr << "Finished writing file " << filename << endl;
  }
  else
  {
    cerr << "There is no model to save" << endl;
    return false;
  }

  return true;
}

//----------------------------------------------------------------------------
// public
void vtkWeightedPCAAnalysisFilter::SetNumberOfInputs(int n)
{
  this->SetNumberOfInputConnections(0, n);
  this->SetNumberOfOutputPorts(n);

  // initialise the outputs
  for (int i=0; i<n; i++)
  {
    vtkPoints *points = vtkPoints::New();
    vtkPolyData *ps = vtkPolyData::New();
    ps->SetPoints(points);
    points->Delete();
    this->GetExecutive()->SetOutputData(i,ps);
    ps->Delete();
  }

  // is this the right thing to be doing here? if we don't initialise the outputs here
  // then the filter crashes but vtkPolyData may not be the type of the inputs
}

//----------------------------------------------------------------------------
void vtkWeightedPCAAnalysisFilter::SetInput(int idx, vtkPointSet *p, vtkFloatArray *weights)
{
  this->SetNthInputConnection(0, idx, p ? p->GetProducerPort() : 0);

  if (Weights.empty())
    Weights.set_size(3*p->GetNumberOfPoints(), this->GetNumberOfInputConnections(0));
//  cerr << "Weights Size: " << Weights.rows() << "x" << Weights.cols() << endl;

  for (int j = 0; j < weights->GetNumberOfTuples(); j++)
  {
    Weights(3*j, idx) = weights->GetValue(j);
    Weights(3*j+1, idx) = weights->GetValue(j);
    Weights(3*j+2, idx) = weights->GetValue(j);
  }
}

//----------------------------------------------------------------------------
void vtkWeightedPCAAnalysisFilter::SetInput(int idx, vtkDataObject* input, vtkFloatArray *weights)
{
  vtkPointSet* p = vtkPointSet::SafeDownCast(input);

  if (p)
  {
    this->SetInput(idx, p, weights);
  }
  else
  {
    vtkErrorMacro(<< this->GetClassName() << " input is a " <<
                  input->GetClassName() << " -- it should be a vtkPointSet");
  }
}

//----------------------------------------------------------------------------
vtkPointSet* vtkWeightedPCAAnalysisFilter::GetInput(int idx)
{
  return vtkPointSet::SafeDownCast(
           this->GetExecutive()->GetInputData(0, idx));
}

//----------------------------------------------------------------------------
int vtkWeightedPCAAnalysisFilter::FillInputPortInformation(int port,
    vtkInformation *info)
{
  int retval = this->Superclass::FillInputPortInformation(port, info);
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return retval;
}

//----------------------------------------------------------------------------
// public
void vtkWeightedPCAAnalysisFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  this->EigenValues->PrintSelf(os,indent.GetNextIndent());
  std::cout << "ConvergenceTolerance: " << ConvergenceTolerance << std::endl;
  std::cout << "Precision: " << Precision << std::endl;
  std::cout << "MaximumIterations: " << MaximumIterations << std::endl;
  std::cout << "ReconstructionIterations: " << ReconstructionIterations << std::endl;
  std::cout << "MissingMode: " << MissingMode << std::endl;
  std::cout << "PreLoaded: " << PreLoaded << std::endl;
  std::cout << "RandomGuess: " << RandomGuess << std::endl;
  std::cout << "ReconstructionError: " << ReconstructionError << std::endl;
}

//----------------------------------------------------------------------------
// public
int vtkWeightedPCAAnalysisFilter::GetModesRequiredFor(double proportion)
{
  int i;

  if(proportion == 1.0)
    return this->GetNumberOfInputConnections(0);

  SumOfEigenValues = 0.0;
  for (i=0; i<this->EigenValues->GetNumberOfTuples(); i++)
  {
    SumOfEigenValues += this->EigenValues->GetValue(i);
  }
  //~ cout << "Sum of Eigenvalues: " << SumOfEigenValues << endl;

  double running_total = 0.0F;
  for (i=0; i<this->EigenValues->GetNumberOfTuples(); i++)
  {
    running_total += this->EigenValues->GetValue(i)/SumOfEigenValues;
    if (running_total>=proportion)
    {
      return i+1;
    }
  }

  return EigenValues->GetNumberOfTuples();
}

vnl_vector<double> vtkWeightedPCAAnalysisFilter::ProjectShape(vtkPointSet *shape, int bsize, vtkFloatArray *weights)
{
//  cerr << "Project Shape" << endl;
  const int n = this->GetOutput(0)->GetNumberOfPoints();
  vnl_vector<double> projectedShape(bsize, 0.0);

  if (shape->GetNumberOfPoints() != n)
  {
    vtkErrorMacro(<<"Input shape does not have the correct number of points");
    vtkErrorMacro(<<"Number of points was: " << shape->GetNumberOfPoints());
    return projectedShape;
  }

  vnl_vector<double> shapevec(n*3, 0.0); //Column vector
  ///Copy shape and subtract mean shape
  if(!weights)
  {
    for (int i = 0; i < n; i++)
    {
      double p[3];
      shape->GetPoint(i, p);

      shapevec(i*3) = (p[0] - MeanShape[i*3]);
      shapevec(i*3+1) = (p[1] - MeanShape[i*3+1]);
      shapevec(i*3+2) = (p[2] - MeanShape[i*3+2]);
    }
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      double p[3];
      shape->GetPoint(i, p);

      shapevec(i*3) = (p[0] - MeanShape[i*3]) * weights->GetValue(i);
      shapevec(i*3+1) = (p[1] - MeanShape[i*3+1]) * weights->GetValue(i);
      shapevec(i*3+2) = (p[2] - MeanShape[i*3+2]) * weights->GetValue(i);
    }
  }

  for (int i = 0; i < bsize; i++) //Project only to first bsize eigenvectors.
  {
    // Project the shape onto eigenvector i
    for (int j = 0; j < 3*n; j++)
      projectedShape(i) += shapevec(j) * InverseWeightedEigenVectors(i, j);
      //In non missing mode, InverseWeightedEigenVectors.transpose() = EigenVectors
  }

//  cerr << "End Project Shape" << endl;
  return projectedShape;
}

void vtkWeightedPCAAnalysisFilter::BackProjectShape(vnl_vector<double> components, vtkPointSet *shape, const int k)
{
//  cerr << "Back-Project Shape" << endl;
  const int n = this->GetOutput(0)->GetNumberOfPoints();

  ///Weight the eigenvectors by the coefficients provided, components is a column vector. Add mean shape to get actual shape
  ///Use unweighted mean if flag set
  vnl_vector<double> shapeVector(3*n, 0.0);

  for (int j = 0; j < n*3; j++) //Back Project only to first bsize eigenvectors.
  {
    shapeVector(j) = MeanShape(j);

    for (int i = 0; i < k; i++) // Back Project the shape onto eigenvector i
      shapeVector(j) += components(i)*EigenVectors(j, i);
  }

  ///Copy shape
  for (int i = 0; i < n; i++)
    shape->GetPoints()->SetPoint(i, shapeVector(i*3), shapeVector(i*3+1), shapeVector(i*3+2));

//  cerr << "End Back-Project Shape" << endl;
}

void vtkWeightedPCAAnalysisFilter::BackProjectWeightedShape(vnl_vector<double> components, vtkPointSet *shape, vtkFloatArray *weights, const int k)
{
//  cerr << "Back-Project Shape" << endl;
  const int n = this->GetOutput(0)->GetNumberOfPoints();

  //Copy weights
  vnl_vector<double> weightsSqrt(3*n, 1.0);
  for(vtkIdType j = 0; j < weights->GetNumberOfTuples(); j ++)
  {
    weightsSqrt(3*j) = sqrt(weights->GetTuple1(j));
    weightsSqrt(3*j+1) = sqrt(weights->GetTuple1(j));
    weightsSqrt(3*j+2) = sqrt(weights->GetTuple1(j));
  }

  ///Weight the eigenvectors by the coefficients provided, components is a column vector. Add mean shape to get actual shape
  ///Use unweighted mean if flag set
  vnl_vector<double> shapeVector(3*n, 0.0);

  for (int j = 0; j < n*3; j++) //Back Project only to first bsize eigenvectors.
  {
    shapeVector(j) = MeanShape(j);

    for (int i = 0; i < k; i++) // Back Project the shape onto eigenvector i
      shapeVector(j) += weightsSqrt(j)*components(i)*EigenVectors(j, i);
  }

  ///Copy shape
  for (int i = 0; i < n; i++)
    shape->GetPoints()->SetPoint(i, shapeVector(i*3), shapeVector(i*3+1), shapeVector(i*3+2));

//  cerr << "End Back-Project Shape" << endl;
}

//=======================================================================
//: Orthogonalise a basis using modified Gram-Schmidt
// Transform basis {vk} to orthonormal basis {ek} with k in range 1..N
// for j = 1 to N
//     ej = vj
//     for k = 1 to j-1
//         ej = ej - <ej,ek>ek   //NB Classical GS has vj in inner product
//     end
//     ej = ej/|ej|
//  end
// \author Martin Roberts
// Code from VXL SVN repository, contrib/mul/mbl

//: Convert input basis {v} to orthonormal basis {e}
// Each basis vector is a column of v, and likewise the orthonormal bases are returned as columns of e
void vtkWeightedPCAAnalysisFilter::Gram_Schmidt_Orth(const vnl_matrix<double>& v, vnl_matrix<double>& e)
{
  unsigned N=v.cols();

  //Note internally it is easier to deal with the basis as a vector of vnl_vectors
  //As matrices are stored row-wise
  vcl_vector<vnl_vector<double > > vbasis(N);
  vcl_vector<vnl_vector<double > > evecs(N);

  //Copy into more convenient holding storage as vector of vectors
  //And also initialise output basis to input
  for (unsigned jcol=0; jcol<N; ++jcol)
  {
    evecs[jcol] = vbasis[jcol] = v.get_column(jcol);
  }
  evecs[0].normalize();

  for (unsigned j=1; j<N; ++j)
  {
    //Loop over previously created bases and subtract off partial projections
    //Thus producing orthogonality

    unsigned n2 = j-1;
    for (unsigned k=0; k<=n2; ++k)
    {
      evecs[j] -= dot_product(evecs[j],evecs[k]) * evecs[k];
    }
    evecs[j].normalize();
  }

  //And copy into column-wise matrix (kth base is the kth column)
  e.set_size(v.rows(),N);
  for (unsigned jcol=0; jcol<N; ++jcol)
  {
    e.set_column(jcol,evecs[jcol]);
  }
}

void vtkWeightedPCAAnalysisFilter::printProgBar( int percent )
{
  std::string bar;

  for (int i = 0; i < 50; i++)
  {
    if ( i < (percent/2))
      bar.replace(i,1,"=");
    else if ( i == (percent/2))
      bar.replace(i,1,">");
    else
      bar.replace(i,1," ");
  }

  std::cerr<< "\r" "[" << bar << "] ";
  std::cerr.width( 3 );
  std::cerr<< percent << "%     " << std::flush;

  if (percent == 100)
    std::cerr << endl;
}
