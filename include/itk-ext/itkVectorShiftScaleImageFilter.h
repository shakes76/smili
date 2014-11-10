/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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
#ifndef __itkVectorShiftScaleImageFilter_h
#define __itkVectorShiftScaleImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

///\todo This file has been extended from ITK to handle scaling of VectorImages. Replace this with ITK equivalent when their bugs are fixed.

namespace itk
{
// This functor class applies a scaling transformation A.x
// to input values assuming a variable length vector is used as TInput.
namespace Functor
{
template< typename TInput, typename  TOutput >
class VectorMagnitudeScaleTransform
{
public:
  typedef typename NumericTraits< typename TInput::ValueType >::RealType RealType;
  VectorMagnitudeScaleTransform() {}
  ~VectorMagnitudeScaleTransform() {}
  void SetFactor(RealType a) { m_Factor = a; }
  void SetDimension(unsigned int n) { m_Dimension = n; }
  bool operator!=(const VectorMagnitudeScaleTransform & other) const
  {
    if ( m_Factor != other.m_Factor )
      {
      return true;
      }
    return false;
  }

  bool operator==(const VectorMagnitudeScaleTransform & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput & x) const
  {
    TOutput result;

    result.SetSize(m_Dimension);
    for ( unsigned int i = 0; i < m_Dimension; i++ )
      {
      const RealType scaledComponent = static_cast< RealType >( x[i] ) * m_Factor;
      result[i] = static_cast< typename TOutput::ValueType >( scaledComponent );
      }
    return result;
  }

private:
  RealType m_Factor;
  unsigned int m_Dimension;
};
}  // end namespace functor

/** \class VectorShiftScaleImageFilter
 * \brief Applies a linear transformation to the magnitude of pixel vectors in a
 * VectorImage.
 *
 * VectorShiftScaleImageFilter applies pixel-wise a linear transformation
 * to the intensity values of input image pixels. The linear transformation is
 * defined by the user in terms of scaling value of the vectors.
 *
 * All computations are performed in the precison of the input pixel's
 * RealType. Before assigning the computed value to the output pixel.
 *
 * \sa ShiftScaleImageFilter
 *
 * \ingroup IntensityImageFilters  MultiThreaded
 *
 * \ingroup ITKImageIntensity
 *
 * \wiki
 * \wikiexample{Images/VectorShiftScaleImageFilter,Apply a transformation to the magnitude of vector valued image pixels}
 * \endwiki
 */
template< typename  TInputImage, typename  TOutputImage = TInputImage >
class ITK_EXPORT VectorShiftScaleImageFilter:
  public
  UnaryFunctorImageFilter< TInputImage, TOutputImage,
                           Functor::VectorMagnitudeScaleTransform<
                             typename TInputImage::PixelType,
                             typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef VectorShiftScaleImageFilter Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    Functor::VectorMagnitudeScaleTransform<
      typename TInputImage::PixelType,
      typename TOutputImage::PixelType > >    Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef typename TOutputImage::PixelType                    OutputPixelType;
  typedef typename TInputImage::PixelType                     InputPixelType;
  typedef typename InputPixelType::ValueType                  InputValueType;
  typedef typename OutputPixelType::ValueType                 OutputValueType;
  typedef typename NumericTraits< InputValueType >::RealType  InputRealType;
  typedef typename NumericTraits< OutputValueType >::RealType OutputRealType;

  typedef typename Superclass::InputImageType    InputImageType;
  typedef typename Superclass::InputImagePointer InputImagePointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(VectorShiftScaleImageFilter, UnaryFunctorImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get the Scale used for the linear transformation
      of magnitude values.*/
  itkGetConstReferenceMacro(Scale, InputRealType);
  itkSetMacro(Scale, InputRealType);

  /** Set/Get the Scale used for the linear transformation
      of magnitude values.*/
  itkGetConstReferenceMacro(Dimension, unsigned int);
  itkSetMacro(Dimension, unsigned int);

  /** Process to execute before entering the multithreaded section */
  void BeforeThreadedGenerateData(void);

  /** Print internal ivars */
  void PrintSelf(std::ostream & os, Indent indent) const;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< InputValueType > ) );
  itkConceptMacro( OutputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< OutputValueType > ) );
  /** End concept checking */
#endif

protected:
  VectorShiftScaleImageFilter();
  virtual ~VectorShiftScaleImageFilter() {}

private:
  VectorShiftScaleImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                    //purposely not implemented

  InputRealType m_Scale;
  unsigned int m_Dimension;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorShiftScaleImageFilter.hxx"
#endif

#endif
