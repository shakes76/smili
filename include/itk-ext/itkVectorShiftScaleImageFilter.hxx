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
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
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
#ifndef __itkVectorShiftScaleImageFilter_hxx
#define __itkVectorShiftScaleImageFilter_hxx

#include "itkVectorShiftScaleImageFilter.h"

namespace itk
{
/**
 *
 */
template< class TInputImage, class TOutputImage >
VectorShiftScaleImageFilter< TInputImage, TOutputImage >
::VectorShiftScaleImageFilter()
{
  m_Scale = 1.0;
  m_Dimension = 3;
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
VectorShiftScaleImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Scale : "
     << static_cast< typename NumericTraits< InputRealType >::PrintType >( m_Scale )
     << std::endl;
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
VectorShiftScaleImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  // set up the functor values
  this->GetFunctor().SetFactor(m_Scale);
  this->GetFunctor().SetDimension(m_Dimension);
//  std::cout << "Dimension for scaling vector image: " << m_Dimension << endl;
}
} // end namespace itk

#endif
