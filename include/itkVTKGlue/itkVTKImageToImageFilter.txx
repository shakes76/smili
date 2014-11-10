/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVTKImageToImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/09/06 20:58:41 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkVTKImageToImageFilter_txx
#define _itkVTKImageToImageFilter_txx

#include "itkVTKImageToImageFilter.h"

namespace itk
{



/**
 * Constructor
 */
template <class TOutputImage>
VTKImageToImageFilter<TOutputImage>
::VTKImageToImageFilter()
{

  m_Exporter = vtkImageExport::New();

  this->SetUpdateInformationCallback( m_Exporter->GetUpdateInformationCallback());
  this->SetPipelineModifiedCallback( m_Exporter->GetPipelineModifiedCallback());
  this->SetWholeExtentCallback( m_Exporter->GetWholeExtentCallback());
  this->SetSpacingCallback( m_Exporter->GetSpacingCallback());
  this->SetOriginCallback( m_Exporter->GetOriginCallback());
  this->SetScalarTypeCallback( m_Exporter->GetScalarTypeCallback());
  this->SetNumberOfComponentsCallback( m_Exporter->GetNumberOfComponentsCallback());
  this->SetPropagateUpdateExtentCallback( m_Exporter->GetPropagateUpdateExtentCallback());
  this->SetUpdateDataCallback( m_Exporter->GetUpdateDataCallback());
  this->SetDataExtentCallback( m_Exporter->GetDataExtentCallback());
  this->SetBufferPointerCallback( m_Exporter->GetBufferPointerCallback());
  this->SetCallbackUserData( m_Exporter->GetCallbackUserData());

}




/**
 * Destructor
 */
template <class TOutputImage>
VTKImageToImageFilter<TOutputImage>
::~VTKImageToImageFilter()
{
  if( m_Exporter )
    {
    m_Exporter->Delete();
    m_Exporter = 0;
    }
}



/**
 * Set a vtkImageData as input 
 */
template <class TOutputImage>
void
VTKImageToImageFilter<TOutputImage>
::SetInput( vtkImageData * inputImage )
{
#if VTK_MAJOR_VERSION <= 5
  m_Exporter->SetInput( inputImage );
#else
  m_Exporter->SetInputData( inputImage );
#endif
}



/**
 * Get the exporter filter
 */
template <class TOutputImage>
vtkImageExport *
VTKImageToImageFilter<TOutputImage>
::GetExporter() const
{
  return m_Exporter;
}



/**
 * Get the importer filter
 */
template <class TOutputImage>
const typename VTKImageToImageFilter<TOutputImage>::Superclass *
VTKImageToImageFilter<TOutputImage>
::GetImporter() const
{
  return this;
}




} // end namespace itk

#endif

