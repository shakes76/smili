/*=========================================================================
  Program: MILX 
  Module: milxWin32Header.h
  Author: David Hellier
  Modified by: 
  Language: C++
  Created: Fri 07 Oct 2005 15:29:37 EST

  Based on itkWin32Header.h

  Copyright: (c) 2005 CSIRO, Australia, 2005.

  This software is protected by international copyright laws. 
  Any unauthorised copying, distribution or reverse engineering is prohibited. 

  Licence: 
  All rights in this Software are reserved to CSIRO. You are only permitted 
  to have this Software in your possession and to make use of it if you have 
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef __milxWIN32Header_h
#define __milxWIN32Header_h

#include "milxConfigure.h"

// add in the Windows variants

#if defined(__CYGWIN__)
#ifndef WIN32
#define WIN32 1
#endif
#ifndef _WIN32
#define _WIN32 1
#endif
#endif

/** Disable some common warnings in MS VC++ */
#if defined(_MSC_VER)

// 'conversion' conversion from 'type1' to 'type2', possible loss of data
#pragma warning ( disable : 4244 )

// 'identifier' : truncation from 'type1' to 'type2'
#pragma warning ( disable : 4305 )

// 'conversion' : truncation of constant value
#pragma warning ( disable : 4309 )

// decorated name length exceeded, name was truncated
#pragma warning ( disable : 4503 )

// 'identifier' : identifier was truncated to 'number' characters in the
// debug information
#pragma warning ( disable : 4786 )

// 'type' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning ( disable : 4800 )

// 'identifier' : class 'type' needs to have dll-interface to be used by
// clients of class 'type2'
#pragma warning ( disable : 4251 )


// non dll-interface class 'type' used as base for dll-interface class 'type2'
#pragma warning ( disable : 4275 )

// C++ exception specification ignored except to indicate a 
// function is not __declspec(nothrow)
#pragma warning ( disable : 4290 )

// 'type' : inconsistent dll linkage.  dllexport assumed.
#pragma warning ( disable : 4273 )

// warning C4521: 'SNRdata' : multiple copy constructors specified
#pragma warning ( disable : 4521 )

// typename keyword in default template arguments is not accepted by
// MSVC.  This macro should only be used in such places.
# if !defined(CABLE_CONFIGURATION) && (_MSC_VER < 1310)
#  define MILX_TYPENAME
# else
#  define MILX_TYPENAME typename
# endif
#else
# define MILX_TYPENAME typename
#endif


#if (defined(_WIN32) || defined(WIN32)) && !defined(MILXSTATIC) 
#define MILX_EXPORT __declspec(dllexport)
#else
/* unix needs nothing */
#define MILX_EXPORT 
#endif

// Already defined mode_t in itk so don't
// redefine in WXWidgets otherwise get the following
// error:
// error C2371: 'mode_t' : redefinition; different basic types
// wxWidgets 2.9.1. Added here as its included from milObject.h
#define wxNEEDS_MODE_T 0

#endif
