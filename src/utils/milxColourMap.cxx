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
#include "milxColourMap.h"

namespace milx
{


ColourMap::ColourMap()
{
  reset();

  mapRange[0] = 0;
  mapRange[1] = 255;

  lookupTable = vtkSmartPointer<vtkLookupTable>::New();
}

ColourMap::~ColourMap()
{

}

void ColourMap::reset()
{
  mapFlag = RAINBOW;
  generated = false;
}

vtkSmartPointer<vtkLookupTable> ColourMap::GetOutput()
{
  if(!generated)
    GenerateData();

  return lookupTable;
}

void ColourMap::GenerateData()
{
  double tblRange[2] = {0, 255};

  lookupTable->SetNanColor(0.85, 0.85, 0.85, 1.0);
  lookupTable->SetTableRange(mapRange[0], mapRange[1]);
//  lookupTable->SetTableRange(tblRange[0], tblRange[1]);

  switch(mapFlag)
  {
    case JET:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(JET_ARRAY[i][j]*tblRange[1]);
//        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case RAINBOW:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(RAINBOW_ARRAY[i][j]*tblRange[1]);
//        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case VTK:
    {
      lookupTable->SetNumberOfColors(256);
      lookupTable->SetHueRange(0.667, 0.0); // 0.667 for blue, and 0 for red
      lookupTable->Build();

      break;
    }

    case NIH:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      unsigned char *p;
      for( unsigned int i=0; i<256; i++ )
      {
        p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(NIH_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case NIH_FIRE:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      unsigned char *p;
      for( unsigned int i=0; i<256; i++ )
      {
        //unsigned char *p;
        //if((alphaStart) && (i < 100))
        //  p = lookupTable->WritePointer(i, 0);
        //else
        p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(NIH_FIRE_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case HOT:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(HOT_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      /*lookupTable->SetRange( 0.0, 255.0 );
      lookupTable->SetHueRange( 0.0, 0.1 );
      lookupTable->SetValueRange( 0.4, 0.8 );
      lookupTable->Build();*/

      break;
    }

    case COOL:

      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(COOL_ARRAY[i][j]*tblRange[1]);
//        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      /*lookupTable->SetRange( 0.0, 255.0 );
      lookupTable->SetHueRange( 0.67, 0.68 );
      lookupTable->SetValueRange( 0.4, 0.8 );
      lookupTable->Build();*/

      break;

    case COOLWARM:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(COOLWARM_ARRAY[i][j]*tblRange[1]);
//        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case KNEE:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      unsigned char *p;
      for( unsigned int i=0; i<256; i++ )
      {
        p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(KNEE_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case AAL:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      unsigned char *p;
      for( unsigned int i=0; i<256; i++ )
      {
        p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(AAL_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

    case FS:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      unsigned char *p;
      for( unsigned int i=0; i<256; i++ )
      {
        p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(FS_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;
    }

//  else if (colourmap == "Muscles")
//  {
//    lookupTable->SetNumberOfTableValues( 256 );
//    lookupTable->SetTableRange(0, tblRange[1]);
//    lookupTable->Build();
//
//    unsigned char *p = lookupTable->WritePointer(0, 1);
//    for (unsigned int j=0; j<3; ++j)
//      *p++ = static_cast<unsigned char>(Muscles[0][j]*tblRange[1]);
//    *p++ = static_cast<unsigned char>(0);
//    for( unsigned int i=1; i<256; i++ )
//    {
//      p = lookupTable->WritePointer(i, 1);
//      for (unsigned int j=0; j<3; ++j)
//        *p++ = static_cast<unsigned char>(Muscles[i][j]*tblRange[1]);
//      *p++ = static_cast<unsigned char>(tblRange[1]);
//    }
//  }
    case LOG_GRAY:

      lookupTable->SetHueRange(0.0, 0.0);
      lookupTable->SetSaturationRange(0, 0);
      lookupTable->SetNumberOfTableValues(256);
      lookupTable->SetValueRange(0, 1);
      lookupTable->SetScaleToLog10();
      lookupTable->Build();

      break;

    case GRAY:

      lookupTable->SetHueRange(0.0, 0.0);
      lookupTable->SetSaturationRange(0, 0);
      lookupTable->SetNumberOfTableValues(256);
      lookupTable->SetValueRange(0.1, 1);
      lookupTable->Build();

      break;

    case BONE:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(BONE_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;

    case SPECTRAL:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(SPECTRAL_ARRAY[i][j]*tblRange[1]);
//        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;

    case GNUPLOT:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(GNUPLOT_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;

    case CUBEHELIX:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      for( unsigned int i=0; i<256; i++ )
      {
        unsigned char *p = lookupTable->WritePointer(i, 1);
        for (unsigned int j=0; j<3; ++j)
          *p++ = static_cast<unsigned char>(CUBEHELIX_ARRAY[i][j]*tblRange[1]);
        *p++ = static_cast<unsigned char>(tblRange[1]);
      }

      break;

    default:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      break;

  }

  generated = true;
}

} //end namespace
