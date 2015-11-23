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

      lookupTable->SetTableValue(0, JET_ARRAY[0][0], JET_ARRAY[0][1], JET_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, JET_ARRAY[i][0], JET_ARRAY[i][1], JET_ARRAY[i][2], 1.0);

      break;
    }

    case RAINBOW:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, RAINBOW_ARRAY[0][0], RAINBOW_ARRAY[0][1], RAINBOW_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, RAINBOW_ARRAY[i][0], RAINBOW_ARRAY[i][1], RAINBOW_ARRAY[i][2], 1.0);

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

      lookupTable->SetTableValue(0, NIH_ARRAY[0][0], NIH_ARRAY[0][1], NIH_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, NIH_ARRAY[i][0], NIH_ARRAY[i][1], NIH_ARRAY[i][2], 1.0);

      break;
    }

    case NIH_FIRE:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, NIH_FIRE_ARRAY[0][0], NIH_FIRE_ARRAY[0][1], NIH_FIRE_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, NIH_FIRE_ARRAY[i][0], NIH_FIRE_ARRAY[i][1], NIH_FIRE_ARRAY[i][2], 1.0);

      break;
    }

    case HOT:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, HOT_ARRAY[0][0], HOT_ARRAY[0][1], HOT_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, HOT_ARRAY[i][0], HOT_ARRAY[i][1], HOT_ARRAY[i][2], 1.0);

      /*lookupTable->SetRange( 0.0, 255.0 );
      lookupTable->SetHueRange( 0.0, 0.1 );
      lookupTable->SetValueRange( 0.4, 0.8 );
      lookupTable->Build();*/

      break;
    }

    case COOL:

      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, COOL_ARRAY[0][0], COOL_ARRAY[0][1], COOL_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, COOL_ARRAY[i][0], COOL_ARRAY[i][1], COOL_ARRAY[i][2], 1.0);

      break;

    case COOLWARM:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, COOLWARM_ARRAY[0][0], COOLWARM_ARRAY[0][1], COOLWARM_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, COOLWARM_ARRAY[i][0], COOLWARM_ARRAY[i][1], COOLWARM_ARRAY[i][2], 1.0);

      break;
    }

    case KNEE:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, KNEE_ARRAY[0][0], KNEE_ARRAY[0][1], KNEE_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, KNEE_ARRAY[i][0], KNEE_ARRAY[i][1], KNEE_ARRAY[i][2], 1.0);

      break;
    }

    case AAL:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, AAL_ARRAY[0][0], AAL_ARRAY[0][1], AAL_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, AAL_ARRAY[i][0], AAL_ARRAY[i][1], AAL_ARRAY[i][2], 1.0);

      break;
    }

    case FS:
    {
      lookupTable->SetNumberOfTableValues( 256 );
      lookupTable->Build();

      lookupTable->SetTableValue(0, FS_ARRAY[0][0], FS_ARRAY[0][1], FS_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, FS_ARRAY[i][0], FS_ARRAY[i][1], FS_ARRAY[i][2], 1.0);

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

      lookupTable->SetTableValue(0, 0.0, 0.0, 0.0, 0.0); //transparent

      break;

    case GRAY:

      lookupTable->SetHueRange(0.0, 0.0);
      lookupTable->SetSaturationRange(0, 0);
      lookupTable->SetNumberOfTableValues(256);
      lookupTable->SetValueRange(0.2, 1);
      lookupTable->Build();

      lookupTable->SetTableValue(0, 0.0, 0.0, 0.0, 0.0); //transparent

      break;

    case SEISMIC:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      lookupTable->SetTableValue(0, SEISMIC_ARRAY[0][0], SEISMIC_ARRAY[0][1], SEISMIC_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, SEISMIC_ARRAY[i][0], SEISMIC_ARRAY[i][1], SEISMIC_ARRAY[i][2], 1.0);

      break;

    case BONE:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      lookupTable->SetTableValue(0, BONE_ARRAY[0][0], BONE_ARRAY[0][1], BONE_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, BONE_ARRAY[i][0], BONE_ARRAY[i][1], BONE_ARRAY[i][2], 1.0);

      break;

    case SPECTRAL:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      lookupTable->SetTableValue(0, SPECTRAL_ARRAY[0][0], SPECTRAL_ARRAY[0][1], SPECTRAL_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, SPECTRAL_ARRAY[i][0], SPECTRAL_ARRAY[i][1], SPECTRAL_ARRAY[i][2], 1.0);

      break;

    case GNUPLOT:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      lookupTable->SetTableValue(0, GNUPLOT_ARRAY[0][0], GNUPLOT_ARRAY[0][1], GNUPLOT_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, GNUPLOT_ARRAY[i][0], GNUPLOT_ARRAY[i][1], GNUPLOT_ARRAY[i][2], 1.0);

      break;

    case CUBEHELIX:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      lookupTable->SetTableValue(0, CUBEHELIX_ARRAY[0][0], CUBEHELIX_ARRAY[0][1], CUBEHELIX_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
        lookupTable->SetTableValue(i, CUBEHELIX_ARRAY[i][0], CUBEHELIX_ARRAY[i][1], CUBEHELIX_ARRAY[i][2], 1.0);

      break;

    case HSV:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      lookupTable->SetTableValue(0, HSV_ARRAY[0][0], HSV_ARRAY[0][1], HSV_ARRAY[0][2], 0.0); //transparent
      for( unsigned int i=1; i<256; i++ )
          lookupTable->SetTableValue(i, HSV_ARRAY[i][0], HSV_ARRAY[i][1], HSV_ARRAY[i][2], 1.0);

      break;

    default:

      lookupTable->SetNumberOfTableValues(256);
      lookupTable->Build();

      break;

  }

  generated = true;
}

} //end namespace
