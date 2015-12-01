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
#ifndef __MILXCOLOURMAP_H
#define __MILXCOLOURMAP_H

#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>

#include "milxGlobal.h"
#include "milxLUT.h"

namespace milx
{

/**
  \class ColourMap
  \brief Represents a the various colour maps available for VTK scalars etc. Default: NIH

  Returns the set colour map as a vtkLookupTable.
  Usage:
  \code
  milxColourMap colourMap = new milxColourMap;
      colourMap->toNIH();

  vtkLookupTable *olut = colourMap->GetOutput();

  vtkImageMapToColors *filterColorsOverlay = vtkImageMapToColors::New();
      filterColorsOverlay->SetLookupTable(olut);

      filterColorsOverlay->SetInput(filterExportVTKOverlay->GetOutput());
      filterColorsOverlay->PassAlphaToOutputOn ();
  \endcode
*/
class SMILI_EXPORT ColourMap
{
public:
  ColourMap();
  ~ColourMap();

  void reset();

  enum ColourMapFlags {JET, RAINBOW, VTK, GRAY, SEISMIC, LOG_GRAY, NIH, NIH_FIRE, AAL, FS, HOT, COOL, COOLWARM, KNEE, BONE, SPECTRAL, GNUPLOT, CUBEHELIX, HSV};

  // Set map as Rainbow
  inline void toJet()
  { reset();  mapFlag = JET; }
  inline void SetJet()
  { toJet();  }
  // Set map as Rainbow
  inline void toRainbow()
  { reset();  mapFlag = RAINBOW; }
  inline void SetRainbow()
  { toRainbow();  }
  // Set map as Inverse Rainbow
  inline void toVTK()
  { reset();  mapFlag = VTK; }
  inline void SetVTK()
  { toVTK();  }
  // Set map as Gray
  inline void toGray()
  { reset();  mapFlag = GRAY; }
  inline void SetGray()
  { toGray();  }
  // Set map as Seismic
  inline void toSeismic()
  { reset();  mapFlag = SEISMIC; }
  inline void SetSeismic()
  { toSeismic();  }
  // Set map as Log Gray
  inline void toLogGray()
  { reset();  mapFlag = LOG_GRAY; }
  inline void SetLogGray()
  { toLogGray();  }
  // Set map as NIH
  inline void toNIH()
  { reset();  mapFlag = NIH; }
  inline void SetNIH()
  { toNIH();  }
  //Set Map as NIH FIRE
  inline void toNIH_FIRE()
  { reset();  mapFlag = NIH_FIRE; }
  inline void SetNIH_FIRE()
  { toNIH_FIRE();  }
  //Set Map as AAL
  inline void toAAL()
  { reset();  mapFlag = AAL; }
  inline void SetAAL()
  { toAAL();  }
  //Set Map as FS (free surfer)
  inline void toFS()
  { reset();  mapFlag = FS; }
  inline void SetFS()
  { toFS();  }
  //Set Map as HOT
  inline void toHOT()
  { reset();  mapFlag = HOT; }
  inline void SetHOT()
  { toHOT();  }
  //Set Map as HOT
  inline void toCOOL()
  { reset();  mapFlag = COOL; }
  inline void SetCOOL()
  { toCOOL();  }
  //Set Map as HOT
  inline void toCOOLWARM()
  { reset();  mapFlag = COOLWARM; }
  inline void SetCOOLWARM()
  { toCOOLWARM();  }
  //Set Map as HOT
  inline void toKnee()
  { reset();  mapFlag = KNEE; }
  inline void SetKnee()
  { toKnee();  }
  //Set Map as Bone
  inline void toBone()
  { reset();  mapFlag = BONE; }
  inline void SetBone()
  { toBone();  }
  //Set Map as Spectral
  inline void toSpectral()
  { reset();  mapFlag = SPECTRAL; }
  inline void SetSpectral()
  { toSpectral();  }
  //Set Map as GNUPlot
  inline void toGNUPlot()
  { reset();  mapFlag = GNUPLOT; }
  inline void SetGNUPlot()
  { toGNUPlot();  }
  //Set Map as CubeHelix
  inline void toCubeHelix()
  { reset();  mapFlag = CUBEHELIX; }
  inline void SetCubeHelix()
  { toCubeHelix();  }
  //Set Map as HSV
  inline void toHSV()
  { reset();  mapFlag = HSV; }
  inline void SetHSV()
  { toHSV();  }

  inline void SetRange(double range[2])
  {
    mapRange[0] = range[0];
    mapRange[1] = range[1];
  }
  inline double* GetRange()
  { return &mapRange[0];  }

  /**
      \brief Get the colourmap set
  */
  vtkSmartPointer<vtkLookupTable> GetOutput();

protected:
  //Flags
  ColourMapFlags mapFlag;
  bool generated;

  void GenerateData();

  coordinateType mapRange[2];

  //Lookup table for map
  vtkSmartPointer<vtkLookupTable> lookupTable;
};

} //end namespace

#endif //__MILXCOLOURMAP_H
