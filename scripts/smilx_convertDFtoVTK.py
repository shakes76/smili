#!/usr/bin/smilx
'''=========================================================================
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
========================================================================='''
'''
This script reads a number of 2D deformation fields converts these to vector fields and outputs them as vtk files.
To use it, run smilx, in python console type "execfile('smilx_convertDFtoVTK.py')"
'''
import os, re, gc

#remove extra chars from image sequence number 
def stripID(caseString):
    newCaseString = caseString.rstrip(".")
    if int(newCaseString) != 0:
        newCaseString = newCaseString.lstrip("0")
    return newCaseString

dfPath = "/export/pelvis-data3/PeterMac_CineMRI/PeterMac_April2013/ResultsFast/JournalPaperFigures/May_C/"
outputPath = "/export/pelvis-data3/PeterMac_CineMRI/PeterMac_April2013/ResultsFast/JournalPaperFigures/May_C_VTK"
count=0
dirs = os.listdir(dfPath)
for file in dirs:
    print file
    case_id_str = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", file)
    print(case_id_str)
    case = int( stripID(case_id_str[0] ))
    outputName = "df_" + str(case) + ".vtk" 
    print(outputPath+outputName) 

    if os.path.isfile(outputPath+outputName):
      print("found")
    else:
      imageName = dfPath + file
      MainWindow.loadFile(imageName) # load image
      defImage = MainWindow.activeImage() # get pointer to image window
      defImage.vectorField(8, 3.0) #subsampling and scaling given
      defImage.hide()
      defField = MainWindow.activeModel() # get pointer to image window
      milxQtFile.saveModel(outputPath+outputName, defField) # save result
      defField.hide()
      defField.close()     
      defImage.close() # close (since its delete on close
      gc.collect() #garbage collection 

