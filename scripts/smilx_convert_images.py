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
This script opens a model one at a time from the path provided and processes it. The processing is 
Clean, Smooth, Decimate by 25%, Smooth, Decimate by 25%, Smooth.
The resulting model is then saved to disk with naming convention provided.
'''
import os

# Load the joints and mydata module into smilx first
execfile("filenames.py")

parent = '../'
outputExt = "vti"
imagePath = parent+'t2maps/'
outputPath = parent+'t2maps_'+outputExt+'/'
outputPrefix = ""
prepend_full_path = True

imageList, caseList = getSortedFileListAndCases(imagePath, 0, '*.vtk', prepend_full_path)

for image, case in zip(imageList, caseList):
    filename, fileExt = os.path.splitext(image)
    if fileExt == 'gz':
        filename, fileExt = os.path.splitext(image)
    
    #load mesh
#    MainWindow.loadFile(image) # load mean mesh
#    currentMesh = MainWindow.activeModel() # get the loaded meshs
    
    outputName = outputPrefix + filename + "." + outputExt
#    milxQtFile.saveImage(outputPath+outputName, currentResult) # save result
    
    MainWindow.closeTabAllWindows()
print "Processing Done"