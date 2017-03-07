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
# Load the joints and mydata module into smilx first
execfile("filenames.py")

meshPath ='segmentations/'
manualMeshPath = 'meshes_prostate_originals/'
outputPath = 'hausdoff_original_meshes/'
outputExt = ".vtk"
outputPrefix = "Hausdoff_"
prepend_full_path = True

manualList, caseList = getSortedFileListAndCases(manualMeshPath, 0, '*.vtk', prepend_full_path)
meshList = getSortedFileListGivenCases(meshPath, caseList, '*.vtk', prepend_full_path)

for manual, mesh, case in zip(manualList, meshList, caseList):
    #load mesh
    MainWindow.loadFile(mesh) # load mean mesh
    currentMesh = MainWindow.activeModel() # get the loaded mesh
    MainWindow.unify() #add active window to unified window
    MainWindow.loadFile(manual) # load mean mesh
    currentManual = MainWindow.activeModel() # get the loaded mesh
    MainWindow.unify() #add active window to unified window
    
    unifiedWindow = MainWindow.getUnifiedWindow() # get the unified window
    unifiedWindow.generateDifference(1e3)
    
    currentResult = MainWindow.activeModel() # get the loaded mesh
    
    outputName = outputPrefix + str(case) + outputExt
    milxQtFile.saveModel(outputPath+outputName, currentResult) # save result
    
    MainWindow.closeTabAllWindows()
print "Processing Done"