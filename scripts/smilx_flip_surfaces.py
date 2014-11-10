#!/usr/bin/dgv2
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
# Load the joints and mydata module into DGV first
execfile("filenames.py")

bone = "fem_L"

affineMeshPath = "surfaces_L/"
outputPath = "surfaces_L_to_R/"

outputExt = ".vtk"
outputPrefix = "flip_" + bone + "_"

affineList, caseList = getSortedFileListAndCases(affineMeshPath, 0, '*.vtk', True)

for affine, case in zip(affineList, caseList):
    MainWindow.loadFile(affine) # load mean mesh
    
    #register meshes
    currentModel = MainWindow.activeModel() # get the loaded mesh
    currentModel.flip(True, False, False) # flip along x
    
    outputName = outputPrefix + str(case) + outputExt
    
    milxQtFile.saveModel(outputPath+outputName, currentModel) # save result
    
    currentModel.setDeletableOnClose(True) # delete on close
    currentModel.close() # close (since its delete on close)
print "Processing Done"
