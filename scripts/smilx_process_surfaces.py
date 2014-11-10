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
import os, re
 
meshPath = "/home/cha588/data/parameterisations/pluim_surf/"
 
outputPath = "/home/cha588/data/parameterisations/pluim_surf_ldc/"
outputExt = ".vtk"
outputPrefix = "prostate_pluim_"
 
dirs = os.listdir(meshPath)
for file in dirs:
    MainWindow.loadFile(meshPath + file) # load mean mesh
     
    currentModel = MainWindow.activeModel() # get the loaded mesh
     
    #find case id in filename and remember it
    case_id_str = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", file)
    case = int( case_id_str[0] )
     
    outputName = outputPrefix + str(case) + outputExt
 
    #process model
    currentModel.clean()
    currentModel.smoothSinc(20) #windowed sinc
    currentModel.decimate(0.25)
    currentModel.smoothSinc(20) #windowed sinc
    currentModel.decimate(0.25)
    currentModel.smoothSinc(20) #windowed sinc
     
    milxQtFile.saveModel(outputPath+outputName, currentModel) # save result
     
    currentModel.close() # close (since its delete on close)
     
print "Processing Done"
